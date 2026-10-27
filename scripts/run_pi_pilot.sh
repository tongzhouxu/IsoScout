#!/usr/bin/env bash
# Run from the prepared pi-benchmark directory with its Python venv activated.
set -euo pipefail
org=${1:-cronobacter}
case "$org" in cronobacter|listeria|campylobacter|salmonella|ecoli_shigella) ;; *) echo 'Unknown bundle name' >&2; exit 2 ;; esac
for tool in isoscout mashpit minimap2 k8 paftools.js datasets; do command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }; done
bundle="$(pwd)/bundles/$org"
test -f "$bundle/inputs.json"
readarray -t selected < <(python - "$bundle/inputs.json" <<'PY'
import json, sys
q = json.load(open(sys.argv[1]))['queries'][0]
print(q['path'])
print(q['query_id'])
PY
)
test "${#selected[@]}" -eq 2
test -f "$bundle/${selected[0]}"
mkdir -p cli-databases results
if [[ -e "cli-databases/$org" && ! -L "cli-databases/$org" ]]; then
    echo "Expected a symlink, found an existing path: cli-databases/$org" >&2; exit 1
fi
ln -sfn "$bundle/database" "cli-databases/$org"
run=$(mktemp -d "$(pwd)/results/${org}-isoscout-full.XXXXXXXX")
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if command -v vcgencmd >/dev/null; then
    { vcgencmd measure_temp; vcgencmd get_throttled; } > "$run/pi-before.txt"
fi
{ isoscout --version; minimap2 --version; datasets version; python -m pip freeze; } > "$run/versions.txt" 2>&1
status=0
/usr/bin/time -v -o "$run/time.txt" \
    isoscout screen --executor native "$bundle/${selected[0]}" \
    --organism "$org" --database-root "$(pwd)/cli-databases" \
    --query-accession "${selected[1]}" --comparison-workers 1 \
    --output "$run/analysis" 2>&1 | tee "$run/run.log" || status=$?
printf '%s\n' "$status" > "$run/exit-code.txt"
if command -v vcgencmd >/dev/null; then
    { vcgencmd measure_temp; vcgencmd get_throttled; } > "$run/pi-after.txt"
fi
printf '\nResults: %s\nExit code: %s\n' "$run" "$status"
if [[ -f "$run/analysis/result.json" ]]; then
    python - "$run/analysis/result.json" <<'PY' | tee "$run/summary.txt"
import json, sys
r = json.load(open(sys.argv[1])); s = r.get('snp_resolution', {})
for name, value in [('analysis_mode', r.get('analysis_mode')), ('status', r.get('status')),
                    ('SNP backend', s.get('backend')), ('SNP status', s.get('status')),
                    ('stop_reason', r.get('stop_reason'))]:
    print(f'{name}: {value}')
print('Stage wall times (seconds; downloads are separate):')
for name, row in r.get('stage_timings', {}).items():
    seconds = row['elapsed_seconds']
    value = 'not run' if seconds is None else f'{seconds:.3f}'
    print(f"  {name}: {value}; calls={row['calls']}; failed_calls={row['failed_calls']}")
PY
fi
exit "$status"
