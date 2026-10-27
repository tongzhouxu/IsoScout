#!/usr/bin/env python3
"""Serial full-IsoScout benchmark over a prepared bundle; never overwrites runs."""
import argparse
import csv
import hashlib
import json
import random
import shutil
import subprocess
from pathlib import Path


def read_metrics(path):
    values = {}
    for line in path.read_text().splitlines():
        if 'Elapsed (wall clock)' in line:
            parts = line.rsplit('): ', 1)[1].strip().split(':')
            seconds = 0.
            for part in parts:
                seconds = seconds * 60 + float(part)
            values['total_seconds'] = seconds
        elif 'Maximum resident set size (kbytes):' in line:
            values['peak_rss_kib'] = int(line.rsplit(':', 1)[1])
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--organism', required=True, choices=['cronobacter','listeria','campylobacter','salmonella','ecoli_shigella'])
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--seed', type=int, default=20261027)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    manifest = Path('bundles') / args.organism / 'inputs.json'
    inputs = json.loads(manifest.read_text())
    queries = inputs['queries']
    if not queries:
        parser.error('Bundle contains no queries')
    script = Path(__file__).resolve().with_name('run_pi_pilot.sh')
    if not script.is_file():
        parser.error('Keep run_pi_pilot.sh beside this script')
    for query in queries:
        source = manifest.parent / query['path']
        if hashlib.sha256(source.read_bytes()).hexdigest() != query['sha256']:
            parser.error('Query checksum mismatch: ' + str(source))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'inputs.json').write_bytes(manifest.read_bytes())
    shutil.copy2(script, output / 'run_pi_pilot.sh')
    shutil.copy2(Path(__file__), output / 'run_pi_batch.py')
    schedule = []
    rng = random.Random(args.seed)
    for repeat in range(1, args.repeats + 1):
        indices = list(range(len(queries)))
        rng.shuffle(indices)
        schedule.extend({'repeat':repeat, 'query_index':i, 'query_id':queries[i]['query_id']} for i in indices)
    (output / 'schedule.json').write_text(json.dumps(schedule, indent=2) + '\n')
    rows = []
    fields = ['repeat','query_id','exit_code','status','snp_status','backend','total_seconds',
              'mashpit_seconds','reference_retrieval_seconds','snp_comparison_seconds',
              'peak_rss_kib','warnings','run_directory']

    def save(state):
        with (output / 'runs.tsv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        (output / 'summary.json').write_text(json.dumps({
            'state':state, 'organism':args.organism, 'seed':args.seed,
            'distinct_queries':len(queries), 'planned_runs':len(schedule),
            'recorded_runs':len(rows), 'successful_runs':sum(r['exit_code']==0 and r['status'].startswith('COMPLETED') and r['backend']=='minimap2' for r in rows),
            'protocol':'Serial full assembly analyses; fresh output/cache per invocation; no dedicated warm-up; network retrieval included; filesystem caches not cleared.',
        }, indent=2) + '\n')

    save('RUNNING')
    for number, item in enumerate(schedule, 1):
        root = output / f"run_{number:03d}"
        root.mkdir()
        print(f"[{number}/{len(schedule)}] {item['query_id']} repeat {item['repeat']}", flush=True)
        try:
            with (root / 'console.log').open('w') as log:
                code = subprocess.call(['bash', str(script), args.organism, str(item['query_index']), str(root)], stdout=log, stderr=subprocess.STDOUT)
        except KeyboardInterrupt:
            save('INTERRUPTED')
            print('Interrupted; partial artifacts retained in ' + str(root), flush=True)
            return 130
        runs = list(root.glob(args.organism + '-isoscout-full.*'))
        row = {'repeat':item['repeat'], 'query_id':item['query_id'], 'exit_code':code,
               'status':'NO_RESULT', 'backend':'', 'run_directory':str(root.relative_to(output))}
        if len(runs) == 1:
            run = runs[0]
            if (run / 'time.txt').exists():
                row.update(read_metrics(run / 'time.txt'))
            if (run / 'analysis/result.json').exists():
                result = json.loads((run / 'analysis/result.json').read_text())
                snp = result.get('snp_resolution', {})
                row.update(status=result['status'], snp_status=snp.get('status'), backend=snp.get('backend',''),
                           warnings=json.dumps(result.get('warnings',[])))
                row.update({name+'_seconds':value['elapsed_seconds'] for name,value in result.get('stage_timings',{}).items()})
        rows.append(row)
        success = code == 0 and row['status'].startswith('COMPLETED') and row['backend']=='minimap2'
        save('RUNNING' if success else 'STOPPED_ON_FAILURE')
        print(f"  {row['status']}; total={row.get('total_seconds')} s; SNP={row.get('snp_comparison_seconds')} s", flush=True)
        if not success:
            print('Stopped; inspect ' + str(root / 'console.log'), flush=True)
            return 2
    save('COMPLETED')
    print('Saved ' + str(output / 'runs.tsv'), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
