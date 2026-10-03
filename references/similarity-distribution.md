# Similarity-distribution diagnostics

Workflow 1.3.0 introduced diagnostics after each successfully parsed Mashpit query,
including screens without `--snp-resolve`. These describe the returned
representatives and preview the configured SNP selection policy. They do not
change candidate selection, expand retrieval, or download genomes.

## Outputs

- `similarity_distribution/summary.json`: ranked scores and adjacent gaps,
  counts within tolerance of the best score, cluster composition, checks at
  ranks 50/100/200, return-limit flags, and the versioned selector's preview.
- `similarity_distribution/rank_similarity.png`: rank versus Jaccard similarity,
  with the top tolerance band and current selection preview highlighted.
- `result.json`: the same diagnostics under `similarity_distribution`.
- `report.md`: numerical summary, cutoff checks, plot, and a link to the JSON.

The JSON records source CSV checksums, sketch size, tolerance, requested return
limit, and selection policy. Full ranked values are preserved rather than rounded
to report precision. Counts describe unique assembly accessions. Identical
duplicate rows are collapsed with a warning; conflicting duplicates, malformed
scores, and missing required metadata produce a diagnostic error. Workflow 1.7.0 uses these same validated, deduplicated records in the global selector.

## Interpretation

The score tolerance is `tie_tolerance_hashes / hash_number`, matching
`generate_cluster_table` in the [pinned Mashpit source](https://github.com/tongzhouxu/mashpit/blob/538d3421302fe6dd129780605b8ff5dedbf4c046c/src/mashpit/query.py).
The actual database sketch size is required; it is never assumed to be 1000.
This is an experimental descriptive sketch-resolution heuristic, not a statistical confidence
interval, ANI estimate, or validated biological cutoff. Tiny floating-point
roundoff at the inclusive boundary is tolerated.

`within_top_tolerance_count` compares each score directly with the best score;
it does not chain together a succession of small adjacent gaps. Rank checkpoint
flags compare ranks N and N+1. If N+1 was not returned, the gap and tie flag are
unknown (`null`), not evidence for a natural cutoff. The checkpoints are descriptive
probes configured in `workflow.json`; none is selected automatically.

`retrieval.limit_reached` means the count of unique returned accessions meets or
exceeds the requested limit. It signals possible censoring, not proof that the
database contains more matching representatives. `additional_matches_exist` is
therefore unknown. If the near-top band reaches that limit, or the final adjacent
scores are within tolerance, `may_split_near_tie` is true. A small tail gap does
not imply that the tail belongs to the top band. Not reaching the limit does not
establish coverage of all isolates, because the database stores representatives.

The selection preview invokes the same mode and query-accession exclusion as the
actual selector described in [snp-resolution.md](snp-resolution.md). Adaptive
mode includes the global top 50, then uses up to 50 more returned-reference
slots for unrepresented cluster leaders and global rank. All-returned mode
requires enough budget for the entire unique returned set. The full audit is
retained under `selection_preview.audit`, including each inclusion or resource
omission. It records selected and omitted counts per cluster, unexamined whole
clusters, and omitted near-top references. `selection_splits_near_tie` is a
descriptive flag; neither it nor any score gap rules out a cluster at SNP level.
The preview never asserts that a reference was downloaded or compared.

Warnings propagate to the screen's overall warning status. Diagnostics errors
and plot failures do not prevent an otherwise valid Mashpit/SNP workflow from
continuing. Plotting uses matplotlib (already a dependency of the pinned
Mashpit toolchain); numerical diagnostics remain available if it cannot render.

## Analyze a previous run

In the configured environment, point the standalone script at an existing
screen's `mashpit/` directory and choose a new output directory:

```bash
python3 scripts/analyze_similarity_distribution.py \
  --mashpit-output-dir /path/to/screen/mashpit \
  --output-dir /path/to/new-diagnostics
```

It reads the original database settings and query flags from `mashpit_run.json`,
and previews the current repository's SNP selection policy (recorded in the new
summary). Pass `--selection-mode all_returned` and `--query-accession` to match
a run using those options. It does not rewrite the original result, report, or provenance. A new
container build or a checkout with the updated scripts is needed; previously
published images do not acquire this feature automatically.

## Validation scope

Synthetic fixtures exercise dense ties, clear gaps, chained adjacent gaps,
return-limit censoring, cluster balancing, duplicates, malformed input, and
empty/singleton results. Integration tests check screening/report wiring and
standalone use of recorded query settings. These verify diagnostic behavior;
they do not establish biological accuracy of strain assignments or the adaptive selection policy.
