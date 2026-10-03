# Refinement verification and benchmarking

`tests/test_ska_integration.py` runs the real pinned SKA2 binary on seeded random
DNA with known, separated differences. It verifies assembly and paired-read
recovery of 1-SNP and 3-SNP distances and nearest-neighbor agreement. It also
exercises the finite-pool benchmark harness. This is a software correctness test,
not evidence of strain-assignment accuracy for real organisms.

Run the standard tests without biological tools:

```bash
PYTHONPYCACHEPREFIX=/tmp/isoscout_pycache python3 -m unittest discover -s tests -v
```

In the configured container, enable the offline SKA2 integration test with:

```bash
ISOSCOUT_TOOL_TESTS=1 python3 -m unittest discover -s tests -v
```

Verification on 2026-10-02: all 89 tests passed in the existing pinned-tool
container with the updated checkout mounted at `/opt/isoscout`, networking
disabled, and `ISOSCOUT_TOOL_TESTS=1`. This includes the real SKA2 assembly/read
test above. The host-only suite passed with that one integration test skipped.
For the current policy change, the 2026-10-03 host-only suite discovered 91 tests: 90 passed and the real-SKA2 test was skipped because `ska` is unavailable on this host. The published container has not been updated with these changes.

The test audit consolidated four redundant tests from the previous 90-test
suite. Genome-cap checks now live with incomplete cluster coverage; full distance
parsing assertions live with the mocked SKA2 runner; cluster preview checks live
with the preview-to-download workflow; saved diagnostic errors are checked by
the workflow failure test. Distinct score/cluster conflicts and parser fields
were retained. Before/after tracing of the default suite covered the same 1,319
production lines and 1,871 execution transitions in process. Three subsequent
database-package tests cover offline membership, tampering, and release-version
mismatch. The optional real SKA2 test remains separate from mocked tests.

## Reproducible policy comparison

The pilot commit `6f4905e4786bfc67f499d0ad5c1acd5aa6f4a2e5` is the
frozen starting point for this change. `scripts/frozen_select_snp_targets_v2.py`
is an unmodified copy of its selector, paired with
`config/frozen-snp-resolution-policy-v2.json`. The current policy is 3.0.0.
Do not alter original pilot output directories. They are inputs to a new
comparison artifact.

For saved pilot Mashpit results, create a manifest with `data_role` set to
`development` and case entries such as:

```json
{
  "data_role": "development",
  "cases": [
    {"case_id": "07", "mashpit_output_dir": "case-07/mashpit", "expected_cluster": "PDS_example"}
  ]
}
```

Run `python3 scripts/compare_selection_policies.py --manifest
/path/to/pilot-manifest.json --output /path/to/new-comparison.json`. Each saved
Mashpit directory must contain its original `mashpit_run.json`, candidate CSV,
and representative CSV. The tool invokes frozen v2, revised adaptive, and
all-returned selection with the original recorded query context. It reports
selected accessions, clusters, resource omissions, exact policy/source checksums,
and optionally how many supplied expected-cluster references each policy
retained. The expected cluster is read **after** all selections and is never
passed into a selector. This is selection coverage, not SNP or label accuracy.
If pilot files are unavailable, do not fill this table from the five stated
summary observations or imply it was run.

`scripts/benchmark_refinement.py` adds an actual finite-pool SKA2 comparison
of the same three policies. Its manifest requires `held_out_query: true`, a
local query, the saved Mashpit output directory, and local reference genome
files. `query_accession` should be supplied when known; remove the query's own
record and duplicate isolates before declaring it held out. The exhaustive
pool must contain every accession selected by each policy. All-returned mode
must fit its recorded resource budget. The finite benchmark pool also has a
recorded `max_reference_genomes` ceiling (defaulting to the versioned total
reference budget of 200); larger pools are rejected rather than run silently. SKA2 compares the finite pool once;
policy-specific subsets reuse the same pairwise matrix. The JSON contains each
policy's selected set, qualifying and excluded comparisons, all nearest ties,
and recall of the finite pool's nearest references. The latter is a search
coverage measure under SKA2, not biological ground truth.

A minimal benchmark manifest shape is:

```json
{
  "held_out_query": true,
  "query": ["query_R1.fastq.gz", "query_R2.fastq.gz"],
  "mashpit_output_dir": "screen/mashpit",
  "references": [
    {"accession": "GCA_000001.1", "cluster": "PDS000001.1", "path": "ref1.fna"}
  ],
  "label_source": "Independent evidence for optional expected_nearest labels"
}
```

Run `python3 scripts/benchmark_refinement.py --manifest
/path/to/benchmark-manifest.json --output-dir /path/to/new-benchmark`.
The pilot cases are development data for checking exclusion and tuning budgets.
Reserve independent isolates and releases, with matched reads and assemblies,
for confirmatory evaluation. Predeclare budget, comparability settings,
reference snapshot, and outcome measures before inspecting their labels.
No improved accuracy is claimed from the software tests or selection audit.

A synthetic regression reconstructs only the five supplied pilot ranks, score
gaps, and expected-cluster reference counts. It is not the saved pilot dataset.
The frozen/revised/all-returned selectors retain these counts of synthetic
expected-cluster references, respectively:

| Case | First rank | Gap | Frozen v2 | Revised adaptive | All returned |
|---|---:|---:|---:|---:|---:|
| 07 | 3 | 0.005 | 0/14 | 14/14 | 14/14 |
| 08 | 96 | 0.016 | 0/1 | 1/1 | 1/1 |
| 11 | 12 | 0.006 | 0/5 | 5/5 | 5/5 |
| 12 | 102 | 0.026 | 0/4 | 1/4 | 4/4 |
| 15 | 3 | 0.004 | 0/74 | 74/74 | 74/74 |

This verifies the selection mechanism only. Actual pilot outputs must be run
through the comparison manifest before reporting their real retained counts,
and broader SKA2 comparison is needed before reporting any label recovery.

## What still requires scientific validation

Before claiming strain-level assignment accuracy, evaluate matched real reads
and assemblies, independent organism-appropriate reference-mapping or cgMLST
results, dense neighborhoods, incomplete genomes, mixed/contaminated inputs,
novel isolates, and queries without a close database neighbor. Compare fixed
candidate counts and the adaptive policy across those cases, retaining all tied
nearest neighbors and reporting unresolved results. Tune computational budgets
and comparability settings on training data and evaluate held-out collections.
The tests supplied here do not substitute for that biological validation.
