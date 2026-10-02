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
The published container has not been updated with these changes.

The test audit consolidated four redundant tests from the previous 90-test
suite. Genome-cap checks now live with incomplete cluster coverage; full distance
parsing assertions live with the mocked SKA2 runner; cluster preview checks live
with the preview-to-download workflow; saved diagnostic errors are checked by
the workflow failure test. Distinct score/cluster conflicts and parser fields
were retained. Before/after tracing of the default suite covered the same 1,319
production lines and 1,871 execution transitions in process. Three subsequent
database-package tests cover offline membership, tampering, and release-version
mismatch. The optional real SKA2 test remains separate from mocked tests.

## Held-out finite-pool benchmark

`scripts/benchmark_refinement.py` compares the adaptive selection with exhaustive
SKA2 comparison in a supplied, finite reference collection. It uses local files
and does not download genomes. Its JSON output includes nearest-neighbor recall,
retention of all tied nearest neighbors, selected/pool counts, comparison
exclusions, checksums, tool commands, and optional agreement with supplied labels.

Example manifest (paths are resolved relative to the manifest):

```json
{
  "held_out_query": true,
  "query": ["query_R1.fastq.gz", "query_R2.fastq.gz"],
  "mashpit_output_dir": "previous-screen/mashpit",
  "references": [
    {"accession": "GCA_000001.1", "cluster": "PDS000001.1", "path": "ref1.fna"},
    {"accession": "GCA_000002.1", "cluster": "PDS000001.1", "path": "ref2.fna"}
  ],
  "expected_nearest": ["GCA_000001.1"],
  "label_source": "Describe the independent evidence for these optional labels"
}
```

`query` may instead be a single assembly path. The historical Mashpit output
must contain `mashpit_run.json`. Include every adaptively selected accession in
the exhaustive pool. Remove the query's own record and duplicate-isolate records
before setting `held_out_query`; a byte-identical query/reference file is rejected,
but byte checks cannot detect the same isolate represented by different files.

```bash
python3 scripts/benchmark_refinement.py \
  --manifest /path/to/manifest.json --output-dir /path/to/new-benchmark
```

The adaptive subset reuses the exhaustive pairwise distances: the fixed SKA2
profile uses no cohort-frequency filter. This evaluates whether the selection
retains nearest neighbors under the same method, not whether SKA2 is biological
ground truth. Optional labels are explicitly attributed to their source.

## What still requires scientific validation

Before claiming strain-level assignment accuracy, evaluate matched real reads
and assemblies, independent organism-appropriate reference-mapping or cgMLST
results, dense neighborhoods, incomplete genomes, mixed/contaminated inputs,
novel isolates, and queries without a close database neighbor. Compare fixed
candidate counts and the adaptive policy across those cases, retaining all tied
nearest neighbors and reporting unresolved results. Tune computational budgets
and comparability settings on training data and evaluate held-out collections.
The tests supplied here do not substitute for that biological validation.
