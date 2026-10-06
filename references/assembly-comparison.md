# Target-to-candidate assembly comparison

Assembly refinement defaults to minimap2 2.31 (`2.31-r1302`) with its distributed
`paftools.js` and k8 1.2. MUMmer4 4.0.0 is an explicit alternative. No assembly
backend calculates candidate-to-candidate distances or constructs a tree.

## Local assemblies and reuse

`--assembly-manifest PATH` accepts a tab-separated file with `accession`, `path`
and `sha256` columns. Paths may be absolute or relative to the manifest. The
mapping is supplied by the user; checksums verify the supplied sequence content,
not the biological correctness of the accession annotation. Duplicate IDs,
missing files and checksum mismatches fail rather than downloading substitutes.
Only accessions absent from this manifest are downloaded.

`--comparison-cache DIR` retains raw alignments under keys containing ordered
input checksums, alignment settings and executable identities. Expansion processes
only new pairs. Cache files are verified before reuse. The minimap2 backend caches
variant interpretation separately, so changing coverage criteria or variant
filters does not require realignment. Input or alignment changes invalidate the
appropriate alignment entry. Interrupted commands never produce a successful
comparison record. Re-run a screen into a new output directory with the same
cache to reuse completed comparisons without overwriting the interrupted run.

To compare every entry in an explicit local candidate set directly, bypassing
Mashpit retrieval and expansion, use:

```bash
python3 scripts/run_minimap.py --target target.fna \
  --candidate-manifest candidates.tsv --output comparisons --cache-dir pair-cache
```

The same arguments are available in `scripts/run_mummer.py`. These standalone
commands return pair measurements; cluster interpretation requires the screening
workflow's accession-to-cluster mapping.

`--comparison-workers` sets concurrent single-thread pair alignments (default 4).
When launching several screens, their combined worker count must fit the available
CPU and memory budget. Each command has a recorded timeout; failed candidates
remain explicit and are not interpreted as zero differences.

## Distance and coverage

The target assembly is always the alignment reference. Minimap2 uses `asm5`,
base-level alignment, a short `cs` tag and no secondary alignments. Eligible PAF
records are sorted in target coordinates and passed to upstream `paftools.js call`.
Only ACGT substitutions with single target coverage contribute to SNP counts.
Indel bases are recorded separately. The target coverage denominator comes from
paftools' uniquely covered regions; candidate coverage is the union of qualifying
PAF query intervals. These are alignment-extent measurements, not exact callable
site counts. Ambiguous bases may be omitted upstream and zero reported ambiguous
difference rows does not establish absence of ambiguous input sequence.

The versioned defaults in `config/assembly-comparison-policy.json` require 85%
target and candidate alignment coverage, a minimum 500 bp alignment and mapping
quality 5. The mapping-quality default follows upstream paftools; the coverage
and alignment-length criteria are explicit IsoScout defaults for assembly
comparison. They are not user-specified criteria, calibrated strain thresholds,
or inherited SKA shared-kmer tests. `--assembly-comparison-policy PATH` records an
explicit alternative policy. Neither density filtering nor contig-edge SNP
masking is applied. Assembly QC remains a separate, unchanged earlier stage.

The MUMmer alternative uses `nucmer --maxmatch`, one-to-one `delta-filter -1`,
`show-coords -rclTH` and `show-snps -rlTHC`. Its policy is
`config/mummer-comparison-policy.json`. It counts ACGT substitutions in unique
mappings and reports alignment-interval coverage and indel bases. Its distances
need not equal minimap2, SKA or another filtered SNP pipeline.

Rank qualifying candidates by observed SNP count, preserving all exact ties.
Also report the minimum SNP rate per aligned target Mb. If the minimum sets do
not overlap, mark the closest-genome conclusion unresolved across these metrics.
Never blend them into an unvalidated composite score or use coverage to silently
break a SNP-count tie. The result is the nearest by the declared metric among
examined qualifying candidates; it is not an exhaustive global nearest-genome
answer, whole-genome identity, or outbreak confirmation.

Detailed per-base calls remain in the pair cache alongside raw alignment files;
round and user-facing summaries contain compact per-candidate measurements.
Every round records backend, effective policy, tool hashes, successful jobs,
cache hits, failures and timings. Reports show candidate SNPs with target and
candidate coverage; the full interpretation retains the pair-level records and
all exclusions. A target-only distance list cannot be used to invent a complete
matrix or phylogenetic tree.

Upstream methods: [minimap2](https://github.com/lh3/minimap2),
[paftools](https://github.com/lh3/minimap2/blob/v2.31/misc/paftools.js), and
[MUMmer4](https://github.com/mummer4/mummer).
