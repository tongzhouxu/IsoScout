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

## Shared-region ranking and unresolved alternatives

Policy 2.0.0 ranks all qualifying references using the **same target positions**:
the intersection of their recorded target alignment intervals. SNPs are recounted
from checksum-verified per-position evidence; no additional alignment or
candidate-to-candidate matrix is constructed. Every reference uses the same
denominator. Raw pair-specific SNP counts and rates remain diagnostic records
because similar coverage percentages do not imply that the same positions were
compared. No reference is removed to improve the intersection or its ranking.

The shared mask must satisfy the policy's existing minimum target-coverage
fraction (85% by default). This conservative default is an explicit software
criterion, not an empirically calibrated strain threshold. If it is not met,
return `INSUFFICIENT_DATA`, with no nearest sample or cluster. Expansion stops at
`insufficient_shared_target_regions`: adding candidates cannot enlarge the
intersection. This can increase unresolved outcomes; there is no silent fallback
to incomparable pair counts. The shared mask describes aligned sequence, not an
exact callable-site mask, recombination-masked core genome or NCBI SNP distance.

A coverage-excluded reference with no larger pair-specific SNP count **or** rate
than a shared-region winner is a potentially competitive alternative. It is not
promoted into the ranking, but prevents a resolved nearest-genome claim. If its
stored label differs from the winning label, or is unknown, the cluster conclusion
also remains ambiguous. The guard identifies uncertainty; it does not establish
that the excluded reference is closer. A same-cluster alternative can leave the
cluster resolved while the genome remains unresolved.

`ranked[].ranking_snp_distance` and `nearest_snp_distance` contain shared-region
counts. `ranked[].snp_distance` retains the original pair-specific count.
`shared_regions` records the exact mask, its checksum, target identity, evidence
checksums, participating references, aligned bases and coverage. Missing, corrupt
or inconsistent position evidence fails rather than assigning zero distance.
`pair_specific_diagnostics` retains the former count/rate minima for audit only.
The compatibility field `nearest_by_aligned_snp_rate` now names the shared-region
rate minima, which necessarily equal the count minima on the same denominator.

`cluster_status` and `genome_status` describe separate conclusions. Preserve exact
genome ties and all cross-cluster ambiguity. `coverage_blockers` records competitive
excluded alternatives, and `cluster_candidates` includes their known labels.
`nearest_cluster` is populated only when the cluster conclusion is resolved;
`nearest_sample` only when the genome conclusion is resolved. `nearest_samples`
and `nearest_clusters` retain the qualifying minimum sets even when an excluded
alternative prevents resolution; never use those fields alone as a resolved call.

Every result is conditional on the examined reference pool. `search_scope`
records retrieval truncation, unexamined returned references and missing downloads.
Expansion searches selected returned clusters; it does not prove that an unseen
cluster is absent or more distant. Sampled stability is reset when the common
mask changes and is blocked by competitive exclusions. Neither a stable result
nor a shared-region minimum establishes global nearest-genome recovery, strain
identity or outbreak membership. Stored labels belong to the recorded database
release and are not automatically translated to newer releases.

Detailed per-base calls remain in the pair cache alongside raw alignment files;
round and user-facing summaries contain compact per-candidate measurements.
Every round records backend, effective policy, tool hashes, successful jobs,
cache hits, failures and timings. Reports separate shared-region and pair-specific SNP counts, with target and
candidate coverage; the full interpretation retains the pair-level records and
all exclusions. A target-only distance list cannot be used to invent a complete
matrix or phylogenetic tree.

Upstream methods: [minimap2](https://github.com/lh3/minimap2),
[paftools](https://github.com/lh3/minimap2/blob/v2.31/misc/paftools.js), and
[MUMmer4](https://github.com/mummer4/mummer).
