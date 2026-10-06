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

The versioned defaults in `config/assembly-comparison-policy.json` use a minimum
500 bp alignment and mapping quality 5. The mapping-quality default follows
upstream paftools. The 85% target/candidate coverage values are now **reporting
references only**: falling below them never removes an aligned candidate.
Neither density filtering nor contig-edge SNP masking is applied. Assembly QC
remains a separate earlier stage. `--assembly-comparison-policy PATH` records an
explicit alternative policy.

The MUMmer alternative uses `nucmer --maxmatch`, one-to-one `delta-filter -1`,
`show-coords -rclTH` and `show-snps -rlTHC`. Its policy is
`config/mummer-comparison-policy.json`. It counts ACGT substitutions in unique
mappings and reports alignment-interval coverage and indel bases. Its distances
need not equal minimap2, SKA or another filtered SNP pipeline.

## Bounded candidate challenges

Policy 3.0.0 compares two candidate references at a time using their saved
**target-relative** evidence. Intersect their aligned target intervals, then count
each reference's SNPs on those same positions. Fewer SNPs wins that challenge;
exact ties remain ties. No candidate-to-candidate alignment is performed.
Different challenges may use different regions, so there is no universal SNP
count, total distance ordering or all-candidate intersection.

A deterministic discovery scan chooses an anchor. Each verification pass compares
that anchor against every usable candidate. A losing anchor can be replaced, with
at most three verification passes and at most `4 * (N - 1)` distinct challenges.
A repeated anchor or exhausted verification budget leaves genome resolution
ambiguous. Accession order schedules exploration; only verified evidence can
support a conclusion. It can affect which unresolved comparisons are explored
within the budget.

Each challenge currently requires at least **100,000 shared aligned target bases**
and **50% of the larger of the two aligned target spans**. These are provisional,
versioned evidence guards, not calibrated strain cutoffs. A failed guard retains
that reference as an unresolved alternative; it never discards the reference to
promote a competing answer. Failed or unaligned candidates also remain unresolved
alternatives. The shared regions describe alignment extent, not an exact callable
mask, recombination-masked core genome or NCBI SNP distance.

Cluster and genome conclusions are separate. An anchor with no observed defeater
supports itself and any exact ties; unresolved alternatives still limit the final
claim. A cluster can also be supported when a checked anchor strictly beats every
examined reference outside its stored cluster, even if comparisons inside that
cluster cycle, tie or have insufficient overlap. Such a cluster witness does not
identify a unique closest genome. Unknown labels and unresolved external-cluster
references prevent that witness. Multiple conflicting witness labels fail closed.

## Result fields and interpretation

- `candidate_challenges` records anchors, passes, pair counts, shared lengths,
  support fractions, mask checksums, evidence hashes, ties and unresolved comparisons.
- `cluster_certificate` records any cluster witnesses. `cluster_status` and
  `genome_status` must be reported separately. `nearest_cluster` and `nearest_sample`
  are populated only when their respective conclusions are resolved.
- `decision_reference_samples` includes the references needed to interpret the
  cluster conclusion. `cluster_candidates` gives their known stored labels.
- `ranked` is retained as a compatibility field containing **display groups**, not
  a total distance order. `challenge_role` labels those groups. Its `snp_distance`
  is a raw pair-specific diagnostic; never rank these counts against each other.
- `nearest_snp_distance` is null and `nearest_by_aligned_snp_rate` is empty: different
  masks cannot support one comparable scalar or rate across all references.
- `coverage_flags` records low coverage without exclusion. The compatibility field
  `coverage_blockers` contains unresolved candidate evidence, including failures.

Missing, corrupt or inconsistent position evidence fails rather than assigning
zero distance. Evidence files and input identities are checksum-verified.

Every conclusion is conditional on the examined reference pool. `search_scope`
records retrieval truncation, unexamined returned references and missing downloads.
Expansion searches selected returned clusters; it cannot recover an unseen cluster.
Unresolved cluster evidence prevents sampled stability; changes in supported
clusters, genomes or resolution states reset it. If only cluster support is stable,
`stability_scope` says so. Stability does not establish a globally closest genome,
strain identity or outbreak membership. Stored labels belong to the recorded
database release and are not automatically translated to newer releases.

Per-base calls remain in the pair cache. Every round records backend, effective
policy, tool hashes, jobs, cache hits, failures and timings. Reports show paired
challenge counts and shared sequence lengths, plus raw diagnostic counts and
coverage. A target-only comparison cannot define a phylogenetic tree.

Upstream methods: [minimap2](https://github.com/lh3/minimap2),
[paftools](https://github.com/lh3/minimap2/blob/v2.31/misc/paftools.js), and
[MUMmer4](https://github.com/mummer4/mummer).
