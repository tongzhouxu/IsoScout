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

## Common finalist regions and bounded verification

Policy 4.1.0 compares saved **target-relative** evidence. A deterministic discovery
scan supplies a seed. Potentially closer references join a growing finalist panel;
all finalists are counted on one identical intersection of their aligned target
positions. Every outside reference is rechecked after the panel changes. At most
three verification passes are allowed; insufficient shared sequence or exhaustion
of this budget leaves genome resolution unresolved.

A panel requires at least **100,000 shared aligned target bases** and **50% of the
largest aligned target span among its members**. These are provisional evidence
guards, not calibrated strain cutoffs. An outside reference contributes an
observed SNP lower bound on the part of the panel region it covers. With at least
100,000 shared bases, a lower bound exceeding the panel minimum rules it out on
that mask, even below 50% overlap: unknown positions cannot reduce the observed
count. Otherwise it joins the panel when both overlap guards pass, or remains an
unresolved alternative. Missing positions are never treated as matches. Failed
or unaligned comparisons also remain unresolved.

Once the panel stops changing, one additional linear pass checks the preferred
reference against each alternative on their full pairwise shared target regions.
An adequately overlapping alternative that ties or beats the preferred reference,
but was not a panel minimum, is retained as unresolved. This detects contradictory
evidence hidden by the smaller common region. It does not select whichever region
produces a desired cluster label, and it is not a guarantee against other regional
or alignment-method effects. Exact panel ties remain ties.

Discovery, up to three verification passes and the sensitivity check require at
most `5 * (N - 1)` comparisons of cached evidence. No candidate-to-candidate
alignment, all-pairs matrix or all-candidate intersection is constructed. Panel
membership can depend on deterministic exploration within the recorded budget;
only verified evidence supports a conclusion. Comparison regions exclude reported target-coordinate deletions and ambiguous
difference spans. Raw alignment extent remains a separate coverage diagnostic.
These regions are not a fully validated callable mask, recombination-masked core
or NCBI SNP distance; unreported ambiguity and method-specific alignments can
still affect comparisons.

Cluster and genome conclusions are separate. A cluster witness requires every
reference that could meet a recorded mask's minimum, including unresolved and
region-sensitive alternatives, to carry the same stored label. A witness can
support that label while individual genomes remain unresolved. Unknown labels,
external-cluster alternatives and conflicting witness labels prevent that claim.

## Result fields and interpretation

- `candidate_challenges` records discovery, verification passes, outside-reference
  lower bounds, mask checksums, evidence hashes and unresolved alternatives.
  Its `status` describes common-panel verification, not the final resolution state.
- `finalist_panel` records the common region, per-finalist counts and tied minima.
  `finalist_snp_distance` is populated only for panel members in the display rows.
- `candidate_challenges.region_sensitivity` records the final linear check and
  contradictory alternatives. Discovery `comparisons` must not be shown as the
  final ranking.
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

Position evidence uses `target-sites-v2`; earlier extent-only evidence must be
regenerated from cached alignments. `target_comparable_bases` records the usable
comparison span; `target_aligned_bases` keeps the original extent diagnostic.
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
policy, tool hashes, jobs, cache hits, failures and timings. Reports show common-panel counts, separately labeled outside-reference lower
bounds, regional contradictions, raw diagnostic counts and coverage. A target-only comparison cannot define a phylogenetic tree.

Upstream methods: [minimap2](https://github.com/lh3/minimap2),
[paftools](https://github.com/lh3/minimap2/blob/v2.31/misc/paftools.js), and
[MUMmer4](https://github.com/mummer4/mummer).
