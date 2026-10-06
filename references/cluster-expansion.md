# Cluster-member expansion and refinement feedback

`--snp-expand` implies `--snp-resolve`. It first compares the selected global
representative set, then explores additional members of every cluster represented
in that selected set, including alternative clusters beyond Mashpit's near-top
score band. This is opt-in because it retrieves public release
metadata and additional assemblies. Query sequences remain local.

## Pinned membership source

The published `databases-v1` Mashpit packages contain representatives but not
complete cluster membership. A new package can include compressed copies of
`Metadata/<PDG>.metadata.tsv` and
`Clusters/<PDG>.reference_target.all_isolates.tsv` from its exact NCBI release.
`database.json` records their relative paths, compressed and original SHA-256
checksums, sizes, release and species. Before a screen, IsoScout verifies both
compressed files. Expansion reads them directly and joins on `target_acc`;
no historical metadata download is needed for a bundled package.

For an older package without bundled tables, `expand_cluster_members.py` still
retrieves the exact release. The database's species and version are required;
the workflow never substitutes `latest` or a different release. Redirects are
rejected. Each download is limited to 1 GiB by the versioned policy. A missing
exact release produces an incomplete-expansion warning while preserving the
initial comparison. If a bundled file fails verification, the database is
rejected rather than silently switching to live NCBI data.

The limit was raised from 512 MiB in policy 1.0.1 after checking the published
database releases on 2026-10-02: Salmonella release `PDG000000002.3864` has a
786,821,557-byte metadata TSV. All five published releases' metadata and
membership TSVs responded with HTTP 200; this checks current availability,
not permanent retention. Historical snapshots still need local preservation.

Versioned assembly accessions are deduplicated. Conflicting cluster membership
is an error. Isolates missing metadata or an assembly, and requested clusters
without members, are recorded. Membership provenance includes URLs, the release,
retrieval timestamp when applicable, sizes, and SHA-256 checksums. Retrieved query/reference files
also have checksums in the final interpretation.

## Focused, bounded rounds (policy 2.0.0)

The initial representative selector remains policy 3.0.1: it includes the
global top 50 unique nonself returned references, then gives other returned
clusters coverage up to its separate 100-reference adaptive cap. The sketch
score tolerance never excludes a cluster from SNP refinement. Mashpit's default
return limit remains 200. Additional rounds use `config/refinement-policy.json`.

The default `focused_balanced_v2` allocation retains the existing limits:
25 new references per round, at most 100 additional attempts, 200 total
reference attempts including the initial set and failed downloads, and at
most eight rounds including the initial comparison. A local accession's sort
order is only a deterministic sampling order, **not** a similarity ranking.

After every comparison, all clusters tied at the minimum qualifying SNP
distance become the leading set. The scheduler normally aims for 75% of a
batch in leading clusters and reserves at least two slots for alternative
selected clusters; one of those slots can favor a cluster with no qualifying
comparison. The remaining alternative slots rotate through selected clusters
using a cursor retained between rounds, so the first clusters cannot repeatedly
consume all alternative capacity. Tied leaders also rotate within their lane.
If the finite leading pool can fit in the remaining attempts and rounds while
preserving the alternative reserve, the scheduler increases focus enough to
complete it. Unused focus or alternative capacity transfers to the other lane.
Priorities are recalculated from qualifying SNP evidence after each round.
No ground-truth cluster or known nearest accession enters scheduling.

The 75% share, two alternative slots, one unresolved-priority slot, and 20
new qualifying comparisons per leading cluster are **experimental computational
settings**. They are not validated biological thresholds or evidence that
unexamined clusters are distant. The versioned policy records them and every
round's `allocation` lists chosen accessions, lane, leading and unresolved
clusters, available capacity, cursors, and rationale. Downloads, comparisons,
exclusions, and newly qualifying references are recorded separately.

A sampled-stability stop requires two unchanged rounds **after** every leading
cluster has either at least 20 qualifying additional comparisons or no
unexamined available member, and after any feasible finite leading pool has
been completed. A round with no newly qualifying comparisons does not count
as stable. A single comparability exclusion in an otherwise comparable cluster
is retained but does not itself force more expansion. A cluster with no
qualifying comparison remains unresolved. If it has untried available members,
stability is blocked and alternative exploration continues. If its finite pool
is exhausted without a qualifying comparison, the unresolved status remains
in the report; it is never treated as SNP-distant.

Other stop reasons distinguish finite-pool exhaustion, insufficient evidence
at pool exhaustion, sampled stability, missing or conflicting membership,
round/tool failure, the additional-attempt budget, and the total ceiling.
`expansion.json` records remaining per-cluster and overall coverage, including
where budget or round limits prevented the focused minimum or exploration of
unresolved alternatives. Stability is only within the sampled set; neither
cluster concordance nor a stable nearest reference proves the globally nearest
genome was found.

## Explicit policy files

`--refinement-policy PATH` can select another versioned policy file for a
screen. Its expansion and comparability settings are applied together. The
selected file's path, version, content, and SHA-256 are recorded in the
expansion audit and provenance. The checked-in default keeps the limits above;
any different resource ceiling must be stated in the supplied policy file.

## Comparability and ties

The default experimental comparability checks require at least 100,000 shared
split-kmer contexts and SKA2's mismatch proportion no greater than 0.10. These
settings are recorded and versioned; they are not organism-validated biological
or strain-assignment cutoffs. SKA2's mismatch proportion measures missing
split-kmer contexts, not nucleotide mismatch rate or ANI. See the pinned
[SKA2 implementation](https://github.com/bacpop/ska.rust/blob/v0.5.1/src/merge_ska_array.rs).

Failed comparisons remain in `excluded_comparisons` with their reasons but cannot win nearest-neighbor
ranking. No qualifying comparison yields `INSUFFICIENT_DATA`. A unique nearest
reference yields `COMPARED`; equally nearest references yield `AMBIGUOUS`, with
all nearest accessions and clusters retained. A cross-cluster tie has no singular
nearest cluster. Zero observed SNPs does not establish whole-genome identity.
Trees require a complete matrix of qualifying comparisons among included tips.

## Read input

For paired-end input, fastp's cleaned R1/R2 files are now used directly for SKA2
refinement, alongside downloaded reference assemblies. The fixed profile uses
`--min-count 3 --min-qual 20 --qual-filter strict`; these are starting error-filter
settings requiring real-data validation. The pinned [SKA2 read-input contract](https://docs.rs/ska/0.5.1/ska/#read-files)
requires two FASTQ paths on the same sample's file-list row. Single FASTQ paths
are rejected to avoid silently treating reads as assemblies.

Mashpit and MLST routing still use the SKESA assembly. Thus users can start from
reads, and the fine-resolution comparison uses reads directly, but the whole
screen is not assembly-free. Assembly failure still stops the initial screen.

## Assembly backend

Assembly inputs default to cached minimap2/paftools target-to-candidate
comparisons. The split-kmer criteria above apply only to SKA. Assembly coverage,
variant rules, and alternative MUMmer settings are documented in
[assembly-comparison.md](assembly-comparison.md). New batches reuse previous
alignments; expansion does not require a pairwise matrix or tree.
