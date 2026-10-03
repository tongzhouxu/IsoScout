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

## Bounded rounds

The initial representative-selection policy remains in
`config/snp-resolution-policy.json`. Additional rounds use
`config/refinement-policy.json` (version 1.1.0):

- Up to 25 new references per round, with at most 100 additional reference
  attempts after initial selection and a total ceiling of 200 attempted
  reference downloads, including the initial set and failed downloads.
- Up to eight rounds, including the initial comparison.
- After each SNP comparison, qualifying SNP-nearest clusters, including all ties,
  come first. Other selected returned clusters follow in best representative
  rank order. Each batch alternates across this order. Within a cluster,
  accessions are ordered deterministically. The configured `cluster_priority`
  is `nearest_then_selected_rank`; `selected_rank_only` is a reproducible
  alternative for a budget comparison. Neither ordering uses unmeasured member
  similarity as evidence.
  This ordering is a reproducible sampling choice, not a similarity ranking of
  unsketched members or a tree-neighborhood search.
- Successfully downloaded references are retained across rounds. Only new
  accessions are downloaded; each SKA2 run compares the cumulative available set.
- Continue when a closer neighbor appears, the nearest set changes, clusters
  remain tied, or some query comparisons fail comparability checks.
- Stop after two expansion rounds with an unchanged nearest set and distance,
  or at pool exhaustion, the separate expansion budget, total resource ceiling,
  the round limit, insufficient
  comparability, missing membership, or a tool failure.
  Rounds with no newly compared references do not count toward stability.
  Stability is a sampled-set stopping rule even when whole returned clusters
  or members remain unexamined; it is never evidence that they lack closer
  SNP neighbors.

Every round retains its requested accessions, downloads, SKA2 output, interpreted
results, and decision. A later failure preserves earlier successful comparisons.
The final `expansion.json` records the stop reason and the number of available
members still unexamined (unknown when discovery was unavailable). The final
`interpretation.json` lists unexamined representative and member accessions,
download failures, and qualifying versus excluded comparisons. Reaching a
stable sampled neighborhood is never reported as exhaustive search or proof
that a closer genome does not exist. Initial representative selection can already
consume the total budget; that condition is reported before membership retrieval.

## Comparability and ties

The default experimental comparability checks require at least 100,000 shared
split-kmer contexts and SKA2's mismatch proportion no greater than 0.10. These
settings are recorded and versioned; they are not organism-validated biological
or strain-assignment cutoffs. SKA2's mismatch proportion measures missing
split-kmer contexts, not nucleotide mismatch rate or ANI. See the pinned
[SKA2 implementation](https://github.com/bacpop/ska.rust/blob/v0.5.1/src/merge_ska_array.rs).

Failed comparisons remain in `excluded_comparisons` but cannot win nearest-neighbor
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
