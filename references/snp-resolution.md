# SNP resolution: selection policy 3.0.1

Optional, opt-in refinement of a Mashpit candidate using ska2 pairwise SNP distances. Governed by `config/snp-resolution-policy.json` (versioned resource policy) and the `snp_resolution` block of `config/workflow.json` (fixed ska2 command profile). Enable with `--snp-resolve`.

## Why this exists

Mashpit's MinHash similarity is a coarse, sketch-resolution screen. A single unambiguous top cluster at Mash resolution is not proof of SNP-level closeness, and two near-tied clusters cannot be told apart by Mash alone. SNP resolution runs ska2 (split k-mer analysis) between the query assembly (or cleaned paired reads) and returned reference genomes across clusters to get an actual pairwise SNP count.

## When it runs

Whenever `--snp-resolve` is set and Mashpit returned a candidate (`mashpit_result.best_candidate` is present) — regardless of whether that candidate was unambiguous. This is deliberate: an unambiguous Mash hit can still be many SNPs away from its representatives, so the refinement step is not limited to the ambiguous case.

## Network dependency

Unlike the rest of the screen, this step is **not** fully local. A Mashpit database only retains sourmash signatures, not the representative assemblies (they are sketched and discarded during `mashpit build`), so resolving SNPs requires re-downloading selected representative genomes from NCBI via the pinned `datasets` CLI. This is why the step is opt-in rather than automatic: default screens keep sensitive query data fully local, and `--snp-resolve` is an explicit choice to reach out to NCBI for public reference genomes (the query sequence itself is never uploaded).

## Target selection and budgets

The default `global-ranked-cluster-coverage-v1` policy (3.0.1) uses an initial
**global** target of 50, a maximum of 100 returned representatives in adaptive
mode, at most 100 additional expansion attempts under the default focused
refinement profile, and a total ceiling of 200
attempted reference downloads including later member expansion. Mashpit still
returns at most 200 representatives by default.
These are separate recorded computational budgets, not biological cutoffs.

1. Validate scores and collapse identical accession rows. Conflicting rows fail.
   Sort all unique returned representatives by descending Mashpit score, then
   accession. Supply `--query-accession` when the query's assembly accession is
   known; an input filename that is itself a versioned `GCA_`/`GCF_` accession
   is also recognized. It is excluded before counting the first 50 and from
   later member expansion. Sequence-identical
   records under another accession cannot be detected from Mashpit scores alone.
2. If Mashpit's top candidate passes its configured screening gate, include the
   first 50 valid, nonself representatives **regardless of cluster or score
   gap**. With fewer than 50, include all available. The sketch tolerance
   (`tie_tolerance_hashes / hash_number`) and upstream `near_top` labels remain
   annotations only; neither excludes a SNP reference.
3. Up to the adaptive returned-reference cap of 100, first include the highest
   ranked representative of each cluster absent from the initial set. Then fill
   any remaining slots in global rank order. Cluster leaders are visited in
   global rank order. This is a deterministic coverage strategy, not a claim
   that a lower scored cluster is biologically implausible.
   The versioned `coverage_strategy` can instead be set to `global_rank_only`.
   Both strategies preserve the first 50; neither uses a score-gap exclusion.
4. `--snp-selection-mode all_returned` requests every valid unique returned
   representative. It requires the total reference
   budget to fit the entire returned set; otherwise selection is `SKIPPED` with
   `all_returned_budget_insufficient` decisions. It never silently turns into a
   partial set. With 200 returned references and the default 200-attempt
   ceiling, it attempts all 200, leaving no member-expansion capacity.

`--snp-resolve` uses adaptive mode by default. `--snp-expand` additionally
explores exact-release members of clusters represented in the selected returned
set. The expansion profile is versioned separately from representative selection;
`--refinement-policy` can select another explicit, versioned policy file without
changing Mashpit retrieval or the initial candidate set. Its expansion and
comparability settings are applied together and recorded. No score-distribution
heuristic changes eligibility or budgets; exploratory plots and near-tie flags
must not be interpreted as evidence that omitted clusters lack closer SNP
neighbors.

`snp_resolution/targets.json` includes the policy and recorded Mashpit query
settings, source checksums, mode, all selected accessions, each unique returned
reference's global rank and inclusion/omission reason, query-self exclusion,
duplicate-row decisions, selected and unexamined clusters, and explicit budgets.
`returned_reference_budget` means a resource omission, never SNP dissimilarity.
The return-limit flag means additional matches *may* be unobserved. The preview
in `similarity_distribution/summary.json` calls the same selector with the same
mode and query accession; it predicts requests, not downloads or comparisons.

For standalone selection, use `select_snp_targets.py --mashpit-output-dir ...
--policy config/snp-resolution-policy.json --mode adaptive|all_returned
--output ...` and optionally `--query-accession`. The saved Mashpit directory
must contain `mashpit_run.json`; no current settings are guessed.

## Genome retrieval and SNP distance

Representative genomes are downloaded with `datasets download genome accession --include genome --dehydrated` followed by `datasets rehydrate` — one batched call for every selected accession, not one request per genome — retried up to `download_attempts` times for genomes that fail. The query assembly or cleaned read pair plus every successfully downloaded reference are built into one ska2 split-kmer file (`ska build -f <name-path list> -k <kmer_size>`, pinned k-mer size in `workflow.json`), then compared with `ska distance`, which reports the number of SNPs differing between every pair — the *entire* pairwise matrix (every representative against every other, not just against the query).

## Interpretation

`interpret_snp_resolution.py` reports qualifying query-to-reference comparisons,
excluded low-overlap comparisons, all tied nearest references, per-cluster
summaries, and an exploratory Neighbor-Joining tree when a complete comparable
matrix is available. Status is `COMPARED`, `AMBIGUOUS`, or `INSUFFICIENT_DATA`;
none establishes strain identity. The structural `confidence` block contains
raw distances and ratios, not statistical confidence. Zero observed differences
are not described as an exact whole-genome match. See
[cluster-expansion.md](cluster-expansion.md) for comparability settings and ties.

The final `coverage` block lists every returned representative left unexamined,
every attempted, verified, and unavailable download, every qualifying and
comparability-excluded query comparison, and discovered but unexamined members.
It also flags when member coverage is unknown because membership was not read.
Per-round allocation and newly qualifying comparisons are in `expansion.json`;
the final report separates cluster-label concordance from recovery of an
individual nearest genome among those examined.
The report separates Mashpit-displayed alternatives from alternatives actually
tested by SKA2. No stable within-cluster result rules out an untested cluster.

`render_snp_tree.py` renders `newick_tree` to `snp_resolution/tree.png` with `QUERY` highlighted in red and bold (via `phytreeviz`, already pulled in transitively by the pinned mashpit commit's own dependencies — no new container pin needed). Rendering is best-effort: a failure only sets `tree_image.status` to `FAIL` with the error message, it never fails SNP resolution or the underlying Mash screen. `generate_report.py` embeds this image in `report.md` when available, falling back to the raw Newick text otherwise.

`generate_report.py` renders all of this into `report.md`, written alongside `result.json` on every run (not just `--snp-resolve` ones): organism determination, the Mashpit candidate table, the ska2 cluster/genome tables and confidence statement, and the Newick tree, in plain language for a reader who doesn't want to open JSON.

## Limitations

- ska2 pairwise SNP distance is itself a screening refinement, not a validated outbreak-confirmation pipeline (unlike, e.g., CFSAN SNP Pipeline or an accredited cgMLST scheme). Always recommend a validated, organism-appropriate high-resolution comparison for actual outbreak confirmation, exactly as for the Mash result.
- `--snp-resolve` compares selected representatives. `--snp-expand` can add available members of selected clusters from the exact release, subject to explicit budgets and stopping rules. Neither mode searches all public genomes or guarantees the nearest isolate globally.
- A failed or partial download (fewer than two usable genomes, including the query) skips SNP resolution with a warning; it does not fail the underlying Mash screen, which remains authoritative on its own.
