# SNP resolution: selection policy 2.0.0

Optional, opt-in refinement of a Mashpit candidate using ska2 pairwise SNP distances. Governed by `config/snp-resolution-policy.json` (selection thresholds) and the `snp_resolution` block of `config/workflow.json` (fixed ska2 command profile). Enable with `--snp-resolve`.

## Why this exists

Mashpit's MinHash similarity is a coarse, sketch-resolution screen. A single unambiguous top cluster at Mash resolution is not proof of SNP-level closeness, and two near-tied clusters cannot be told apart by Mash alone. SNP resolution runs ska2 (split k-mer analysis) between the query assembly (or cleaned paired reads) and reference genomes from the relevant cluster(s) to get an actual pairwise SNP count.

## When it runs

Whenever `--snp-resolve` is set and Mashpit returned a candidate (`mashpit_result.best_candidate` is present) — regardless of whether that candidate was unambiguous. This is deliberate: an unambiguous Mash hit can still be many SNPs away from its representatives, so the refinement step is not limited to the ambiguous case.

## Network dependency

Unlike the rest of the screen, this step is **not** fully local. A Mashpit database only retains sourmash signatures, not the representative assemblies (they are sketched and discarded during `mashpit build`), so resolving SNPs requires re-downloading the relevant representative genomes from NCBI via the pinned `datasets` CLI. This is why the step is opt-in rather than automatic: default screens keep sensitive query data fully local, and `--snp-resolve` is an explicit choice to reach out to NCBI for public reference genomes (the query sequence itself is never uploaded).

## Target selection

The default policy is `adaptive-boundary-v1`, version `2.0.0`. It starts with
50 references and has a hard ceiling of 200 references (excluding the query).
These are computational defaults awaiting biological benchmarking, not
validated strain-assignment cutoffs. The earlier per-cluster cap is removed.

1. Validate scores and deduplicate exact assembly accessions. Distinct accessions
   remain distinct even when scores are identical. Conflicting records for the
   same accession fail selection rather than choosing one silently.
2. Retain the top cluster, upstream `near_top` clusters, and clusters with a
   representative within the inclusive top score tolerance. The tolerance uses
   the actual database hash count and recorded query's hash tolerance. Candidate
   clusters with no available returned representative are reported explicitly.
3. Sort eligible representatives by decreasing similarity, then accession. Start
   with up to 50, extending the ranked prefix as needed to include the leading
   representative of every plausible cluster.
4. Extend the boundary while each adjacent score gap is at most the Mashpit
   tolerance. Stop at a larger gap or the end of the returned eligible pool.
   This boundary expansion intentionally follows adjacent near-ties; it can span
   more than one tolerance from the best score. The diagnostics' top-band count
   still compares every score directly with the best score.
5. If the desired set exceeds 200, reserve cluster leaders first, then prioritize
   higher score bands. Share slots round-robin across clusters within a band cut
   by the ceiling. Cluster order is deterministic (leader score, then cluster
   identifier), as is within-cluster order (score, then accession). If even one
   leader per cluster cannot fit, list the omitted clusters explicitly. Never
   present a resource-driven tie break as biological evidence.

`snp_resolution/targets.json` records policy and query settings, source checksums,
all selected accessions, the desired count before capping, the boundary gap,
missing clusters, and an inclusion/exclusion reason for every unique accession.
Duplicate input rows have separate `duplicate_accession` decisions pointing to
retained source records. Other reason codes distinguish initial ranked targets,
cluster coverage, near-tie expansion, the hard ceiling, scores beyond the boundary,
and clusters outside the plausible set.

`desired_set_complete` refers only to this policy's desired set of returned
representatives. A separate retrieval-limit flag warns when additional candidates
may be unobserved, even if every desired returned genome fits. The initial selector does not re-query Mashpit. Optional `--snp-expand` adds
cluster-member exploration after this initial stage; see [cluster-expansion.md](cluster-expansion.md).
No optimal count or nearest-neighbor guarantee is claimed.

The similarity diagnostics preview the same selector before downloading. The SNP
result retains the selection audit and lists selected references that could not
be downloaded, so selection coverage and actual comparison coverage remain visible.

For standalone selection, use the existing `select_snp_targets.py` CLI with
`--mashpit-output-dir`, `--policy`, and `--output`. The input directory must include
`mashpit_run.json` from the original screen, which supplies the recorded query
flags and database settings. Version 1 policies are rejected explicitly.

## Genome retrieval and SNP distance

Representative genomes are downloaded with `datasets download genome accession --include genome --dehydrated` followed by `datasets rehydrate` — one batched call for every selected accession, not one request per genome — retried up to `download_attempts` times for genomes that fail. The query assembly or cleaned read pair plus every successfully downloaded reference are built into one ska2 split-kmer file (`ska build -f <name-path list> -k <kmer_size>`, pinned k-mer size in `workflow.json`), then compared with `ska distance`, which reports the number of SNPs differing between every pair — the *entire* pairwise matrix (every representative against every other, not just against the query).

Historical measurement with the former fixed 100/100 caps against a previously tested 714-representative Listeria cluster: ~42s total (versus ~14s at the old 5/20 caps), ~299MB downloaded (versus ~15MB). This timing does not benchmark the new adaptive policy or its 200-reference ceiling.

## Interpretation

`interpret_snp_resolution.py` reports qualifying query-to-reference comparisons,
excluded low-overlap comparisons, all tied nearest references, per-cluster
summaries, and an exploratory Neighbor-Joining tree when a complete comparable
matrix is available. Status is `COMPARED`, `AMBIGUOUS`, or `INSUFFICIENT_DATA`;
none establishes strain identity. The structural `confidence` block contains
raw distances and ratios, not statistical confidence. Zero observed differences
are not described as an exact whole-genome match. See
[cluster-expansion.md](cluster-expansion.md) for comparability settings and ties.

`render_snp_tree.py` renders `newick_tree` to `snp_resolution/tree.png` with `QUERY` highlighted in red and bold (via `phytreeviz`, already pulled in transitively by the pinned mashpit commit's own dependencies — no new container pin needed). Rendering is best-effort: a failure only sets `tree_image.status` to `FAIL` with the error message, it never fails SNP resolution or the underlying Mash screen. `generate_report.py` embeds this image in `report.md` when available, falling back to the raw Newick text otherwise.

`generate_report.py` renders all of this into `report.md`, written alongside `result.json` on every run (not just `--snp-resolve` ones): organism determination, the Mashpit candidate table, the ska2 cluster/genome tables and confidence statement, and the Newick tree, in plain language for a reader who doesn't want to open JSON.

## Limitations

- ska2 pairwise SNP distance is itself a screening refinement, not a validated outbreak-confirmation pipeline (unlike, e.g., CFSAN SNP Pipeline or an accredited cgMLST scheme). Always recommend a validated, organism-appropriate high-resolution comparison for actual outbreak confirmation, exactly as for the Mash result.
- `--snp-resolve` compares selected representatives. `--snp-expand` can add available members of plausible clusters from the exact release, subject to explicit budgets and stopping rules. Neither mode searches all public genomes or guarantees the nearest isolate globally.
- A failed or partial download (fewer than two usable genomes, including the query) skips SNP resolution with a warning; it does not fail the underlying Mash screen, which remains authoritative on its own.
