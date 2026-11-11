# SNP resolution: selection policy 3.0.1

Default SNP comparison uses target-to-candidate assembly alignments (minimap2/paftools by default). Paired-read input uses the QC-checked SKESA assembly; supplied assemblies enter the same comparison stage. See [assembly-comparison.md](assembly-comparison.md) for the assembly method, local inputs, caching and criteria. Governed by `config/snp-resolution-policy.json` (versioned resource policy) and the `snp_resolution` block of `config/workflow.json` (default assembly backend). Runs automatically; use `--screen-only` explicitly to omit SNP comparison.

## Why this exists

Mashpit's MinHash similarity is a coarse, sketch-resolution screen. A single unambiguous top cluster at Mash resolution is not proof of SNP-level closeness, and two near-tied clusters cannot be told apart by Mash alone. SNP resolution compares the target with candidates using the selected backend. The distance definition and alignment coverage must accompany any closest-genome interpretation.

## When it runs

By default, unless `--screen-only` is requested, when Mashpit returned an eligible candidate (`mashpit_result.best_candidate` is present) — regardless of whether that candidate was unambiguous. This is deliberate: an unambiguous Mash hit can still be many SNPs away from its representatives, so the refinement step is not limited to the ambiguous case.

## Network dependency

This step can use checksum-verified local assemblies through `--assembly-manifest`. Missing candidates require public downloads. A Mashpit database only retains sourmash signatures, not the representative assemblies (they are sketched and discarded during `mashpit build`), so resolving SNPs requires re-downloading selected representative genomes from NCBI via the pinned `datasets` CLI. The full workflow may therefore contact NCBI for public reference genomes. The query sequence itself is never uploaded. Explicit `--screen-only` mode performs no candidate downloads.

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

The default full analysis uses adaptive selection. `--snp-expand` additionally
explores exact-release members of clusters represented in the selected returned
set. The expansion profile is versioned separately from representative selection;
`--refinement-policy` can select another explicit, versioned policy file without
changing Mashpit retrieval or the initial candidate set. Its expansion settings are recorded. Assembly comparability comes from the separately recorded assembly policy. No score-distribution
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

Local manifest entries are checked first; only missing accessions are downloaded
with the pinned NCBI datasets CLI. Assemblies default to minimap2/paftools
comparisons against the target alone, with new expansion candidates processed
incrementally. See [assembly-comparison.md](assembly-comparison.md).

Paired-read runs compare the generated assembly using the same backend and
criteria. `--snp-backend mummer` selects the alternative assembly backend for
either input mode. There is no automatic backend fallback.

## Interpretation

Assembly results report qualifying target–candidate comparisons, all ties,
alignment coverage, per-pair failures, cache use, common-finalist counts, outside-reference lower bounds, regional contradictions and unresolved alternatives. Read [assembly-comparison.md](assembly-comparison.md) for the common-region comparisons and insufficient-evidence rules.
Target-to-candidate comparisons do not generate a phylogenetic tree. Results
retain unresolved alternatives and excluded comparisons; zero observed SNPs
does not establish whole-genome identity.

The final `coverage` block lists every returned representative left unexamined,
every attempted, verified, and unavailable download, every qualifying and
comparability-excluded query comparison, and discovered but unexamined members.
It also flags when member coverage is unknown because membership was not read.
Per-round allocation and newly qualifying comparisons are in `expansion.json`;
the final report separates cluster-label concordance from recovery of an
individual nearest genome among those examined.
The report separates Mashpit-displayed alternatives from alternatives actually
tested by the selected backend. No stable within-cluster result rules out an untested cluster.

`generate_report.py` renders the results into `report.md`, written alongside
`result.json` on every run (including `--screen-only` runs): organism determination,
the Mashpit candidate table, the selected backend's cluster/genome results,
alignment coverage, comparison plots and uncertainty. The SNP `query_input` and
file checksums identify the actual assembly compared; the top-level `input_type`
retains whether the original input was an assembly or paired reads.

## Limitations

- SNP comparison with minimap2/paftools or MUMmer4 is a screening refinement, not a validated outbreak-confirmation pipeline (unlike, e.g., CFSAN SNP Pipeline or an accredited cgMLST scheme). Always recommend a validated, organism-appropriate high-resolution comparison for actual outbreak confirmation, exactly as for the Mash result.
- The default full analysis compares selected representatives. `--snp-expand` can add available members of selected clusters from the exact release, subject to explicit budgets and stopping rules. Neither mode searches all public genomes or guarantees the nearest isolate globally.
- Missing required comparisons make a full run `INCOMPLETE` with a nonzero exit status. Completed comparisons and the Mashpit screen remain available as partial evidence.
