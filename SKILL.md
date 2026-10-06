---
name: isoscout
description: Reproducibly screen bacterial isolate assemblies or paired-end Illumina FASTQ/FASTQ.GZ reads with Mashpit. Use when Codex needs to validate isolate inputs, run fixed read QC and assembly, assess bacterial assembly quality, select a Salmonella, Escherichia coli/Shigella, Listeria, Campylobacter, or Cronobacter database from a user-provided organism or local MLST classification, run Mashpit, interpret structured candidate-cluster results, or explain conservative failure and setup states.
---

# IsoScout

Use the bundled scripts for computation. Do not construct bioinformatics commands, select alternate tools, change thresholds, parse raw output mentally, or infer an organism from sequence text.

## Run a screen

1. Read [references/setup.md](references/setup.md) when checking the environment or resolving a missing database.
2. Run:

   ```bash
     python3 scripts/screen_isolate.py INPUT [INPUT_R2] --output OUTPUT_DIR \
     [--database-root DATABASE_ROOT] [--organism ORGANISM_KEY] [--snp-resolve | --snp-expand] [--snp-selection-mode adaptive|all_returned] [--query-accession ACCESSION] \
     [--refinement-policy POLICY_JSON] [--snp-backend minimap2|mummer|ska] \
     [--assembly-manifest ASSEMBLIES_TSV] [--comparison-cache CACHE_DIR] [--comparison-workers 4]
   ```

3. Read `OUTPUT_DIR/result.json`. Use `status`, `stop_reason`, and `user_summary` as the authoritative result.
4. Report the concise summary first. State that a candidate is a screening result, not proof of outbreak relatedness. Recommend a validated high-resolution SNP comparison when a candidate is present.
5. Link the user to `result.json`, `provenance.json`, `report.md`, and retained logs. Never invent a result if a file is missing or a stage failed.

`OUTPUT_DIR/report.md` is a plain-language rendering of `result.json` for a non-technical reader — organism determination, Mashpit's candidate clusters and scores, and (when `--snp-resolve` ran) the target-to-candidate SNP tables, alignment coverage and interpretation. It is written on every run, not only `--snp-resolve` ones. Prefer linking to it over reciting JSON when the audience isn't reading code.

Pass `--snp-resolve` to compare the target assembly directly with selected Mashpit representatives using minimap2 and its upstream paftools variant caller once a candidate is found (read [references/snp-resolution.md](references/snp-resolution.md) first). This is opt-in because, unlike the rest of the screen, it may download candidate genomes from NCBI; `--assembly-manifest` supplies checksum-verified local assemblies first. Report `result.json`'s `snp_resolution` block alongside the Mash result when present, including whether it agrees with Mashpit's top candidate.

With `--snp-expand`, additionally explore members of selected returned clusters from the exact NCBI release under the versioned round/resource policy. Read [references/cluster-expansion.md](references/cluster-expansion.md). Report cluster-label concordance separately from the nearest examined reference set (all ties), per-cluster attempted/downloaded/qualifying/excluded/unexamined coverage, allocation choices, unresolved clusters, stop reason, and failed downloads. The versioned focused scheduler reserves rotating alternative slots and requires focused qualifying comparisons before sampled stability. `--refinement-policy` selects an explicit policy file; do not silently raise default budgets. Never equate correct cluster concordance, sampled stability, or zero observed SNPs with global nearest-genome recovery or strain identity. Assembly refinement compares only target–candidate pairs, reuses completed alignments, and does not construct a tree. Paired-read screens retain the legacy SKA2 path using cleaned reads; initial Mashpit screening uses the generated assembly. Use `--snp-backend mummer` only when that explicit alternative is desired, and never silently change backends on failure.

When available, use the database package's verified membership snapshot for
expansion. Legacy packages may retrieve the exact-release tables; a corrupt
snapshot must stop the analysis rather than be replaced by a newer release.

The script accepts exactly one assembly (`.fa`, `.fasta`, `.fna`) or one recognized R1/R2 FASTQ pair. The FASTQ path uses the fixed workflow described in [references/workflow.md](references/workflow.md). Long reads, hybrid reads, interleaved reads, and metagenomes are unsupported.

Use `--organism` when the organism is already known. It accepts `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, or `cronobacter` and directly selects that database. When omitted, the fixed workflow runs local `mlst --full --csv` against its pinned bundled PubMLST schemes and maps the detected scheme to a supported database. It does not upload sequence data.

Every successfully parsed Mashpit query also writes `similarity_distribution/summary.json` and, when plotting is available, `similarity_distribution/rank_similarity.png`. Report any return-limit or omitted-near-tie warnings from `result.json`. These diagnostics preview policy 3.0.1 (global first 50, up to 100 returned references in adaptive mode, 200 total reference attempts). Sketch tolerance is descriptive and never excludes a cluster from SNP consideration. With `--snp-resolve`, the same selector supplies the download list and records decisions in `snp_resolution/targets.json`. Report unexamined returned clusters and members, resource omissions, failed downloads, comparability exclusions, and all nearest ties. Distinguish Mashpit-displayed alternatives from SNP-tested alternatives. `all_returned` requires sufficient recorded budget and does not silently select a partial set. Read [references/similarity-distribution.md](references/similarity-distribution.md) to interpret the score tolerance and unknown boundaries.

For local assembly manifests, cache identity and candidate comparisons, read [references/assembly-comparison.md](references/assembly-comparison.md). Assembly refinement compares finalists on one common set of target positions, rechecks all outside candidates after that set changes, and retains conflicting evidence from larger aligned regions. Reported target deletion gaps and ambiguous differences are removed from the comparison regions. Partial-overlap SNP counts are lower bounds, never complete distances. Coverage below 85% is flagged, not excluded. Report `cluster_status` and `genome_status` separately; a cluster witness can support a cluster despite uncertainty between genomes within it. Preserve ties, failed comparisons, insufficient overlap and cross-cluster ambiguity. `ranked` contains display groups, not a total SNP-distance ordering; `nearest_snp_distance` is null. Inspect `finalist_panel`, `candidate_challenges` (including `region_sensitivity`), `cluster_certificate` and `decision_reference_samples` for evidence. Common-panel verification alone does not establish cluster or genome resolution. Conclusions apply only to the examined pool: expansion within selected clusters cannot recover an unseen cluster. Missing comparisons are never zero-SNP matches.

## Apply stop rules

Stop before Mashpit when input validation or broad assembly QC is `FAIL`, local MLST routing is unsupported or uncertain, the selected database is absent or invalid, or an external command fails. Proceed with visible caveats when QC or MLST routing is `WARN`. Query only the selected Mashpit database and never substitute another database.

Read references only as needed:

- [references/qc-policy.md](references/qc-policy.md): deterministic PASS/WARN/FAIL rules and evidence
- [references/database-routing.md](references/database-routing.md): user/local-MLST routing contract
- [references/mashpit-interpretation.md](references/mashpit-interpretation.md): result categories and wording
- [references/snp-resolution.md](references/snp-resolution.md): assembly or read SNP-distance refinement of a Mashpit candidate
- [references/limitations.md](references/limitations.md): scope and scientific limitations

## Preserve reproducibility

Keep the generated output directory intact. It records commands, logs, checksums, software versions, database metadata, workflow/profile versions, and timestamps. Do not edit result or provenance files after the run.
