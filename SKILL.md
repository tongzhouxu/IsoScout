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
     [--database-root DATABASE_ROOT] [--organism ORGANISM_KEY] [--snp-resolve]
   ```

3. Read `OUTPUT_DIR/result.json`. Use `status`, `stop_reason`, and `user_summary` as the authoritative result.
4. Report the concise summary first. State that a candidate is a screening result, not proof of outbreak relatedness. Recommend a validated high-resolution SNP comparison when a candidate is present.
5. Link the user to `result.json`, `provenance.json`, `report.md`, and retained logs. Never invent a result if a file is missing or a stage failed.

`OUTPUT_DIR/report.md` is a plain-language rendering of `result.json` for a non-technical reader — organism determination, Mashpit's candidate clusters and scores, and (when `--snp-resolve` ran) the ska2 SNP distance tables, confidence statement, and Newick tree. It is written on every run, not only `--snp-resolve` ones. Prefer linking to it over reciting JSON when the audience isn't reading code.

Pass `--snp-resolve` to additionally compute ska2 pairwise SNP distances between the query and the relevant Mashpit representatives once a candidate is found (read [references/snp-resolution.md](references/snp-resolution.md) first). This is opt-in because, unlike the rest of the screen, it downloads representative genomes from NCBI. Report `result.json`'s `snp_resolution` block alongside the Mash result when present, including whether it agrees with Mashpit's top candidate.

The script accepts exactly one assembly (`.fa`, `.fasta`, `.fna`) or one recognized R1/R2 FASTQ pair. The FASTQ path uses the fixed workflow described in [references/workflow.md](references/workflow.md). Long reads, hybrid reads, interleaved reads, and metagenomes are unsupported.

Use `--organism` when the organism is already known. It accepts `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, or `cronobacter` and directly selects that database. When omitted, the fixed workflow runs local `mlst --full --csv` against its pinned bundled PubMLST schemes and maps the detected scheme to a supported database. It does not upload sequence data.

## Apply stop rules

Stop before Mashpit when input validation or broad assembly QC is `FAIL`, local MLST routing is unsupported or uncertain, the selected database is absent or invalid, or an external command fails. Proceed with visible caveats when QC or MLST routing is `WARN`. Query only the selected Mashpit database and never substitute another database.

Read references only as needed:

- [references/qc-policy.md](references/qc-policy.md): deterministic PASS/WARN/FAIL rules and evidence
- [references/database-routing.md](references/database-routing.md): user/local-MLST routing contract
- [references/mashpit-interpretation.md](references/mashpit-interpretation.md): result categories and wording
- [references/snp-resolution.md](references/snp-resolution.md): optional ska2 SNP-distance refinement of a Mashpit candidate
- [references/limitations.md](references/limitations.md): scope and scientific limitations

## Preserve reproducibility

Keep the generated output directory intact. It records commands, logs, checksums, software versions, database metadata, workflow/profile versions, and timestamps. Do not edit result or provenance files after the run.
