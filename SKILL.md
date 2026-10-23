---
name: isoscout
description: Screen a bacterial isolate assembly or paired-end Illumina reads against supported NCBI Pathogen Detection cluster databases with IsoScout. Use for local, reproducible Mashpit screening of Salmonella, E. coli/Shigella, Listeria, Campylobacter, or Cronobacter, including input validation, local MLST routing, result interpretation, and optional candidate refinement.
---

# IsoScout

Run the bundled workflow; do not assemble bioinformatics commands or choose new thresholds. IsoScout is a screening tool. A candidate cluster is evidence for follow-up, not proof of outbreak relatedness.

## Screen an isolate

1. Confirm Docker is running and the current image was built from this checkout. Read [references/setup.md](references/setup.md) only when setup or a database is missing.
2. Confirm the user supplied either one FASTA/FNA assembly or an R1/R2 Illumina FASTQ pair. Long reads, hybrid reads, interleaved reads, and metagenomes are unsupported.
3. Run the host launcher from this skill directory:

   ```bash
   python3 <skill-directory>/scripts/launch_isoscout.py INPUT [INPUT_R2] \
     --output OUTPUT_DIR [--organism ORGANISM_KEY]
   ```

   Use a new output directory. `ORGANISM_KEY` is one of `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, or `cronobacter`. Omit it only when the organism is unknown; local MLST then routes among those five supported databases. This is not general species identification. The routed database must already be installed.
4. Read `OUTPUT_DIR/result.json`. Treat `status`, `stop_reason`, and `user_summary` as authoritative. Never infer a result from incomplete logs.
5. Give the concise result first. Explain the evidence, alternatives or ties, what references were examined or left unexamined, material limitations, and the next appropriate validation step. Link `report.md`, `result.json`, and `provenance.json`.

## CDC REP designations

Databases from the [`databases-rep-2026-10`](https://github.com/tongzhouxu/IsoScout/releases/tag/databases-rep-2026-10) release carry a `REP_CLUSTER` table. When it is present, `result.json` field `rep_designations` and the report give the REP designation and source URL for the final called cluster (the SNP-resolved cluster, or Mashpit's top candidate when refinement did not run) and for Mashpit's displayed candidates, matched by exact PDS accession. Report a labelled cluster as "the called NCBI SNP cluster is linked to REP strain X", never as confirmation that the isolate belongs to the outbreak strain. Report a missing table or row (`NO_LABEL`, `TABLE_ABSENT`) as "no REP label in this database", never as "not a REP strain". `UNREADABLE` means no REP status is known.

The default screen uses Mashpit only and does not download candidate genomes. It writes the plain-language `report.md` on every run. It also writes score-distribution diagnostics when Mashpit returns a valid query result; report any return-limit or omitted-near-tie warnings.

## Optional refinement

Use `--snp-resolve` only when the user requests candidate refinement or authorizes its NCBI genome downloads. It compares the target with selected returned representatives. Read [references/snp-resolution.md](references/snp-resolution.md) before interpreting it, preserve all ties and failures, and state whether the refined evidence agrees with the Mashpit candidate.

Use `--snp-expand` only when the user explicitly requests bounded exploration of additional members from the selected clusters. Read [references/cluster-expansion.md](references/cluster-expansion.md). Report examined and unexamined coverage, allocation and stop reasons, unavailable references, and all equally supported results. Sampled stability, cluster concordance, or zero observed SNPs does not establish the globally nearest genome or strain identity.

For checksum-verified local candidate assemblies, comparison caches, or alternate configured backends, read [references/assembly-comparison.md](references/assembly-comparison.md). Do not silently change a backend, policy, or resource budget after a failure. Assembly comparisons use target-relative regions that can differ between references, so raw pairwise SNP counts may be lower bounds and must not be presented as a global ranking. Report cluster and genome resolution separately.

## Stop conditions

Stop and report the recorded reason when input validation or broad assembly QC fails, MLST routing is unsupported or uncertain, the selected database is absent or invalid, or an external command fails. Continue with visible caveats for `WARN`. Never substitute another organism database.

Read these only as needed:

- [references/database-routing.md](references/database-routing.md) for known-organism and local-MLST routing
- [references/qc-policy.md](references/qc-policy.md) for fixed PASS/WARN/FAIL rules
- [references/mashpit-interpretation.md](references/mashpit-interpretation.md) for candidate wording
- [references/similarity-distribution.md](references/similarity-distribution.md) for score and retrieval diagnostics
- [references/limitations.md](references/limitations.md) for scientific scope

Keep the output directory intact. It contains input and database checksums, actual source and container identities when launched through the host launcher, observed tool versions, commands, logs, and policy versions. Do not edit generated result or provenance files.
