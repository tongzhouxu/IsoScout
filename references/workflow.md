# Fixed workflow

`scripts/screen_isolate.py` owns routing; subordinate scripts own computation.

Assembly input:

```text
validate FASTA → broad assembly QC → user organism or local MLST scheme detection
→ select one supported database → organism-specific QC → run Mashpit once
→ structured Mashpit parser → provenance
```

Paired Illumina input adds:

```text
validate complete R1/R2 FASTQ → fastp → SKESA → QUAST → assembly-input path
```

Commands and parameters come from versioned `config/*.json` policies. Mashpit is pinned to upstream commit `538d3421302fe6dd129780605b8ff5dedbf4c046c`. The fixed profile is:

```text
mashpit query ASSEMBLY DATABASE --number 200 --threshold 0.85 --tie-tolerance-hashes 2
```

Mashpit writes its selected-database query within the retained `mashpit/` directory. The parser consumes the pinned commit's native `<sample>_representative_matches.csv` and `<sample>_cluster_candidates.csv` files. Missing, multiple, or malformed result files are errors.

The assembly component returns:

```json
{
  "status": "PASS",
  "assembly_path": "assembly.fasta",
  "read_paths": ["trimmed_R1.fastq.gz", "trimmed_R2.fastq.gz"],
  "read_qc": {},
  "commands": [],
  "workflow_version": "1.1.0"
}
```

Any future replacement assembly Skill must satisfy this contract and must not weaken input validation or QC stop rules.

`--organism` directly selects a supported database when the organism is already known. Without it, local `mlst` 2.35.0 auto-detects a bundled PubMLST scheme. The classifier's complete CSV, stderr log, normalized routing JSON, and command are retained.

## Similarity diagnostics

After parsing Mashpit, the workflow runs `analyze_similarity_distribution.py` to describe the returned scores and audit the versioned SNP selector. It writes a JSON summary and a rank–similarity plot, included in the report even without `--snp-resolve`. Errors or possible truncation add visible warnings without stopping the screening/refinement workflow. See [similarity-distribution.md](similarity-distribution.md). The current workflow uses policy 3.0.1: the global first 50 are guaranteed when available, with up to 100 returned representatives in adaptive mode and 200 total reference attempts. Query limits and the screening gate are unchanged.

## Optional SNP resolution

With `--snp-resolve`, once Mashpit returns a candidate, an additional stage selects selected returned representatives across clusters, downloads them from NCBI, and compares assemblies directly with the target using minimap2/paftools, or uses the legacy SKA path for reads. See [references/snp-resolution.md](snp-resolution.md). This stage is opt-in and its failure only downgrades `result.json` to a warning — it never turns an otherwise-successful Mash screen into a stopped one.

With `--snp-expand` (workflow 1.9.0), refine in evidence-led focused rounds over the exact-release member pool with reserved rotating alternative coverage. Cleaned read pairs feed SKA2 directly when the input was reads. See [cluster-expansion.md](cluster-expansion.md) for resource limits, feedback, and stop states.
