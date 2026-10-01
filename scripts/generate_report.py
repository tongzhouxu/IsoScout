#!/usr/bin/env python3
"""Render result.json into a plain-language Markdown report.

Every fact in the report is a direct read of an already-computed field in
result.json - no new computation, no invented labels. Intended for someone
who does not want to read JSON: what organism was this screened against, what
did Mashpit find, what did ska2's SNP distances say, and what should happen
next.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any


TOP_N_GENOMES = 10

INPUT_TYPE_LABELS = {
    "assembly": "an existing genome assembly",
    "illumina_paired_fastq": "raw paired-end Illumina reads",
}


def _mashpit_section(result: dict[str, Any]) -> list[str]:
    lines = ["## Step 2: Screened against the Mashpit database", ""]
    database = result.get("database")
    if database:
        lines.append(
            f"Database: **{database['name']}** (NCBI Pathogen Detection release `{database['version']}`)."
        )
    match = result.get("mashpit_result") or {}
    best = match.get("best_candidate")
    if not best:
        lines.append("")
        lines.append("Mashpit found no candidate cluster for this query at all.")
        return lines
    lines.append("")
    if match.get("below_threshold"):
        lines.append(
            f"**Mashpit's best hit was only score {best['score']:.3f}, below its own screening "
            "threshold (0.85). This is not a meaningful match — treat it as noise, not a real candidate.**"
        )
        return lines
    lines.append(f"Best-matching cluster: **{best['cluster']}**, similarity score **{best['score']:.3f}** "
                 "(1.000 = identical sketch, 0 = no shared content).")
    alternatives = match.get("alternative_candidates") or []
    if alternatives:
        lines.append("")
        lines.append("Other clusters that scored close enough to also consider:")
        lines.append("")
        lines.append("| Cluster | Score |")
        lines.append("|---|---|")
        for alt in alternatives:
            lines.append(f"| {alt['cluster']} | {alt['score']:.3f} |")
    if match.get("ambiguous"):
        lines.append("")
        lines.append(
            "**Mashpit could not cleanly separate the top clusters** - more than one scored close enough "
            "to be a plausible match at this coarse screening resolution. This is exactly what SNP-level "
            "resolution (below) is for."
        )
    tree_image = match.get("tree_image") or {}
    lines.append("")
    if tree_image.get("status") == "PASS":
        lines.append(
            "**Mash tree** (coarse, sketch-based resolution - the query alongside its closest database "
            "representatives above Mashpit's 0.85 similarity threshold; built by Mashpit itself during "
            "the query, not by this skill):"
        )
        lines.append("")
        lines.append("![Mashpit tree](mashpit_tree.png)")
    else:
        lines.append(f"(No Mash tree available: {tree_image.get('reason', 'unknown reason')})")
    return lines


def _snp_section(result: dict[str, Any]) -> list[str]:
    snp = result.get("snp_resolution")
    lines = ["## Step 3: Confirmed with SNP-level detail (ska2)", ""]
    if not snp:
        lines.append("Not run for this screen (pass `--snp-resolve` to enable it).")
        return lines
    status = snp.get("status")
    if status == "SKIPPED":
        lines.append(f"Skipped: {snp.get('reason', 'no reason recorded')}")
        return lines
    if status == "ERROR":
        lines.append(f"Failed: {snp.get('error', 'unknown error')}")
        return lines
    if status != "RESOLVED":
        lines.append(f"Status: {status}")
        return lines

    lines.append(
        "Mashpit's similarity score is a coarse, sketch-based estimate. To get an exact count of genetic "
        "differences, the query and the representative genomes from the candidate cluster(s) above were "
        "compared base-by-base with ska2."
    )
    lines.append("")
    lines.append("**By cluster, closest to farthest:**")
    lines.append("")
    lines.append("| Cluster | Genomes compared | Closest match (SNPs) | Typical distance (SNPs) | Farthest (SNPs) |")
    lines.append("|---|---|---|---|---|")
    for row in snp.get("cluster_summary", []):
        lines.append(
            f"| {row['cluster']} | {row['genomes_compared']} | {row['min_snp_distance']:.2f} | "
            f"{row['median_snp_distance']:.2f} | {row['max_snp_distance']:.2f} |"
        )

    ranked = snp.get("ranked", [])
    lines.append("")
    lines.append(f"**{min(TOP_N_GENOMES, len(ranked))} closest individual genomes** "
                 f"(out of {len(ranked)} compared):")
    lines.append("")
    lines.append("| Genome | Cluster | SNP distance |")
    lines.append("|---|---|---|")
    for item in ranked[:TOP_N_GENOMES]:
        lines.append(f"| {item['sample']} | {item['cluster']} | {item['snp_distance']:.2f} |")
    if len(ranked) > TOP_N_GENOMES:
        lines.append("")
        lines.append(
            f"({len(ranked) - TOP_N_GENOMES} more in `snp_resolution/interpretation.json` under `ranked`.)"
        )

    confidence = snp.get("confidence", {})
    lines.append("")
    lines.append(f"**Bottom line:** {confidence.get('statement', 'No confidence statement available.')}")

    lines.append("")
    lines.append(
        "**SNP tree** (exact base-by-base resolution, built from the ska2 distances above, not from "
        "Mash sketches - your query highlighted):"
    )
    lines.append("")
    tree_image = snp.get("tree_image") or {}
    if tree_image.get("status") == "PASS":
        lines.append("![SNP tree with query highlighted](snp_resolution/tree.png)")
    else:
        lines.append(
            "(Image rendering "
            + ("was skipped: " + tree_image.get("reason", "") if tree_image.get("status") == "SKIPPED"
               else "failed: " + tree_image.get("error", "unknown error"))
            + " - the Newick text below can still be pasted into a tree viewer such as "
            "[iTOL](https://itol.embl.de/) or FigTree.)"
        )
    lines.append("")
    lines.append(
        "Newick text (also saved as `snp_resolution/interpretation.json`'s `newick_tree` field):"
    )
    lines.append("")
    lines.append("```")
    lines.append(snp.get("newick_tree") or "(tree unavailable)")
    lines.append("```")
    return lines


def generate_report(result: dict[str, Any]) -> str:
    lines = [
        f"# IsoScout screening report: {result['sample']}",
        "",
        f"**Overall status: {result['status']}**",
        "",
        "## Step 1: Sample and organism",
        "",
        f"Input: {INPUT_TYPE_LABELS.get(result.get('input_type'), result.get('input_type', 'unknown'))}.",
    ]
    read_qc = result.get("read_qc")
    if read_qc:
        coverage = read_qc.get("estimated_coverage")
        q30 = read_qc.get("q30_fraction")
        lines.append(
            f"Read quality: **{read_qc['status']}** — "
            + (f"{coverage:.0f}x estimated coverage" if coverage is not None else "coverage unavailable")
            + (f", {q30:.1%} of bases at Q30 or better" if q30 is not None else "")
            + f" ({read_qc.get('read_pairs_validated', 'unknown')} read pairs)."
        )
    if result.get("assembly_qc"):
        qc = result["assembly_qc"]
        lines.append(f"Assembly quality check: **{qc['status']}**"
                     + (f" ({'; '.join(qc['warnings'])})" if qc.get("warnings") else "") + ".")
    routing = result.get("routing") or {}
    if routing.get("source") == "user":
        lines.append(f"Organism: **{routing.get('organism')}** (told to us directly, not auto-detected).")
    elif routing.get("source") == "mlst":
        mlst = routing.get("mlst", {})
        lines.append(
            f"Organism: **{routing.get('organism')}** (auto-detected from the assembly's DNA using the "
            f"`{mlst.get('scheme')}` typing scheme, confidence status `{mlst.get('status')}`)."
        )
    if result.get("stop_reason"):
        lines.append("")
        lines.append(f"**Screening stopped: {result['stop_reason']}**")
        return "\n".join(lines) + "\n"
    lines.append("")
    lines.extend(_mashpit_section(result))
    lines.append("")
    lines.extend(_snp_section(result))
    lines.append("")
    lines.append("## What this does and does not mean")
    lines.append("")
    lines.append(
        "A close genetic match means this sample and the matched genome(s) likely share a recent common "
        "source. It does **not** by itself prove they are from the same outbreak, food source, or event - "
        "that requires a public health investigation using this result as a starting point, plus a "
        "validated, accredited confirmation pipeline (e.g. an organism-appropriate reference-based SNP "
        "pipeline or cgMLST scheme)."
    )
    warnings = result.get("warnings", [])
    if warnings:
        lines.append("")
        lines.append("**Warnings raised during this screen:**")
        for warning in warnings:
            lines.append(f"- {warning}")
    return "\n".join(lines) + "\n"


def main() -> int:
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("--result", required=True, help="Path to result.json")
    parser.add_argument("--output", required=True, help="Path to write report.md")
    args = parser.parse_args()
    result = json.loads(Path(args.result).read_text(encoding="utf-8"))
    Path(args.output).write_text(generate_report(result), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
