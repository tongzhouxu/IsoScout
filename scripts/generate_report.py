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
        lines.append("Other Mashpit-displayed clusters (sketch alternatives; SNP testing coverage is reported separately):")
        lines.append("")
        lines.append("| Cluster | Score |")
        lines.append("|---|---|")
        for alt in alternatives:
            lines.append(f"| {alt['cluster']} | {alt['score']:.3f} |")
    if match.get("ambiguous"):
        lines.append("")
        lines.append(
            "**Mashpit could not cleanly separate the top clusters** - more than one scored close enough "
            "to be a plausible match at this coarse screening resolution. The SNP comparison below only covers "
            "references actually downloaded and compared."
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


def _similarity_section(result: dict[str, Any]) -> list[str]:
    diagnostics = result.get("similarity_distribution")
    if not diagnostics:
        return []
    lines = ["### Candidate similarity distribution", ""]
    if diagnostics["status"] == "ERROR":
        return lines + [f"Diagnostics unavailable: {diagnostics['error']}", ""]
    lines.append(
        f"Returned **{diagnostics['returned_genomes']} unique representatives**; "
        f"**{diagnostics['within_top_tolerance_count']}** are within Mashpit's sketch tolerance "
        f"of the best score (tolerance {diagnostics['tolerance']['score_tolerance']:.6g}). "
        "This describes the returned subset, not all isolates in the database."
    )
    lines.extend(["", "| Rank | Score | Gap to next | Cutoff splits a near-tie? |", "|---|---|---|---|"])
    for point in diagnostics["rank_checkpoints"]:
        score = f"{point['score']:.6f}" if point["score"] is not None else "Not returned"
        gap = f"{point['gap_to_next']:.6g}" if point["gap_to_next"] is not None else "Unknown"
        tied = "Unknown" if point["splits_near_tie"] is None else ("Yes" if point["splits_near_tie"] else "No")
        lines.append(f"| {point['rank']} | {score} | {gap} | {tied} |")
    preview = diagnostics["selection_preview"]
    lines.append("")
    if preview["eligible"]:
        lines.append(
            f"The current SNP selection policy would select **{preview['selected_genomes']}** unique genomes "
            f"and omit **{preview['omitted_within_top_tolerance_count']}** returned genomes within tolerance "
            "of the best score. This is a policy preview, not a count of downloaded or compared genomes."
        )
        audit = preview.get("audit", {})
        if audit.get("policy", {}).get("strategy") == "global-ranked-cluster-coverage-v1":
            lines.extend(["", (
                f"Mode **{audit['mode']}**: first {audit['policy']['initial_target_genomes']} global ranked references, "
                f"up to {audit['policy']['max_returned_representatives']} returned references in adaptive mode, "
                f"and {audit['policy']['max_total_genomes']} total reference attempts. "
                f"**{audit['unexamined_returned_representatives']}** returned representatives and "
                f"**{len(audit['unexamined_clusters'])}** whole returned clusters are outside this selection preview. "
                "These omissions reflect the budget, not SNP evidence."
            )])
    else:
        lines.append("SNP selection preview is inactive because no candidate meets the current screening gate.")
    lines.append("")
    if diagnostics["retrieval"]["limit_reached"]:
        lines.append("**Mashpit's return limit was reached.** Additional matches may be unobserved; their existence and scores are unknown.")
    else:
        lines.append("Mashpit's return limit was not reached.")
    lines.extend(["", "The tolerance is an experimental descriptive sketch-resolution heuristic, not a SNP eligibility gate, confidence interval, ANI estimate, or strain-assignment cutoff.", ""])
    if diagnostics.get("plot", {}).get("status") == "PASS":
        lines.extend(["![Rank versus Mashpit similarity](similarity_distribution/rank_similarity.png)", ""])
    else:
        lines.extend(["(Similarity plot unavailable; numerical diagnostics are retained.)", ""])
    lines.append("[Full ranked scores, cluster composition, and selection audit](similarity_distribution/summary.json)")
    return lines


def _snp_section(result: dict[str, Any]) -> list[str]:
    snp = result.get("snp_resolution")
    backend = (snp or {}).get("backend","ska")
    lines = [f"## Step 3: Compared at SNP resolution ({backend})", ""]
    if not snp:
        lines.append("Not run for this screen (pass `--snp-resolve` to enable it).")
        return lines
    expansion = snp.get("expansion", {})
    if expansion.get("enabled"):
        lines.extend([
            f"Cluster expansion: **{expansion['stop_reason']}**, {len(expansion['rounds'])} comparison round(s), "
            f"{expansion['attempted_genomes']} reference downloads attempted. "
            f"Additional attempts: {expansion.get('additional_reference_attempts', 0)}/{expansion['policy'].get('max_additional_reference_attempts', 'unknown')}. "
            f"Unexamined members in the available pool: {expansion['pending_available_members'] if expansion['pending_available_members'] is not None else 'unknown'}.",
            "", "[Expansion decisions and pinned membership provenance](snp_resolution/expansion.json)", "",
            "A stable sampled neighborhood does not establish that every closer genome has been found.", "",
        ])
        uncertainty = expansion.get("remaining_uncertainty") or {}
        if uncertainty:
            lines.append(
                f"Remaining uncertainty: {uncertainty.get('unexamined_additional_pool') if uncertainty.get('unexamined_additional_pool') is not None else 'unknown'} "
                f"additional references unexamined; {len(uncertainty.get('unresolved_selected_clusters') or [])} "
                "selected clusters without qualifying SNP evidence; "
                f"{len(uncertainty.get('leading_clusters_below_focused_minimum') or [])} leading clusters below the configured focused minimum. "
                f"Budget or round limit prevented required exploration: {'yes' if uncertainty.get('budget_or_round_limit_prevented_required_exploration') else 'no'}."
            )
            lines.append("")
        lines.extend(["Expansion decisions after each comparison:", "",
                      "| Round | Leading SNP clusters | Next leading requests | Next alternative requests | New qualifying comparisons | Decision |",
                      "|---:|---|---:|---:|---:|---|"])
        for round_record in expansion["rounds"]:
            allocation = round_record.get("allocation") or {}
            chosen = allocation.get("chosen") or []
            leading = sum(row.get("lane") in {"leading_snp_cluster", "unused_alternative_capacity_to_leader"} for row in chosen)
            alternatives = len(chosen) - leading
            feedback = round_record.get("feedback") or {}
            clusters = ", ".join(allocation.get("leading_clusters", [])) or "None recorded"
            lines.append(f"| {round_record['round']} | {clusters} | {leading} | {alternatives} | "
                         f"{len(feedback.get('new_qualifying_comparisons', []))} | {feedback.get('reason', 'unknown')} |")
        lines.extend(["", "The allocation audit records every requested accession, lane, cursor, and rationale in `snp_resolution/expansion.json`.", ""])
    query_input = snp.get("query_input", {})
    if query_input:
        lines.extend(["Target input: **" + ("cleaned paired-end reads" if query_input["type"] == "paired_reads" else "assembly") + "**.", ""])
    status = snp.get("status")
    if snp.get("selection") or snp.get("decisions"):
        lines.extend(["[Reference selection and inclusion/exclusion reasons](snp_resolution/targets.json)", ""])
    coverage = snp.get("coverage") or {}
    if coverage:
        lines.extend([
            f"SNP reference coverage: {coverage['mashpit_returned_unique']} unique Mashpit representatives returned; "
            f"{coverage['returned_selected']} selected; {len(coverage['reference_downloads_attempted'])} downloads attempted; "
            f"{len(coverage['reference_downloads_verified'])} verified; "
            f"{len(coverage['query_comparisons_qualifying'])} qualifying query comparisons; "
            f"{len(coverage['query_comparisons_excluded'])} comparability exclusions; "
            f"{len(coverage['returned_unexamined'])} returned representatives unexamined.",
            f"Unexamined available members: {len(coverage['members_unexamined']) if coverage['members_unexamined'] is not None else 'unknown'}. "
            "Mashpit-displayed alternatives are not necessarily SNP-tested alternatives.",
            "", "[Full reference and member coverage](snp_resolution/interpretation.json)", "",
        ])
        displayed = result.get("mashpit_result") or {}
        displayed_clusters = [row["cluster"] for row in ([displayed["best_candidate"]] if displayed.get("best_candidate") else []) + displayed.get("alternative_candidates", [])]
        if displayed_clusters:
            by_cluster = {row["cluster"]: row for row in coverage["clusters"]}
            lines.extend(["Mashpit-displayed clusters versus SNP-tested references:", "",
                          "| Cluster | Returned | Selected | Qualifying SNP comparisons | Unexamined returned |",
                          "|---|---:|---:|---:|---:|"])
            for cluster in dict.fromkeys(displayed_clusters):
                row = by_cluster.get(cluster, {})
                lines.append(f"| {cluster} | {row.get('returned', 0)} | {row.get('selected', 0)} | {row.get('qualifying_query_comparisons', 0)} | {row.get('returned_unexamined', 0)} |")
            lines.append("")
        lines.extend(["SNP comparison coverage by returned cluster (including any explored members):", "",
                      "| Cluster | Attempted | Downloaded | Qualifying | Excluded | Returned unexamined | Members unexamined |",
                      "|---|---:|---:|---:|---:|---:|---:|"])
        for row in coverage.get("clusters", []):
            members_left = row.get("members_unexamined")
            lines.append(f"| {row['cluster']} | {row['download_attempted']} | {row['download_verified']} | "
                         f"{row['qualifying_query_comparisons']} | {row['excluded_query_comparisons']} | "
                         f"{row['returned_unexamined']} | {members_left if members_left is not None else 'Unknown'} |")
        lines.append("")
    if status == "SKIPPED":
        lines.append(f"Skipped: {snp.get('reason', 'no reason recorded')}")
        return lines
    if status == "ERROR":
        lines.append(f"Failed: {snp.get('error', 'unknown error')}")
        return lines
    if status not in {"COMPARED", "AMBIGUOUS", "RESOLVED"}:
        lines.append(f"Status: {status}")
        return lines

    mash_best = (result.get("mashpit_result") or {}).get("best_candidate") or {}
    nearest_clusters = snp.get("nearest_clusters", [])
    agreement = snp.get("agrees_with_mash_top_candidate")
    agreement_text = "agrees" if agreement is True else ("disagrees" if agreement is False else "unresolved")
    lines.append(f"**Cluster-label result among examined qualifying references:** nearest SNP cluster(s): "
                 f"{', '.join(nearest_clusters) if nearest_clusters else 'none'}; Mashpit top cluster: "
                 f"{mash_best.get('cluster', 'unknown')}; comparison: {agreement_text}. "
                 "This is not agreement with an independent cluster label or proof of the globally nearest genome.")
    lines.append("")

    if backend in {"mummer","minimap2"}:
        lines.extend([
            f"Each candidate was aligned directly to the target assembly using {backend}. No candidate-to-candidate matrix was calculated.", "",
            snp.get("confidence",{}).get("statement","No qualifying comparisons."), "",
            "| Candidate | Cluster | SNPs | Indel bases | Target aligned (%) | Candidate aligned (%) |",
            "|---|---|---:|---:|---:|---:|",
        ])
        for row in snp.get("ranked",[])[:TOP_N_GENOMES]:
            lines.append(f"| {row['sample']} | {row['cluster']} | {row['snp_distance']} | {row['indel_bases']} | "
                         f"{100*row['target_aligned_fraction']:.2f} | {100*row['candidate_aligned_fraction']:.2f} |")
        lines.extend(["", "All minimum-SNP-count ties: " + ", ".join(snp.get("nearest_samples",[])),
                      "Minimum SNP rate per aligned target Mb: " + ", ".join(snp.get("nearest_by_aligned_snp_rate",[])), ""])
        if snp.get("ranking_basis_conflict"):
            lines.append("**Closest-genome conclusion unresolved:** SNP-count and aligned-SNP-rate minima disagree.")
        if (snp.get("ranking_image") or {}).get("status")=="PASS":
            lines.extend(["", "![Candidate SNP distances and alignment coverage](snp_resolution/candidate_distances.png)"])
        lines.extend(["", "Full pair-level metrics, failures, criteria and cache records: [comparison audit](snp_resolution/interpretation.json).",
                      "A target-only comparison does not define a phylogenetic tree."])
        return lines

    lines.append(
        "SKA2 counts differences in comparable split-kmer contexts. Ranking includes only references "
        "that pass the configured shared-kmer and missingness checks. These experimental comparability "
        "checks are not strain-assignment thresholds."
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
                 f"(out of {len(ranked)} qualifying references):")
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
    lines.append(f"**Interpretation:** {confidence.get('statement', 'No comparison statement available.')}")
    if snp.get("nearest_samples"):
        lines.append("All equally nearest qualifying references: " + ", ".join(snp["nearest_samples"]) + ".")
    lines.append("A stable result within sampled clusters does not rule out unexamined returned clusters or members.")

    lines.append("")
    lines.append(
        "**Exploratory SNP-distance tree** (built from qualifying SKA2 comparisons; your query highlighted):"
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
    lines.extend(_similarity_section(result))
    lines.append("")
    lines.extend(_snp_section(result))
    lines.append("")
    lines.append("## What this does and does not mean")
    lines.append("")
    lines.append(
        "A close match identifies candidate genomic neighbors within the examined reference set. "
        "It does **not** by itself establish a common source, outbreak, or event - "
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
