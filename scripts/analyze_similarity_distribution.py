#!/usr/bin/env python3
"""Describe returned Mashpit scores and audit the current selection policy.

Diagnostics only: never expand a query, download genomes, or change selection.
The tolerance follows generate_cluster_table in pinned Mashpit commit
538d3421302fe6dd129780605b8ff5dedbf4c046c: tie_tolerance_hashes / hash_number.
This is a sketch-resolution heuristic, not a confidence interval or SNP cutoff.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, WorkflowError, load_json, sha256_file, write_json
from parse_mashpit_results import (
    load_candidates,
    locate_candidate_file, locate_representative_file,
)
from select_snp_targets import select_targets
from mashpit_similarity import read_scores, within_tolerance, recorded_query_context


def analyze(
    mashpit_output_dir: Path, database: dict[str, Any], profile: dict[str, Any],
    selection_policy: dict[str, Any], rank_checkpoints: list[int],
    selection_mode: str | None = None, query_accession: str | None = None,
) -> dict[str, Any]:
    source = locate_representative_file(mashpit_output_dir)
    candidate_source = locate_candidate_file(mashpit_output_dir)
    rows, duplicate_rows = read_scores(source)
    hash_number = database.get("mashpit_database_settings", {}).get("hash_number")
    if not isinstance(hash_number, int) or isinstance(hash_number, bool) or hash_number <= 0:
        raise WorkflowError("A positive database hash_number is required; no sketch size is assumed.")
    limit = profile["number"]
    tolerance_hashes = profile["tie_tolerance_hashes"]
    if not isinstance(limit, int) or isinstance(limit, bool) or limit <= 0:
        raise WorkflowError("Mashpit return limit must be a positive integer.")
    if not isinstance(tolerance_hashes, (int, float)) or not math.isfinite(tolerance_hashes) or tolerance_hashes < 0:
        raise WorkflowError("Mashpit hash tolerance must be finite and nonnegative.")
    if any(not isinstance(rank, int) or isinstance(rank, bool) or rank <= 0 for rank in rank_checkpoints):
        raise WorkflowError("Rank checkpoints must be positive integers.")
    tolerance = tolerance_hashes / hash_number
    warnings = []
    if duplicate_rows:
        warnings.append(
            f"Collapsed {duplicate_rows} duplicate accession rows for diagnostics and selection."
        )
    candidates = load_candidates(candidate_source)
    eligible = bool(candidates) and candidates[0]["score"] >= profile["threshold"]
    selected = select_targets(mashpit_output_dir, selection_policy, database, profile,
                              mode=selection_mode, query_accession=query_accession)
    relevant = selected["selected_clusters"]
    warnings.extend(selected.get("warnings", []))
    selected_accessions = {item["accession"] for item in selected["targets"]}
    best = rows[0]["score"] if rows else None
    for index, row in enumerate(rows):
        row.update({
            "rank": index + 1,
            "gap_from_best": best - row["score"],
            "gap_to_next": row["score"] - rows[index + 1]["score"] if index + 1 < len(rows) else None,
            "within_top_tolerance": within_tolerance(best - row["score"], tolerance),
            "selected_by_current_policy": row["accession"] in selected_accessions,
        })
    clusters = []
    for cluster in sorted({row["cluster"] for row in rows}):
        members = [row for row in rows if row["cluster"] == cluster]
        kept = [row for row in members if row["selected_by_current_policy"]]
        omitted = [row for row in members if not row["selected_by_current_policy"]]
        gap = min(row["score"] for row in kept) - max(row["score"] for row in omitted) if kept and omitted else None
        clusters.append({
            "cluster": cluster, "returned_genomes": len(members),
            "best_score": members[0]["score"], "lowest_score": members[-1]["score"],
            "within_top_tolerance_count": sum(row["within_top_tolerance"] for row in members),
            "represented_in_selection_preview": cluster in relevant,
            "selected_genomes": len(kept), "omitted_genomes": len(omitted),
            "selected_to_omitted_gap": gap,
            "selection_splits_near_tie": within_tolerance(gap, tolerance) if gap is not None else False,
        })
    checkpoints = []
    for rank in sorted(set(rank_checkpoints)):
        observed = rank <= len(rows)
        gap = rows[rank - 1]["gap_to_next"] if observed else None
        checkpoints.append({
            "rank": rank, "observed": observed,
            "score": rows[rank - 1]["score"] if observed else None,
            "gap_to_next": gap,
            "splits_near_tie": within_tolerance(gap, tolerance) if gap is not None else None,
        })
    limit_reached = len(rows) >= limit
    tail_gap = rows[-2]["gap_to_next"] if len(rows) > 1 else None
    tail_near_tie = within_tolerance(tail_gap, tolerance) if tail_gap is not None else None
    top_band_reaches_tail = rows[-1]["within_top_tolerance"] if rows else None
    omitted_top = sum(row["within_top_tolerance"] and not row["selected_by_current_policy"] for row in rows)
    if limit_reached:
        warnings.append(
            "Mashpit's return limit was reached; additional database matches may be unobserved. "
            "This does not prove that additional matches exist."
        )
        if top_band_reaches_tail:
            warnings.append("The near-top score band reaches the return limit; its full size is unknown.")
        elif tail_near_tie:
            warnings.append("The last two returned scores are within tolerance; the return limit may split a near-tie.")
    if eligible and omitted_top:
        warnings.append(f"The current selection policy omits {omitted_top} returned genomes within tolerance of the best score.")
    if eligible and any(cluster["selection_splits_near_tie"] for cluster in clusters):
        warnings.append("The current selection policy splits a near-tie within at least one cluster.")
    return {
        "schema_version": "2.0.0", "profile": "distribution-v2",
        "status": "WARN" if warnings else ("PASS" if rows else "EMPTY"),
        "scope": "Returned database representatives only; not all isolates or a strain assignment.",
        "sources": [{"path": str(path), "sha256": sha256_file(path)} for path in (source, candidate_source)],
        "tolerance": {
            "hash_number": hash_number, "tie_tolerance_hashes": tolerance_hashes,
            "score_tolerance": tolerance, "formula": "tie_tolerance_hashes / hash_number",
            "interpretation": "Experimental descriptive Mashpit sketch-resolution heuristic; not a confidence interval or biological cutoff or cluster eligibility gate.",
        },
        "returned_genomes": len(rows), "duplicate_rows": duplicate_rows,
        "best_score": best, "lowest_score": rows[-1]["score"] if rows else None,
        "within_top_tolerance_count": sum(row["within_top_tolerance"] for row in rows),
        "retrieval": {
            "requested_limit": limit, "limit_reached": limit_reached,
            "additional_matches_exist": None,
            "tail_adjacent_gap": tail_gap, "tail_is_near_tie": tail_near_tie,
            "top_band_reaches_tail": top_band_reaches_tail,
            "may_split_near_tie": limit_reached and (bool(top_band_reaches_tail) or bool(tail_near_tie)),
        },
        "selection_preview": {
            "scope": "Current policy if SNP resolution is requested; not a download or comparison result.",
            "eligible": eligible, "policy": selection_policy,
            "mode": selected["mode"],
            "audit": selected,
            "selected_genomes": len(selected_accessions),
            "unexamined_returned_representatives": selected["unexamined_returned_representatives"],
            "unexamined_clusters": selected["unexamined_clusters"],
            "omitted_within_top_tolerance_count": omitted_top if eligible else None,
        },
        "rank_checkpoints": checkpoints, "clusters": clusters, "ranked": rows, "warnings": warnings,
    }


def render_plot(result: dict[str, Any], path: Path) -> dict[str, Any]:
    if not result["ranked"]:
        return {"status": "SKIPPED", "reason": "No representative scores to plot."}
    try:
        # Keep numerical diagnostics usable without plotting dependencies.
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        figure = Figure(figsize=(10, 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        rows = result["ranked"]
        axis.plot([row["rank"] for row in rows], [row["score"] for row in rows], color="#475569", linewidth=1)
        for selected, color, label in ((False, "#94a3b8", "Outside selection preview"), (True, "#0369a1", "In selection preview")):
            group = [row for row in rows if row["selected_by_current_policy"] == selected]
            if group:
                axis.scatter([row["rank"] for row in group], [row["score"] for row in group], s=16, color=color, label=label, zorder=3)
        best = result["best_score"]
        band_floor = max(0, best - result["tolerance"]["score_tolerance"])
        axis.axhspan(band_floor, best, color="#f59e0b", alpha=0.16, label="Within tolerance of best")
        for checkpoint in result["rank_checkpoints"]:
            if checkpoint["observed"]:
                axis.axvline(checkpoint["rank"], color="#cbd5e1", linestyle=":", linewidth=1)
        span = max(best - result["lowest_score"], best - band_floor, 0.001)
        axis.set_ylim(max(0, min(result["lowest_score"], band_floor) - span * .12), min(1.001, best + span * .15))
        axis.set_xlim(0.5, max(1.5, len(rows) + .5))
        axis.set_xlabel("Representative rank (highest similarity first)")
        axis.set_ylabel("Mashpit Jaccard similarity (zoomed scale; not ANI)")
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
        axis.set_title("Similarity distribution of returned representatives")
        axis.grid(axis="y", alpha=.18)
        axis.legend(loc="best", fontsize=9)
        figure.supxlabel(
            f"{len(rows)} genomes | {result['within_top_tolerance_count']} within top tolerance | "
            + ("Return limit reached; unseen matches unknown" if result["retrieval"]["limit_reached"] else "Return limit not reached"),
            fontsize=10,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=160)
        return {"status": "PASS", "path": str(path)}
    except Exception as error:
        return {"status": "FAIL", "error": str(error)}


def run_diagnostics(
    mashpit_output_dir: Path, database: dict[str, Any], workflow: dict[str, Any],
    selection_policy: dict[str, Any], output_dir: Path,
    selection_mode: str | None = None, query_accession: str | None = None,
) -> dict[str, Any]:
    try:
        result = analyze(
            mashpit_output_dir, database, workflow["mashpit"], selection_policy,
            workflow["similarity_diagnostics"]["rank_checkpoints"],
            selection_mode, query_accession,
        )
        result["plot"] = render_plot(result, output_dir / "rank_similarity.png")
        if result["plot"]["status"] == "FAIL":
            result["warnings"].append("Similarity plot could not be rendered: " + result["plot"]["error"])
            result["status"] = "WARN"
    except (WorkflowError, OSError, ValueError, KeyError, TypeError) as error:
        result = {"status": "ERROR", "error": str(error), "warnings": [f"Similarity diagnostics unavailable: {error}"]}
    write_json(output_dir / "summary.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mashpit-output-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--selection-mode", choices=("adaptive", "all_returned"))
    parser.add_argument("--query-accession")
    args = parser.parse_args()
    source = Path(args.mashpit_output_dir).expanduser().resolve()
    output = Path(args.output_dir).expanduser().resolve()
    if output.exists():
        parser.error("Output directory already exists; choose a new directory to preserve earlier results.")
    # Read the actual run's limit/tolerance/database rather than current defaults.
    database, profile = recorded_query_context(source)
    workflow = load_json(CONFIG_DIR / "workflow.json")
    workflow["mashpit"] = profile
    result = run_diagnostics(source, database, workflow, load_json(CONFIG_DIR / "snp-resolution-policy.json"), output,
                             args.selection_mode, args.query_accession)
    print(f"Similarity diagnostics: {result['status']}; {output / 'summary.json'}")
    return 2 if result["status"] == "ERROR" else 0


if __name__ == "__main__":
    raise SystemExit(main())
