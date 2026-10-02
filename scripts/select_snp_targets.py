#!/usr/bin/env python3
"""Select an adaptive, bounded reference set and retain a decision audit."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from common import WorkflowError, load_json, sha256_file, write_json
from mashpit_similarity import read_scores, recorded_query_context, score_tolerance, within_tolerance
from parse_mashpit_results import load_candidates, locate_candidate_file, locate_representative_file


def _positive_integer(policy: dict[str, Any], key: str) -> int:
    value = policy.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise WorkflowError(f"Selection policy {key} must be a positive integer.")
    return value


def _bounded_selection(rows: list[dict[str, Any]], clusters: list[str], cap: int, tolerance: float) -> set[str]:
    """Reserve cluster leaders, prioritize higher score bands, share a cut band."""
    selected: set[str] = set()
    for cluster in clusters:
        leader = next((row for row in rows if row["cluster"] == cluster), None)
        if leader is not None and len(selected) < cap:
            selected.add(leader["accession"])
    remaining = [row for row in rows if row["accession"] not in selected]
    index = 0
    while index < len(remaining) and len(selected) < cap:
        end = index + 1
        while end < len(remaining) and within_tolerance(remaining[end - 1]["score"] - remaining[end]["score"], tolerance):
            end += 1
        band = remaining[index:end]
        # Fair deterministic allocation only within a score band; this is a
        # resource-budget decision, not evidence of a preferred biological match.
        queues = {cluster: [row for row in band if row["cluster"] == cluster] for cluster in clusters}
        for offset in range(max((len(queue) for queue in queues.values()), default=0)):
            for cluster in clusters:
                if offset < len(queues[cluster]) and len(selected) < cap:
                    selected.add(queues[cluster][offset]["accession"])
        index = end
    return selected


def select_targets(
    mashpit_output_dir: Path, policy: dict[str, Any],
    database: dict[str, Any] | None = None, profile: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if policy.get("strategy") != "adaptive-boundary-v1" or policy.get("policy_version") != "2.0.0":
        raise WorkflowError("Expected version 2.0.0 adaptive-boundary-v1 selection policy.")
    baseline = _positive_integer(policy, "initial_target_genomes")
    cap = _positive_integer(policy, "max_total_genomes")
    if baseline > cap:
        raise WorkflowError("initial_target_genomes cannot exceed max_total_genomes.")
    if database is None or profile is None:
        if database is not None or profile is not None:
            raise WorkflowError("Supply both database and profile, or use the recorded run context.")
        database, profile = recorded_query_context(mashpit_output_dir)
    tolerance = score_tolerance(database, profile)
    limit = _positive_integer(profile, "number")
    threshold = profile.get("threshold")
    if not isinstance(threshold, (int, float)) or not 0 <= threshold <= 1:
        raise WorkflowError("Recorded screening threshold must lie between zero and one.")
    source = locate_representative_file(mashpit_output_dir)
    candidate_source = locate_candidate_file(mashpit_output_dir)
    rows, duplicates = read_scores(source)
    candidates = load_candidates(candidate_source)
    if any(not 0 <= candidate["score"] <= 1 for candidate in candidates):
        raise WorkflowError("Candidate scores must be finite and lie between zero and one.")
    candidates.sort(key=lambda row: (-row["score"], row["cluster"]))
    best = rows[0]["score"] if rows else None
    candidate_best = candidates[0]["score"] if candidates else None
    eligible = best is not None and candidate_best is not None and candidate_best >= threshold
    # Retain upstream near_top alternatives, plus the inclusive score band to
    # avoid losing a cluster through floating-point rounding at its boundary.
    relevant = {candidate["cluster"] for candidate in candidates if candidate["near_top"]} if eligible else set()
    if eligible:
        relevant.add(candidates[0]["cluster"])
        relevant.update(row["cluster"] for row in rows if within_tolerance(best - row["score"], tolerance))
    leaders = {cluster: max((row["score"] for row in rows if row["cluster"] == cluster), default=-1) for cluster in relevant}
    clusters = sorted(relevant, key=lambda cluster: (-leaders[cluster], cluster))
    pool = [row for row in rows if row["cluster"] in relevant]
    leader_indices = [next((i for i, row in enumerate(pool) if row["cluster"] == cluster), -1) for cluster in clusters]
    initial_end = min(len(pool), baseline)
    end = max([initial_end] + [index + 1 for index in leader_indices])
    seed_end = end
    # Extend through adjacent near ties until a gap exceeds the tolerance.
    # This intentionally can span more than one tolerance from the first score;
    # the diagnostics' separate top-band count always anchors to the best score.
    while 0 < end < len(pool) and within_tolerance(pool[end - 1]["score"] - pool[end]["score"], tolerance):
        end += 1
    desired = pool[:end]
    selected = _bounded_selection(desired, clusters, cap, tolerance)
    desired_accessions = {row["accession"] for row in desired}
    baseline_accessions = {row["accession"] for row in pool[:initial_end]}
    leader_accessions = {pool[index]["accession"] for index in leader_indices if index >= 0}
    cap_omitted = len(desired) - len(selected)
    missing_clusters = [cluster for cluster in clusters if not any(row["cluster"] == cluster and row["accession"] in selected for row in pool)]
    warnings = []
    if cap_omitted:
        warnings.append(f"Adaptive selection requires {len(desired)} returned genomes but the hard ceiling is {cap}; {cap_omitted} were omitted. Candidate coverage is incomplete.")
    if missing_clusters:
        warnings.append("Plausible clusters lack selected representatives: " + ", ".join(missing_clusters))
    if len(rows) >= limit:
        warnings.append("Mashpit's return limit was reached; adaptive selection covers only returned representatives and additional candidates may be unobserved.")
    if duplicates:
        warnings.append(f"Collapsed {duplicates} duplicate accession rows before adaptive selection.")
    decisions = []
    for row in rows:
        accession = row["accession"]
        reasons = []
        if accession in selected:
            if accession in leader_accessions:
                reasons.append("plausible_cluster_representative")
            if accession in baseline_accessions:
                reasons.append("initial_ranked_target")
            if accession not in baseline_accessions and accession not in leader_accessions:
                reasons.append("boundary_near_tie_expansion" if pool.index(row) >= seed_end else "cluster_coverage_expansion")
        elif not eligible:
            reasons.append("no_candidate_above_screening_gate")
        elif row["cluster"] not in relevant:
            reasons.append("outside_plausible_clusters")
        elif accession in desired_accessions:
            reasons.append("hard_resource_ceiling")
        else:
            reasons.append("beyond_boundary_score_gap")
        decisions.append({**row, "selected": accession in selected, "reasons": reasons})
    return {
        "schema_version": "2.0.0", "status": "SELECTED" if selected else "SKIPPED",
        **({"reason": "No eligible representative genomes to resolve."} if not selected else {}),
        "policy": dict(policy), "scope": "Returned Mashpit representatives; selection is not strain assignment.",
        "query_context": {"mashpit_profile": profile, "database_settings": database.get("mashpit_database_settings"), "database_version": database.get("version")},
        "sources": [{"path": str(path), "sha256": sha256_file(path)} for path in (source, candidate_source)],
        "score_tolerance": tolerance, "relevant_clusters": clusters,
        "returned_genomes": len(rows), "duplicate_rows_removed": duplicates,
        "eligible_genomes": len(pool), "desired_genomes": len(desired),
        "selected_genomes": len(selected), "hard_ceiling_omitted": cap_omitted,
        "missing_plausible_clusters": missing_clusters,
        "desired_set_complete": cap_omitted == 0 and not missing_clusters,
        "retrieval_limit_reached": len(rows) >= limit,
        "boundary": {
            "initial_end_rank_in_eligible_pool": initial_end,
            "cluster_coverage_end_rank": seed_end,
            "expanded_end_rank": end,
            "score": pool[end - 1]["score"] if end else None,
            "gap_to_next": pool[end - 1]["score"] - pool[end]["score"] if 0 < end < len(pool) else None,
            "stop_reason": "score_gap" if end < len(pool) else "end_of_returned_eligible_pool",
        },
        "targets": [row for row in rows if row["accession"] in selected],
        "decisions": decisions,
        "duplicate_decisions": [{"source_row": line, "accession": row["accession"], "selected": False, "reason": "duplicate_accession", "canonical_source_row": row["source_rows"][0]} for row in rows for line in row["source_rows"][1:]],
        "warnings": warnings,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mashpit-output-dir", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        result = select_targets(Path(args.mashpit_output_dir), load_json(Path(args.policy)))
    except (WorkflowError, OSError, ValueError, KeyError, TypeError) as error:
        result = {"status": "ERROR", "error": str(error), "targets": []}
    write_json(Path(args.output), result)
    return 0 if result["status"] != "ERROR" else 2


if __name__ == "__main__":
    raise SystemExit(main())
