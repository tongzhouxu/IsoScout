#!/usr/bin/env python3
"""Select returned Mashpit representatives without a score-gap eligibility gate."""
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


def select_targets(
    mashpit_output_dir: Path, policy: dict[str, Any],
    database: dict[str, Any] | None = None, profile: dict[str, Any] | None = None,
    mode: str | None = None, query_accession: str | None = None,
) -> dict[str, Any]:
    if policy.get("strategy") != "global-ranked-cluster-coverage-v1" or policy.get("policy_version") != "3.0.0":
        raise WorkflowError("Expected version 3.0.0 global-ranked-cluster-coverage-v1 selection policy.")
    initial = _positive_integer(policy, "initial_target_genomes")
    returned_cap = _positive_integer(policy, "max_returned_representatives")
    total_cap = _positive_integer(policy, "max_total_genomes")
    if not initial <= returned_cap <= total_cap:
        raise WorkflowError("Require initial_target_genomes <= max_returned_representatives <= max_total_genomes.")
    coverage_strategy = policy.get("coverage_strategy")
    if coverage_strategy not in {"unrepresented_cluster_leaders_then_global_rank", "global_rank_only"}:
        raise WorkflowError("Unknown cluster coverage strategy.")
    mode = mode or policy.get("default_mode")
    if mode not in {"adaptive", "all_returned"}:
        raise WorkflowError("Selection mode must be adaptive or all_returned.")
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
    candidates.sort(key=lambda row: (-row["score"], row["cluster"]))
    if any(not 0 <= candidate["score"] <= 1 for candidate in candidates):
        raise WorkflowError("Candidate scores must be finite and lie between zero and one.")
    eligible = bool(rows and candidates and candidates[0]["score"] >= threshold)
    available = [row for row in rows if row["accession"] != query_accession]
    best = rows[0]["score"] if rows else None
    near_top = sorted({row["cluster"] for row in candidates if row["near_top"]})
    selected: set[str] = set()
    reasons: dict[str, str] = {}
    insufficient_budget = mode == "all_returned" and len(available) > total_cap
    if eligible and not insufficient_budget:
        if mode == "all_returned":
            for row in available:
                selected.add(row["accession"])
                reasons[row["accession"]] = "all_returned_benchmark"
        else:
            for row in available[:initial]:
                selected.add(row["accession"])
                reasons[row["accession"]] = "global_initial_rank"
            # Give each as-yet-unrepresented returned cluster one chance before
            # filling the remaining representative budget by global rank.
            if coverage_strategy == "unrepresented_cluster_leaders_then_global_rank":
                represented = {row["cluster"] for row in available if row["accession"] in selected}
                for row in available:
                    if len(selected) >= returned_cap:
                        break
                    if row["cluster"] not in represented:
                        selected.add(row["accession"])
                        reasons[row["accession"]] = "cluster_coverage_leader"
                        represented.add(row["cluster"])
            for row in available:
                if len(selected) >= returned_cap:
                    break
                if row["accession"] not in selected:
                    selected.add(row["accession"])
                    reasons[row["accession"]] = "additional_global_rank"
    selected_rows = [row for row in available if row["accession"] in selected]
    selected_clusters = sorted({row["cluster"] for row in selected_rows}, key=lambda cluster: (next(i for i, row in enumerate(available) if row["cluster"] == cluster), cluster))
    unexamined_clusters = sorted({row["cluster"] for row in available} - set(selected_clusters))
    warnings = []
    if insufficient_budget:
        warnings.append(f"All-returned mode requires {len(available)} reference attempts, exceeding the recorded total budget of {total_cap}; no partial baseline was selected.")
    elif eligible and len(selected) < len(available):
        warnings.append(f"Resource budget omitted {len(available) - len(selected)} returned representatives, including {len(unexamined_clusters)} entire clusters. No SNP dissimilarity is implied.")
    if len(rows) >= limit:
        warnings.append("Mashpit return limit reached; additional representatives may be unobserved.")
    if duplicates:
        warnings.append(f"Collapsed {duplicates} duplicate accession rows before selection.")
    if query_accession and not any(row["accession"] == query_accession for row in rows):
        warnings.append("Supplied query accession was not among returned representatives.")
    decisions = []
    for rank, row in enumerate(rows, 1):
        accession = row["accession"]
        if accession == query_accession:
            reason = "query_self_excluded"
        elif not eligible:
            reason = "no_candidate_above_screening_gate"
        elif insufficient_budget:
            reason = "all_returned_budget_insufficient"
        else:
            reason = reasons.get(accession, "returned_reference_budget")
        decisions.append({**row, "global_rank": rank, "selected": accession in selected,
                          "reason": reason, "gap_from_best": best - row["score"] if best is not None else None,
                          "within_top_tolerance": within_tolerance(best - row["score"], tolerance) if best is not None else None,
                          "upstream_near_top_cluster": row["cluster"] in near_top})
    return {
        "schema_version": "3.0.0", "status": "SELECTED" if selected else "SKIPPED",
        **({"reason": "All-returned budget is insufficient." if insufficient_budget else "No eligible representatives to resolve."} if not selected else {}),
        "policy": dict(policy), "mode": mode,
        "scope": "Returned Mashpit representatives only; no omitted cluster is ruled out by sketch score.",
        "query_context": {"mashpit_profile": profile, "database_settings": database.get("mashpit_database_settings"), "database_version": database.get("version"), "query_accession": query_accession},
        "sources": [{"path": str(path), "sha256": sha256_file(path)} for path in (source, candidate_source)],
        "score_tolerance": tolerance, "upstream_near_top_clusters": near_top,
        "returned_genomes": len(rows), "valid_nonself_genomes": len(available), "duplicate_rows_removed": duplicates,
        "retrieval_limit_reached": len(rows) >= limit, "initial_target_genomes": min(initial, len(available)),
        "selected_genomes": len(selected), "unexamined_returned_representatives": len(available) - len(selected),
        "selected_clusters": selected_clusters, "unexamined_clusters": unexamined_clusters,
        "resource_budget": {"mashpit_return_limit": limit, "initial_reference_target": initial,
                            "max_returned_representatives": returned_cap, "max_total_reference_attempts": total_cap,
                            "all_returned_required": len(available) if mode == "all_returned" else None,
                            "all_returned_budget_sufficient": not insufficient_budget if mode == "all_returned" else None},
        "targets": selected_rows, "decisions": decisions,
        "duplicate_decisions": [{"source_row": line, "accession": row["accession"], "selected": False,
                                  "reason": "duplicate_accession", "canonical_source_row": row["source_rows"][0]}
                                 for row in rows for line in row["source_rows"][1:]],
        "warnings": warnings,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mashpit-output-dir", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--mode", choices=("adaptive", "all_returned"))
    parser.add_argument("--query-accession", help="Known assembly accession for excluding the query itself.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        result = select_targets(Path(args.mashpit_output_dir), load_json(Path(args.policy)),
                                mode=args.mode, query_accession=args.query_accession)
    except (WorkflowError, OSError, ValueError, KeyError, TypeError) as error:
        result = {"status": "ERROR", "error": str(error), "targets": []}
    write_json(Path(args.output), result)
    return 0 if result["status"] != "ERROR" else 2


if __name__ == "__main__":
    raise SystemExit(main())
