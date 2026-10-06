"""Versioned deterministic scheduling and evidence-aware stopping for SNP expansion.

Reference accession order is a reproducible sampling order, not a similarity
ranking. All eligibility comes from returned clusters and qualifying SNP data.
"""
from __future__ import annotations

from collections import defaultdict
import math
from typing import Any

from common import WorkflowError


def validate_expansion_policy(policy: dict[str, Any]) -> None:
    if policy.get("allocation_strategy") != "focused_balanced_v2":
        raise WorkflowError("Expected focused_balanced_v2 expansion allocation strategy.")
    for key in ("batch_size", "max_additional_reference_attempts", "max_rounds", "max_total_genomes",
                "stable_rounds_to_stop", "min_focused_qualifying_per_leader", "alternative_min_slots",
                "unresolved_priority_slots"):
        value = policy.get(key)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0 or (key not in {"unresolved_priority_slots"} and value == 0):
            raise WorkflowError(f"Expansion {key} must be a positive integer (zero allowed only for unresolved_priority_slots).")
    share = policy.get("leading_share")
    if not isinstance(share, (int, float)) or isinstance(share, bool) or not 0 < share < 1:
        raise WorkflowError("Expansion leading_share must be between zero and one.")
    if not isinstance(policy.get("complete_finite_leader_pool"), bool):
        raise WorkflowError("Expansion complete_finite_leader_pool must be Boolean.")
    if policy["alternative_min_slots"] >= policy["batch_size"]:
        raise WorkflowError("Expansion alternative_min_slots must be smaller than batch_size.")


def evidence_by_cluster(interpretation: dict[str, Any], selected_clusters: list[str]) -> dict[str, dict[str, int]]:
    evidence = {cluster: {"qualifying": 0, "excluded": 0} for cluster in selected_clusters}
    for row in interpretation.get("ranked", []):
        if row.get("cluster") in evidence:
            evidence[row["cluster"]]["qualifying"] += 1
    for row in interpretation.get("excluded_comparisons", []):
        if row.get("cluster") in evidence:
            evidence[row["cluster"]]["excluded"] += 1
    return evidence


def member_coverage(
    members: list[dict[str, Any]], initial_accessions: set[str], attempted: set[str],
    verified: set[str], interpretation: dict[str, Any], selected_clusters: list[str],
) -> list[dict[str, Any]]:
    """Count the finite discovered pool, including returned unselected references."""
    by_cluster: dict[str, set[str]] = {cluster: set() for cluster in selected_clusters}
    for row in members:
        if row["cluster"] in by_cluster:
            by_cluster[row["cluster"]].add(row["accession"])
    qualifying = {row["sample"] for row in interpretation.get("ranked", [])}
    excluded = {row["sample"] for row in interpretation.get("excluded_comparisons", [])}
    output = []
    for cluster in selected_clusters:
        pool = by_cluster[cluster] - initial_accessions
        output.append({
            "cluster": cluster, "available_additional": len(pool),
            "attempted_additional": len(pool & attempted),
            "downloaded_additional": len(pool & verified),
            "qualifying_additional": len(pool & qualifying),
            "excluded_additional": len(pool & excluded),
            "unexamined_additional": len(pool - attempted),
        })
    return output


def _rotate(clusters: list[str], after: str | None) -> list[str]:
    if after not in clusters:
        return list(clusters)
    index = clusters.index(after) + 1
    return clusters[index:] + clusters[:index]


def _draw(
    queues: dict[str, list[dict[str, Any]]], clusters: list[str], count: int,
    after: str | None, lane: str,
) -> tuple[list[dict[str, Any]], str | None, list[dict[str, str]]]:
    chosen: list[dict[str, Any]] = []
    audit: list[dict[str, str]] = []
    while len(chosen) < count:
        picked = False
        for cluster in _rotate(clusters, after):
            if queues.get(cluster):
                row = queues[cluster].pop(0)
                chosen.append(row)
                audit.append({"accession": row["accession"], "cluster": cluster, "lane": lane})
                after = cluster
                picked = True
                break
        if not picked:
            break
    return chosen, after, audit


def schedule_focused_batch(
    members: list[dict[str, Any]], attempted: set[str], selected_clusters: list[str],
    interpretation: dict[str, Any], policy: dict[str, Any], remaining_budget: int,
    remaining_rounds: int, state: dict[str, Any] | None = None,
    query_accession: str | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    """Reserve alternative capacity and rotate both lanes across rounds."""
    validate_expansion_policy(policy)
    state = dict(state or {})
    if remaining_budget < 0 or remaining_rounds < 0:
        raise WorkflowError("Expansion budget and remaining rounds cannot be negative.")
    seen: dict[str, dict[str, Any]] = {}
    for row in members:
        accession = row["accession"]
        if accession == query_accession or accession in attempted:
            continue
        if accession in seen and seen[accession]["cluster"] != row["cluster"]:
            raise WorkflowError(f"Conflicting cluster membership for {accession}.")
        seen[accession] = row
    queues: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in seen.values():
        if row["cluster"] in selected_clusters:
            queues[row["cluster"]].append(row)
    for queue in queues.values():
        queue.sort(key=lambda row: row["accession"])
    capacity = min(policy["batch_size"], remaining_budget, sum(len(queue) for queue in queues.values()))
    leaders = [cluster for cluster in selected_clusters if cluster in interpretation.get("nearest_clusters", [])]
    alternatives = [cluster for cluster in selected_clusters if cluster not in leaders]
    evidence = evidence_by_cluster(interpretation, selected_clusters)
    unresolved = [cluster for cluster in alternatives if evidence[cluster]["qualifying"] == 0]
    leader_pending = sum(len(queues[cluster]) for cluster in leaders)
    alternative_pending = sum(len(queues[cluster]) for cluster in alternatives)
    has_both = leader_pending > 0 and alternative_pending > 0
    alternate_reserve = min(policy["alternative_min_slots"], max(0, capacity - 1)) if has_both else 0
    total_future_capacity = min(remaining_budget, remaining_rounds * policy["batch_size"])
    # Completion is feasible only after reserving at least the stated alternative
    # slots in each remaining round. The policy may then use more focus slots
    # than its normal share, but never consumes the alternative reserve.
    rounds_needed = math.ceil(leader_pending / max(1, policy["batch_size"] - policy["alternative_min_slots"])) if has_both else 0
    future_reserve = min(alternative_pending, policy["alternative_min_slots"] * rounds_needed) if has_both else 0
    finite_completion_feasible = bool(leaders and leader_pending and policy["complete_finite_leader_pool"]
                                      and rounds_needed <= remaining_rounds
                                      and leader_pending <= total_future_capacity - future_reserve)
    if not leaders or not leader_pending:
        focus_target = 0
    elif not alternative_pending:
        focus_target = capacity
    elif finite_completion_feasible:
        focus_target = capacity - alternate_reserve
    else:
        focus_target = min(capacity - alternate_reserve, max(1, round(capacity * policy["leading_share"])))
    chosen: list[dict[str, Any]] = []
    decisions: list[dict[str, str]] = []
    leaders_chosen, state["leader_after"], decisions_part = _draw(queues, leaders, focus_target, state.get("leader_after"), "leading_snp_cluster")
    chosen.extend(leaders_chosen)
    decisions.extend(decisions_part)
    remaining = capacity - len(chosen)
    unresolved_quota = min(policy["unresolved_priority_slots"], remaining)
    unresolved_chosen, state["unresolved_after"], decisions_part = _draw(queues, unresolved, unresolved_quota, state.get("unresolved_after"), "unresolved_alternative")
    chosen.extend(unresolved_chosen)
    decisions.extend(decisions_part)
    remaining = capacity - len(chosen)
    alternative_chosen, state["alternative_after"], decisions_part = _draw(queues, alternatives, remaining, state.get("alternative_after"), "alternative_cluster")
    chosen.extend(alternative_chosen)
    decisions.extend(decisions_part)
    remaining = capacity - len(chosen)
    if remaining:
        extra_focus, state["leader_after"], decisions_part = _draw(queues, leaders, remaining, state.get("leader_after"), "unused_alternative_capacity_to_leader")
        chosen.extend(extra_focus)
        decisions.extend(decisions_part)
    audit = {
        "strategy": policy["allocation_strategy"], "leading_clusters": leaders,
        "alternative_clusters": alternatives, "unresolved_alternatives": unresolved,
        "qualifying_evidence_by_cluster": evidence, "leading_pending_before": leader_pending,
        "alternative_pending_before": alternative_pending, "capacity": capacity,
        "normal_leading_share": policy["leading_share"], "alternative_min_slots": alternate_reserve,
        "finite_leader_completion_feasible": finite_completion_feasible,
        "focus_target": focus_target, "chosen": decisions,
        "cursor_after": {key: state.get(key) for key in ("leader_after", "alternative_after", "unresolved_after")},
        "rationale": "Qualifying SNP minima define all tied leaders; reserved alternative slots rotate across selected clusters. Accession order is deterministic sampling, not a similarity ranking.",
    }
    return chosen, state, audit


def assess_progress(
    previous: dict[str, Any] | None, current: dict[str, Any],
    coverage: list[dict[str, Any]], policy: dict[str, Any],
    completion_feasible: bool, stable_rounds: int,
) -> dict[str, Any]:
    """An isolated exclusion does not block stability; a cluster without any
    qualifying evidence and with pending references does."""
    validate_expansion_policy(policy)
    present = {row["sample"] for row in current.get("ranked", [])}
    prior = {row["sample"] for row in previous.get("ranked", [])} if previous else set()
    new_qualifying = present - prior
    leaders = set(current.get("nearest_clusters", []))
    focused = {row["cluster"]: row["qualifying_additional"] for row in coverage}
    pending = {row["cluster"]: row["unexamined_additional"] for row in coverage}
    unresolved_pending = [row["cluster"] for row in coverage
                          if row["cluster"] not in leaders and row["unexamined_additional"]
                          and not any(item.get("cluster") == row["cluster"] for item in current.get("ranked", []))]
    below_min = [cluster for cluster in leaders
                 if focused.get(cluster, 0) < policy["min_focused_qualifying_per_leader"] and pending.get(cluster, 0)]
    if not current.get("ranked"):
        reason = "no_qualifying_evidence_continue_broad"
    elif current.get("coverage_blockers"):
        reason = "unresolved_coverage_exclusion"
    elif previous is None:
        reason = "initial_focused_exploration"
    elif (current.get("shared_regions") or {}).get("mask_sha256") != (previous.get("shared_regions") or {}).get("mask_sha256"):
        reason = "shared_target_regions_changed"
    elif current.get("nearest_snp_distance") != previous.get("nearest_snp_distance") or set(current.get("nearest_samples", [])) != set(previous.get("nearest_samples", [])):
        reason = "nearest_result_changed"
    elif not new_qualifying:
        reason = "no_new_qualifying_comparisons"
    elif below_min:
        reason = "focused_coverage_below_minimum"
    elif completion_feasible and any(pending.get(cluster, 0) for cluster in leaders):
        reason = "finite_leader_pool_completion_pending"
    elif unresolved_pending:
        reason = "unresolved_alternative_with_pending_references"
    else:
        reason = "nearest_set_unchanged_after_required_exploration"
    stable = reason == "nearest_set_unchanged_after_required_exploration"
    next_stable = stable_rounds + 1 if stable else 0
    return {
        "action": "STOP" if next_stable >= policy["stable_rounds_to_stop"] else "EXPAND",
        "reason": "stable_sampled_neighborhood" if next_stable >= policy["stable_rounds_to_stop"] else reason,
        "stable": stable, "stable_rounds": next_stable,
        "new_qualifying_comparisons": sorted(new_qualifying),
        "focused_qualifying_by_leader": {cluster: focused.get(cluster, 0) for cluster in sorted(leaders)},
        "leading_clusters_below_minimum": sorted(below_min),
        "unresolved_alternatives_with_pending_references": sorted(unresolved_pending),
        "finite_leader_completion_feasible": completion_feasible,
        "isolated_exclusions_do_not_block_stability": True,
    }
