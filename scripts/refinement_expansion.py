"""Deterministic expansion scheduling and feedback; no hidden confidence claims."""
from __future__ import annotations
from typing import Any


def next_batch(members: list[dict[str, Any]], attempted: set[str], clusters: list[str], batch_size: int, remaining_budget: int) -> list[dict[str, Any]]:
    queues = {cluster: sorted((row for row in members if row["cluster"] == cluster and row["accession"] not in attempted), key=lambda row: row["accession"]) for cluster in clusters}
    batch = []
    for offset in range(max((len(rows) for rows in queues.values()), default=0)):
        for cluster in clusters:
            if offset < len(queues[cluster]) and len(batch) < min(batch_size, remaining_budget):
                batch.append(queues[cluster][offset])
    return batch


def feedback(previous: dict[str, Any] | None, current: dict[str, Any]) -> dict[str, Any]:
    if not current.get("ranked"):
        return {"action": "STOP", "reason": "insufficient_comparability", "stable": False}
    if previous is None:
        return {"action": "EXPAND", "reason": "initial_member_exploration", "stable": False}
    reasons = []
    if previous.get("nearest_snp_distance") is None or current["nearest_snp_distance"] < previous["nearest_snp_distance"]:
        reasons.append("closer_neighbor_found")
    elif current["nearest_snp_distance"] != previous["nearest_snp_distance"]:
        reasons.append("nearest_distance_changed")
    if set(current.get("nearest_samples", [])) != set(previous.get("nearest_samples", [])):
        reasons.append("nearest_set_changed")
    if len(current.get("nearest_clusters", [])) > 1:
        reasons.append("clusters_remain_tied")
    if current.get("excluded_comparisons"):
        reasons.append("incomplete_comparability")
    return {"action": "EXPAND", "reason": ";".join(reasons) if reasons else "nearest_set_unchanged", "stable": not reasons}
