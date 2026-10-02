#!/usr/bin/env python3
"""Interpret ska2 pairwise SNP distances: rank, group by cluster, and build a tree.

Reports the full pairwise matrix ska2 already computes (not just query-vs-
representative), a per-cluster distance summary, a Neighbor-Joining tree over
every compared genome, and a structural confidence comparison (nearest
cluster's closest genome vs. the next-nearest cluster's closest genome).
No strength label ("high/medium/low confidence") is invented: only the raw
distances and their ratio are reported, consistent with
references/mashpit-interpretation.md's policy on the Mash score itself.
"""

from __future__ import annotations

import argparse
import json
import statistics
import math
from pathlib import Path
from typing import Any

from build_snp_tree import neighbor_joining
from common import CONFIG_DIR, WorkflowError, load_json, write_json


QUERY_SAMPLE_ID = "QUERY"


def query_distances(distance_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    results = []
    for row in distance_rows:
        if row["sample1"] == QUERY_SAMPLE_ID:
            other = row["sample2"]
        elif row["sample2"] == QUERY_SAMPLE_ID:
            other = row["sample1"]
        else:
            continue
        results.append({
            "sample": other,
            "snp_distance": row["snp_distance"],
            "mismatch_proportion": row["mismatch_proportion"],
            "match_count": row["match_count"],
            "mismatch_count": row["mismatch_count"],
        })
    results.sort(key=lambda item: item["snp_distance"])
    return results


def _display_label(sample: str, accession_to_cluster: dict[str, str]) -> str:
    cluster = accession_to_cluster.get(sample)
    return sample if cluster is None else f"{sample}_{cluster}"


def _relabel_rows(
    distance_rows: list[dict[str, Any]], accession_to_cluster: dict[str, str]
) -> list[dict[str, Any]]:
    relabeled = []
    for row in distance_rows:
        relabeled.append({
            **row,
            "sample1": _display_label(row["sample1"], accession_to_cluster),
            "sample2": _display_label(row["sample2"], accession_to_cluster),
        })
    return relabeled


def cluster_summary(ranked: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_cluster: dict[str, list[float]] = {}
    for item in ranked:
        by_cluster.setdefault(item["cluster"], []).append(item["snp_distance"])
    summary = [
        {
            "cluster": cluster,
            "genomes_compared": len(distances),
            "min_snp_distance": min(distances),
            "median_snp_distance": statistics.median(distances),
            "max_snp_distance": max(distances),
        }
        for cluster, distances in by_cluster.items()
    ]
    summary.sort(key=lambda item: item["min_snp_distance"])
    return summary


def build_confidence(summary: list[dict[str, Any]]) -> dict[str, Any]:
    nearest = summary[0]
    runner_up = summary[1] if len(summary) > 1 else None
    tied = [row["cluster"] for row in summary if row["min_snp_distance"] == nearest["min_snp_distance"]]
    ratio = (runner_up["min_snp_distance"] / nearest["min_snp_distance"] if runner_up and nearest["min_snp_distance"] > 0 else None)
    statement = (
        f"Minimum observed SKA2 distance: {nearest['min_snp_distance']:.2f} SNPs across comparable "
        f"split-kmer contexts; nearest cluster(s): {', '.join(tied)}. "
    )
    if len(tied) > 1:
        statement += "Multiple clusters are tied; no unique cluster assignment is supported. "
    elif runner_up:
        statement += f"Next cluster minimum: {runner_up['min_snp_distance']:.2f} SNPs. "
    else:
        statement += "Only one cluster had qualifying comparisons. "
    statement += "Zero observed differences do not establish whole-genome identity. Distances and ratios are not statistical confidence or outbreak confirmation."
    return {"nearest_cluster": tied[0] if len(tied) == 1 else None, "nearest_clusters": tied,
            "nearest_cluster_min_snp_distance": nearest["min_snp_distance"],
            "next_cluster": runner_up["cluster"] if runner_up else None,
            "next_cluster_min_snp_distance": runner_up["min_snp_distance"] if runner_up else None,
            "separation_ratio": ratio, "statement": statement}


def comparison_qualifies(row: dict[str, Any], policy: dict[str, Any]) -> bool:
    return (math.isfinite(row["snp_distance"]) and row["snp_distance"] >= 0
            and math.isfinite(row["mismatch_proportion"]) and 0 <= row["mismatch_proportion"] <= policy["max_mismatch_proportion"]
            and row["match_count"] >= policy["min_shared_split_kmers"] and row["mismatch_count"] >= 0)


def interpret(
    distance_rows: list[dict[str, Any]], targets: list[dict[str, Any]],
    mash_best_cluster: str | None, comparability: dict[str, Any] | None = None,
) -> dict[str, Any]:
    policy = comparability or load_json(CONFIG_DIR / "refinement-policy.json")["comparability"]
    all_ranked = query_distances(distance_rows)
    accession_to_cluster = {item["accession"]: item["cluster"] for item in targets}
    ranked, excluded = [], []
    for item in all_ranked:
        item["cluster"] = accession_to_cluster.get(item["sample"])
        if item["cluster"] is not None and comparison_qualifies(item, policy):
            ranked.append(item)
        else:
            excluded.append({**item, "reason": "insufficient_comparability_or_unknown_reference"})
    ranked.sort(key=lambda item: (item["snp_distance"], item["sample"]))
    warnings = []
    if excluded:
        warnings.append(f"{len(excluded)} query comparisons failed the configured overlap/validity checks and were excluded from nearest-neighbor ranking.")
    result = {"comparability_policy": policy, "excluded_comparisons": excluded,
              "pairwise_distances": distance_rows, "ranked": ranked, "warnings": warnings}
    if not ranked:
        return {**result, "status": "INSUFFICIENT_DATA", "nearest_samples": [], "nearest_clusters": [],
                "nearest_snp_distance": None, "newick_tree": None,
                "warnings": warnings + ["No query-to-reference comparisons meet the configured comparability checks."]}
    minimum = ranked[0]["snp_distance"]
    nearest = [item for item in ranked if item["snp_distance"] == minimum]
    clusters = sorted({item["cluster"] for item in nearest})
    unique_cluster = clusters[0] if len(clusters) == 1 else None
    agrees = (unique_cluster == mash_best_cluster) if unique_cluster and mash_best_cluster else None
    if len(nearest) > 1:
        warnings.append(f"{len(nearest)} references share the minimum observed SNP distance; no unique nearest isolate is resolved.")
    if agrees is False:
        warnings.append(f"The SNP-nearest cluster ({unique_cluster}) differs from Mashpit's top candidate ({mash_best_cluster}).")
    summary = cluster_summary(ranked)
    newick = None
    tree_status = "SKIPPED"
    allowed = {item["sample"] for item in ranked} | {QUERY_SAMPLE_ID}
    tree_rows = [row for row in distance_rows if row["sample1"] in allowed and row["sample2"] in allowed]
    if all(comparison_qualifies(row, policy) for row in tree_rows):
        try:
            newick = neighbor_joining(_relabel_rows(tree_rows, accession_to_cluster))
            tree_status = "PASS"
        except (WorkflowError, KeyError) as error:
            tree_status = str(error)
    if newick is None:
        warnings.append("SNP tree unavailable: a complete matrix of qualifying comparisons is required.")
    return {**result, "status": "AMBIGUOUS" if len(nearest) > 1 else "COMPARED",
        "nearest_sample": nearest[0]["sample"] if len(nearest) == 1 else None,
        "nearest_samples": [item["sample"] for item in nearest],
        "nearest_cluster": unique_cluster, "nearest_clusters": clusters,
        "nearest_snp_distance": minimum, "agrees_with_mash_top_candidate": agrees,
        "cluster_summary": summary, "confidence": build_confidence(summary),
        "newick_tree": newick, "tree_status": tree_status,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--distances-json", required=True)
    parser.add_argument("--targets-json", required=True)
    parser.add_argument("--mash-best-cluster")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    distances = json.loads(Path(args.distances_json).read_text(encoding="utf-8"))
    targets = json.loads(Path(args.targets_json).read_text(encoding="utf-8"))
    result = interpret(distances, targets, args.mash_best_cluster)
    write_json(Path(args.output), result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
