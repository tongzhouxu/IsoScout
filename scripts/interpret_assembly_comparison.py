#!/usr/bin/env python3
"""Interpret target-only assembly distances without inventing a complete matrix."""
from __future__ import annotations
from typing import Any
from interpret_snp_resolution import cluster_summary


def interpret(comparisons: list[dict], targets: list[dict], best_cluster: str | None, policy: dict) -> dict[str,Any]:
    mapping = {row["accession"]:row["cluster"] for row in targets}
    ranked, excluded, warnings = [], [], []
    criteria = policy["comparability"]
    for raw in comparisons:
        row = {**raw,"cluster":mapping.get(raw["sample"])}
        reasons = []
        if row["status"] != "COMPARED": reasons.append("comparison_failed")
        elif row["target_aligned_bases"] == 0 or row["candidate_aligned_bases"] == 0: reasons.append("no_alignment")
        else:
            if row["target_aligned_fraction"] < criteria["min_target_aligned_fraction"]: reasons.append("insufficient_target_coverage")
            if row["candidate_aligned_fraction"] < criteria["min_candidate_aligned_fraction"]: reasons.append("insufficient_candidate_coverage")
        if row["cluster"] is None: reasons.append("unknown_reference_cluster")
        if reasons: excluded.append({**row,"reasons":reasons})
        else: ranked.append(row)
    ranked.sort(key=lambda row:(row["snp_distance"],row["sample"]))
    result = {"backend":policy["backend"],"comparison_scope":"target_to_candidate","metric":policy["metric"],
              "comparability_policy":criteria,"criteria_source":policy["criteria_source"],
              "ranked":ranked,"excluded_comparisons":excluded,"query_comparisons":comparisons,
              "newick_tree":None,"tree_status":"NOT_REQUESTED","warnings":warnings}
    if excluded: warnings.append(f"{len(excluded)} candidate comparisons were excluded or failed; see per-candidate reasons.")
    if not ranked:
        return {**result,"status":"INSUFFICIENT_DATA","nearest_samples":[],"nearest_clusters":[],"nearest_snp_distance":None}
    minimum = ranked[0]["snp_distance"]
    nearest = [row for row in ranked if row["snp_distance"]==minimum]
    nearest_ids = [row["sample"] for row in nearest]
    rate_min = min(row["snps_per_target_aligned_mb"] for row in ranked)
    rate_ids = [row["sample"] for row in ranked if row["snps_per_target_aligned_mb"]==rate_min]
    # A disagreement between minima is evidence of comparison-basis sensitivity,
    # not a reason to blend the two metrics into an unvalidated composite score.
    conflict = not (set(nearest_ids) & set(rate_ids))
    clusters = sorted({row["cluster"] for row in nearest})
    unique_cluster = clusters[0] if len(clusters)==1 and not conflict else None
    if conflict: warnings.append("SNP-count and SNP-rate minima disagree because aligned sequence differs; the nearest genome is unresolved across these metrics.")
    if len(nearest)>1: warnings.append("Several references share the minimum observed SNP count; no unique nearest genome is resolved.")
    statement = (f"Minimum observed assembly-comparison distance: {minimum} ACGT substitutions among examined qualifying candidates. "
                 "Aligned positions vary between pairs; aligned length is not a callable-site count. "
                 "Indels and ambiguous bases are reported separately; density and contig-edge SNP filters are not applied. "
                 "Zero observed SNPs do not establish whole-genome identity or outbreak membership.")
    if conflict: statement += " SNP-count and SNP-rate rankings disagree; report an unresolved closest-genome conclusion."
    return {**result,"status":"AMBIGUOUS" if len(nearest)>1 or conflict else "COMPARED",
            "nearest_samples":nearest_ids,"nearest_sample":nearest_ids[0] if len(nearest_ids)==1 and not conflict else None,
            "nearest_clusters":clusters,"nearest_cluster":unique_cluster,"nearest_snp_distance":minimum,
            "nearest_by_aligned_snp_rate":rate_ids,"ranking_basis_conflict":conflict,
            "agrees_with_mash_top_candidate":unique_cluster==best_cluster if unique_cluster and best_cluster else None,
            "cluster_summary":cluster_summary(ranked),"confidence":{"statement":statement}}
