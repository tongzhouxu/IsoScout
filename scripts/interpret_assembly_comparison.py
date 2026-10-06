#!/usr/bin/env python3
"""Interpret cached target comparisons on one common set of target positions."""
from __future__ import annotations
from common import WorkflowError
from shared_target_regions import rank_on_shared_regions


def interpret(comparisons, targets, best_cluster, policy):
    if (policy.get('ranking') or {}).get('basis') != 'shared_target_regions':
        raise WorkflowError('Assembly policy must explicitly select shared_target_regions ranking')
    if policy['ranking'].get('coverage_exclusion_guard') is not True:
        raise WorkflowError('Coverage-excluded alternatives must remain visible to interpretation')
    mapping = {row['accession']: row['cluster'] for row in targets}
    eligible, excluded, warnings = [], [], []
    criteria = policy['comparability']
    for raw in comparisons:
        row = {**raw, 'cluster': mapping.get(raw['sample'])}
        reasons = []
        if row['status'] != 'COMPARED':
            reasons.append('comparison_failed')
        elif not row['target_aligned_bases'] or not row['candidate_aligned_bases']:
            reasons.append('no_alignment')
        else:
            if row['target_aligned_fraction'] < criteria['min_target_aligned_fraction']:
                reasons.append('insufficient_target_coverage')
            if row['candidate_aligned_fraction'] < criteria['min_candidate_aligned_fraction']:
                reasons.append('insufficient_candidate_coverage')
        if row['cluster'] is None:
            reasons.append('unknown_reference_cluster')
        if reasons:
            excluded.append({**row, 'reasons': reasons})
        else:
            eligible.append(row)
    if excluded:
        warnings.append(f'{len(excluded)} candidate comparisons were excluded or failed; see per-candidate reasons.')
    result = {'backend': policy['backend'], 'comparison_scope': 'target_to_candidate',
              'metric': 'ACGT substitutions on the same shared target regions',
              'ranking_basis': 'shared_target_regions', 'comparability_policy': criteria,
              'criteria_source': policy['criteria_source'], 'ranked': [],
              'excluded_comparisons': excluded, 'query_comparisons': comparisons,
              'newick_tree': None, 'tree_status': 'NOT_REQUESTED', 'warnings': warnings,
              'status': 'INSUFFICIENT_DATA', 'nearest_samples': [], 'nearest_sample': None,
              'nearest_clusters': [], 'nearest_cluster': None, 'nearest_snp_distance': None,
              'nearest_by_aligned_snp_rate': [], 'nearest_by_aligned_snp_rate_clusters': [],
              'cluster_status': 'INSUFFICIENT_DATA', 'genome_status': 'INSUFFICIENT_DATA',
              'cluster_candidates': [], 'ranking_basis_conflict': False,
              'coverage_blockers': [], 'agrees_with_mash_top_candidate': None,
              'conclusion_scope': 'examined_reference_pool_only'}
    if not eligible:
        result['shared_regions'] = {'status': 'INSUFFICIENT_DATA', 'reason': 'no_qualifying_references'}
        result['confidence'] = {'statement': 'No references satisfy the recorded comparison criteria; no nearest match is assigned.'}
        return result
    raw_min = min(r['snp_distance'] for r in eligible)
    rate_min = min(r['snps_per_target_aligned_mb'] for r in eligible)
    result['pair_specific_diagnostics'] = {
        'minimum_count_samples': sorted(r['sample'] for r in eligible if r['snp_distance'] == raw_min),
        'minimum_rate_samples': sorted(r['sample'] for r in eligible if r['snps_per_target_aligned_mb'] == rate_min),
        'interpretation': 'Counts and rates use different target positions across pairs; neither set determines the cluster call.'}
    ranked, shared = rank_on_shared_regions(eligible, criteria['min_target_aligned_fraction'])
    result.update(ranked=ranked, shared_regions=shared)
    if shared['status'] != 'SUFFICIENT':
        warnings.append('The common target regions do not meet the recorded coverage requirement. No candidate was removed to enlarge the common region set.')
        result['confidence'] = {'statement': 'Insufficient shared target sequence for a resolved comparison. Pair-specific counts remain diagnostic only.'}
        return result
    minimum = ranked[0]['ranking_snp_distance']
    nearest = [r for r in ranked if r['ranking_snp_distance'] == minimum]
    ids = sorted(r['sample'] for r in nearest)
    clusters = sorted({r['cluster'] for r in nearest})
    # A low-coverage candidate is not promoted. Its lower raw count or rate can
    # nevertheless invalidate a confident conclusion obtained by excluding it.
    blockers = []
    for row in excluded:
        if row['status'] != 'COMPARED' or not row.get('target_aligned_bases'):
            continue
        if not set(row['reasons']) & {'insufficient_target_coverage', 'insufficient_candidate_coverage', 'unknown_reference_cluster'}:
            continue
        if any(row['snp_distance'] <= winner['snp_distance'] or
               row['snps_per_target_aligned_mb'] <= winner['snps_per_target_aligned_mb'] for winner in nearest):
            blockers.append({'sample': row['sample'], 'cluster': row['cluster'],
                             'snp_distance': row['snp_distance'],
                             'snps_per_target_aligned_mb': row['snps_per_target_aligned_mb'],
                             'exclusion_reasons': row['reasons'],
                             'reason': 'excluded_alternative_has_no_larger_pair_specific_count_or_rate'})
    candidates = sorted(set(clusters) | {r['cluster'] for r in blockers if r['cluster']})
    missing_label = any(r['cluster'] is None for r in blockers)
    cluster_resolved = len(candidates) == 1 and not missing_label
    genome_resolved = len(nearest) == 1 and not blockers
    if blockers:
        warnings.append('A potentially competitive reference was excluded by comparison criteria. It is not a validated match, but the exclusion prevents a resolved genome claim and, when its label differs or is unknown, a resolved cluster claim.')
    statement = (f'Minimum {minimum} observed ACGT substitutions on {shared["shared_target_bases"]} identical target positions '
                 'among examined qualifying references. All exact ties are retained. These are aligned regions, '
                 'not a complete callable-site mask or NCBI SNP distances. The candidate pool may omit a closer match.')
    if blockers:
        statement += ' Coverage-excluded alternatives remain unresolved; do not report the qualifying minimum as a confirmed closest match.'
    # With a common denominator, count and rate minima are identical. Legacy
    # fields name these common-region minima; raw minima live only in diagnostics.
    result.update(status='COMPARED' if genome_resolved and cluster_resolved else 'AMBIGUOUS',
                  nearest_samples=ids, nearest_sample=ids[0] if genome_resolved else None,
                  nearest_clusters=clusters, nearest_cluster=candidates[0] if cluster_resolved else None,
                  nearest_snp_distance=minimum, nearest_by_aligned_snp_rate=ids,
                  nearest_by_aligned_snp_rate_clusters=clusters,
                  cluster_status='RESOLVED' if cluster_resolved else 'AMBIGUOUS',
                  genome_status='RESOLVED' if genome_resolved else 'AMBIGUOUS',
                  cluster_candidates=candidates, coverage_blockers=blockers,
                  agrees_with_mash_top_candidate=candidates[0] == best_cluster if cluster_resolved and best_cluster else None,
                  confidence={'statement': statement})
    return result
