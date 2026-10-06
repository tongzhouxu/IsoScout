#!/usr/bin/env python3
"""Interpret target-relative candidates using bounded shared-region challenges."""
from __future__ import annotations
from common import WorkflowError
from shared_target_regions import challenge_candidates


def interpret(comparisons, targets, best_cluster, policy):
    ranking = policy.get('ranking') or {}
    if ranking.get('basis') != 'common_finalist_regions':
        raise WorkflowError('Assembly policy must select common_finalist_regions')
    mapping = {row['accession']: row['cluster'] for row in targets}
    usable, excluded, warnings = [], [], []
    criteria = policy['comparability']
    for raw in comparisons:
        row = {**raw, 'cluster': mapping.get(raw['sample'])}
        if row['status'] != 'COMPARED':
            excluded.append({**row,'reasons':['comparison_failed']})
        elif not row.get('target_aligned_bases') or not row.get('candidate_aligned_bases'):
            excluded.append({**row,'reasons':['no_alignment']})
        else:
            flags = []
            if row['target_aligned_fraction'] < criteria['min_target_aligned_fraction']:
                flags.append('target_coverage_below_reporting_reference')
            if row['candidate_aligned_fraction'] < criteria['min_candidate_aligned_fraction']:
                flags.append('candidate_coverage_below_reporting_reference')
            usable.append({**row,'coverage_flags':flags})
    if any(row['coverage_flags'] for row in usable):
        warnings.append('Some references have lower alignment coverage. They remain in candidate challenges; the coverage reporting reference does not exclude them.')
    if excluded:
        warnings.append(f'{len(excluded)} failed or unaligned references remain unresolved alternatives.')
    result = {'backend':policy['backend'],'comparison_scope':'target_to_candidate',
              'metric':'SNP counts on common finalist regions with lower bounds for outside candidates',
              'ranking_basis':'common_finalist_regions',
              'comparability_policy':criteria,'criteria_source':policy['criteria_source'],
              'ranked':usable,'excluded_comparisons':excluded,'query_comparisons':comparisons,
              'newick_tree':None,'tree_status':'NOT_REQUESTED','warnings':warnings,
              'status':'INSUFFICIENT_DATA','nearest_samples':[],'nearest_sample':None,
              'nearest_clusters':[],'nearest_cluster':None,'nearest_snp_distance':None,
              'nearest_by_aligned_snp_rate':[],'nearest_by_aligned_snp_rate_clusters':[],
              'cluster_status':'INSUFFICIENT_DATA','genome_status':'INSUFFICIENT_DATA',
              'cluster_candidates':[],'ranking_basis_conflict':False,'coverage_blockers':[],
              'decision_reference_samples':[],'agrees_with_mash_top_candidate':None,
              'conclusion_scope':'examined_reference_pool_only',
              'ordering':'display_groups_not_total_distance_order'}
    if not usable:
        result['confidence']={'statement':'No usable target comparisons; no closest reference or cluster is assigned.'}
        return result
    audit = challenge_candidates(usable, ranking)
    result['candidate_challenges'] = audit
    by = {r['sample']:r for r in usable + excluded}
    supported = audit['supported_samples']
    uncertain = sorted(set(audit['unresolved_samples']) | {r['sample'] for r in excluded})
    decision = sorted(set(supported) | set(uncertain))
    # A cluster witness uses one recorded mask: every reference that could
    # meet its minimum has the same stored label. No expected label is provided.
    witnesses = []
    steps = audit['verification_passes'][-1:] if audit['status'] == 'VERIFIED' else audit['verification_passes']
    for step in steps:
        alternatives = sorted(set(step['potential_samples']) | set(audit.get('region_sensitivity',{}).get('alternatives',[])) | {r['sample'] for r in excluded})
        label = by[step['anchor']]['cluster']
        if label and all(by[name]['cluster'] == label for name in alternatives):
            witnesses.append({'anchor':step['anchor'],'cluster':label,'mask_sha256':step['mask_sha256'],
                              'shared_target_bases':step['shared_target_bases'],
                              'same_cluster_unresolved_or_competing_samples':alternatives,
                              'external_references_strictly_beaten':sum(r['cluster'] != label for r in usable)})
    witness_labels = {w['cluster'] for w in witnesses}
    if len(witness_labels) > 1:
        witnesses = []  # Conflicting masks cannot certify one cluster.
    if witnesses and audit['status'] != 'VERIFIED':
        decision = sorted({name for w in witnesses for name in [w['anchor'], *w['same_cluster_unresolved_or_competing_samples']]})
    labels = sorted({by[name]['cluster'] for name in decision if by[name]['cluster']})
    unknown_label = any(not by[name]['cluster'] for name in decision)
    certified = audit['status'] == 'VERIFIED'
    cluster_resolved = (certified or bool(witnesses)) and len(labels) == 1 and not unknown_label
    genome_resolved = certified and len(supported) == 1 and not uncertain
    minimum_labels = sorted({by[name]['cluster'] for name in supported if by[name]['cluster']})
    if witnesses and not certified:
        minimum_labels = sorted(witness_labels)
    if audit.get('region_sensitivity',{}).get('alternatives'):
        warnings.append('The nearest-reference preference changes when additional aligned regions are included. Those alternatives remain in the final decision; one favorable region mask cannot establish a robust closest match.')
    if uncertain:
        warnings.append('Some candidate comparisons are unresolved. These alternatives are retained in the cluster conclusion rather than discarded.')
    if not certified:
        warnings.append('Common-region verification did not reach a supported fixed point. Genome resolution remains ambiguous; cluster support is reported separately.')
    rank_order = {name:0 if name in supported else 1 if name in uncertain else 2 for name in by}
    displayed = sorted(usable,key=lambda r:(rank_order[r['sample']],r['sample']))
    for row in displayed:
        row['challenge_role']='supported' if row['sample'] in supported else 'unresolved' if row['sample'] in uncertain else 'challenged_alternative'
    result.update(status='COMPARED' if genome_resolved and cluster_resolved else 'AMBIGUOUS',
                  ranked=displayed,nearest_samples=supported,
                  nearest_sample=supported[0] if genome_resolved else None,
                  nearest_clusters=minimum_labels,
                  nearest_cluster=labels[0] if cluster_resolved else None,
                  cluster_status='RESOLVED' if cluster_resolved else 'AMBIGUOUS',
                  genome_status='RESOLVED' if genome_resolved else 'AMBIGUOUS',
                  cluster_candidates=labels,decision_reference_samples=decision,
                  cluster_certificate={'status':'VERIFIED' if witnesses else 'NOT_ESTABLISHED',
                                       'basis':'all_potential_minima_on_one_recorded_mask_have_one_cluster_label',
                                       'witnesses':witnesses},
                  coverage_blockers=[{'sample':name,'cluster':by[name]['cluster'],'reason':'unresolved_candidate_evidence'} for name in uncertain],
                  agrees_with_mash_top_candidate=labels[0] == best_cluster if cluster_resolved and best_cluster else None)
    panel = audit.get('finalist_panel') or {}
    result['finalist_panel'] = panel
    for row in result['ranked']:
        row['finalist_snp_distance'] = panel.get('snp_counts',{}).get(row['sample'])
    result['confidence']={'statement':
        ('The finalists are compared on one identical set of target positions. Every outside candidate has been rechecked against that mask. '
         if certified else 'The bounded common-region comparison did not establish a unique nearest-genome conclusion. ')+
        ('A recorded common-region witness supports the stored cluster while genome alternatives remain unresolved. '
         if cluster_resolved and not genome_resolved else '')+
        'Observed SNP counts for partly overlapping alternatives are lower bounds on the panel region, not complete distances. '
        'A lower bound greater than the panel minimum rules out that alternative on these positions. '
        'A check on the larger pairwise shared regions retains contradictory candidates as unresolved alternatives. '
        'Other insufficient overlaps remain unresolved. Conclusions apply only to examined references and do not establish NCBI membership, strain identity or a globally closest genome.'}
    return result
