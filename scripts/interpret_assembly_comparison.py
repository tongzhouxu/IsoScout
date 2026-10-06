#!/usr/bin/env python3
"""Interpret target-relative candidates using bounded shared-region challenges."""
from __future__ import annotations
from common import WorkflowError
from shared_target_regions import challenge_candidates


def interpret(comparisons, targets, best_cluster, policy):
    ranking = policy.get('ranking') or {}
    if ranking.get('basis') != 'candidate_anchor_shared_regions':
        raise WorkflowError('Assembly policy must select candidate_anchor_shared_regions')
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
              'metric':'paired target SNP counts on shared positions for each anchor challenge',
              'ranking_basis':'candidate_anchor_shared_regions',
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
    # Cluster and genome resolution are distinct. A verified pass may lose only
    # to references carrying the same stored cluster label while strictly
    # beating every reference outside that cluster. That witnesses a cluster
    # conclusion even when the within-cluster genome comparisons form a cycle.
    witnesses = []
    for step in audit['verification_passes']:
        label = by[step['anchor']]['cluster']
        alternatives = sorted(set(step['defeaters'] + step['tied'] + step['unresolved']) |
                              {r['sample'] for r in excluded})
        if label and by[step['anchor']]['target_aligned_bases'] >= ranking['min_shared_bases'] and all(by[name]['cluster'] == label for name in alternatives):
            witnesses.append({'anchor':step['anchor'],'cluster':label,
                              'same_cluster_unresolved_or_competing_samples':alternatives,
                              'external_references_strictly_beaten':sum(r['cluster'] != label for r in usable)})
    witness_labels = {w['cluster'] for w in witnesses}
    if len(witness_labels) > 1:
        raise WorkflowError('Conflicting cluster verification witnesses')
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
    if uncertain:
        warnings.append('Some candidate comparisons are unresolved. These alternatives are retained in the cluster conclusion rather than discarded.')
    if not certified:
        warnings.append('No individual reference passed bounded verification against every candidate. Genome resolution remains ambiguous; cluster support is reported separately.')
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
                                       'basis':'witness_strictly_beats_every_external_cluster_reference',
                                       'witnesses':witnesses},
                  coverage_blockers=[{'sample':name,'cluster':by[name]['cluster'],'reason':'unresolved_candidate_evidence'} for name in uncertain],
                  agrees_with_mash_top_candidate=labels[0] == best_cluster if cluster_resolved and best_cluster else None)
    result['confidence']={'statement':
        ('The supported reference has no observed challenger with fewer SNPs on the target regions shared by that pair. '
         if certified else 'No closest reference was verified against all examined candidates. ')+
        ('A cluster witness strictly beats every examined reference outside its stored cluster; within-cluster alternatives remain unresolved. '
         if cluster_resolved and not certified else '')+
        'Each comparison records both SNP counts and its shared sequence length. Low coverage alone does not exclude a candidate. '
        'Small or poorly overlapping comparisons remain unresolved. There is no single SNP distance comparable across all rows, '
        'and this result does not establish NCBI cluster membership, strain identity or a globally closest genome.'}
    return result
