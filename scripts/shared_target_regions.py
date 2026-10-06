"""Compare saved target-relative SNP calls on identical positions per challenge.

No reference-to-reference alignment or pairwise distance matrix is constructed.
Intervals and substitution positions use one-based inclusive target coordinates.
"""
from __future__ import annotations
from bisect import bisect_right
import hashlib
import json
from pathlib import Path
from common import WorkflowError, load_json, sha256_file, write_json


def merge(spans):
    merged = []
    for lo, hi in sorted(spans):
        if not isinstance(lo, int) or not isinstance(hi, int) or not 1 <= lo <= hi:
            raise WorkflowError("Invalid shared-region interval")
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
        else:
            merged.append((lo, hi))
    return merged


def intersect(left, right):
    result = {}
    for contig in sorted(left.keys() & right.keys()):
        a, b = left[contig], right[contig]
        i = j = 0
        spans = []
        while i < len(a) and j < len(b):
            lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if lo <= hi:
                spans.append((lo, hi))
            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1
        if spans:
            result[contig] = spans
    return result


def within(regions, contig, position):
    spans = regions.get(contig, [])
    i = bisect_right(spans, (position, float("inf"))) - 1
    return i >= 0 and position <= spans[i][1]


def subtract(regions, excluded):
    """Remove closed intervals without expanding them into per-base sets."""
    result = {}
    for contig, spans in regions.items():
        cuts = merge(excluded.get(contig, []))
        kept = []
        j = 0
        for lo, hi in merge(spans):
            while j < len(cuts) and cuts[j][1] < lo:
                j += 1
            start, k = lo, j
            while k < len(cuts) and cuts[k][0] <= hi:
                left, right = cuts[k]
                if start < left:
                    kept.append((start, left - 1))
                start = max(start, right + 1)
                if start > hi:
                    break
                k += 1
            if start <= hi:
                kept.append((start, hi))
        if kept:
            result[contig] = kept
    return result


def save_evidence(path, metrics, target_hash, candidate_hash, lengths):
    evidence = {"schema": "target-sites-v2", "coordinates": "one_based_inclusive",
                "target_sha256": target_hash, "candidate_sha256": candidate_hash,
                "target_lengths": lengths,
                "region_basis": "unique_alignment_without_reported_target_gaps_or_ambiguous_differences",
                "regions": metrics["target_comparison_intervals"],
                "snps": metrics["snp_positions"]}
    write_json(path, evidence)
    return {"path": str(path), "sha256": sha256_file(path), "schema": evidence["schema"]}


def load_evidence(row):
    pointer = row.get("site_evidence") or {}
    if not pointer.get("path") or not pointer.get("sha256"):
        raise WorkflowError("Per-position evidence is missing; raw pair counts cannot resolve a shared-region ranking")
    path = Path(pointer["path"])
    if sha256_file(path) != pointer["sha256"]:
        raise WorkflowError("Shared-region evidence checksum mismatch")
    data = load_json(path)
    if data.get("schema") != "target-sites-v2" or data.get("coordinates") != "one_based_inclusive":
        raise WorkflowError("Unsupported per-position evidence format; regenerate from cached alignments to exclude target gaps")
    if data.get("target_sha256") != row.get("target_sha256") or data.get("candidate_sha256") != row.get("candidate_sha256"):
        raise WorkflowError("Shared-region evidence belongs to different inputs")
    lengths = data["target_lengths"]
    if any(type(n) is not int or n <= 0 for n in lengths.values()) or sum(lengths.values()) != row["target_bases"]:
        raise WorkflowError("Shared-region target lengths differ")
    regions = {name: merge(spans) for name, spans in data["regions"].items()}
    for name, spans in regions.items():
        if name not in lengths or any(hi > lengths[name] for lo, hi in spans):
            raise WorkflowError("Shared-region coordinates outside target")
    bases = sum(hi-lo+1 for spans in regions.values() for lo, hi in spans)
    if bases != row.get("target_comparable_bases") or bases > row["target_aligned_bases"]:
        raise WorkflowError("Shared-region coverage differs from pair metrics")
    snps = [tuple(site) for site in data["snps"]]
    if len(set(snps)) != len(snps) or len(snps) != row["snp_distance"]:
        raise WorkflowError("Shared-region SNP count differs from pair metrics")
    if any(type(pos) is not int or not within(regions, name, pos) for name, pos in snps):
        raise WorkflowError("SNP outside recorded target regions")
    return data, regions, snps


def challenge_candidates(rows, policy):
    """Verify finalists on one region mask with lower bounds for every outsider.

    A discovery scan only supplies a seed. The finalist panel grows monotonically;
    its common mask can shrink, so every outsider is rechecked after each change.
    A partial comparison can rule a reference out when its observed SNPs already
    exceed the panel minimum; missing positions cannot decrease that lower bound.
    """
    from array import array
    from bisect import bisect_left
    min_bases = policy['min_shared_bases']
    fraction = policy['min_shared_fraction_of_larger_aligned_span']
    passes = policy['max_verification_passes']
    if type(min_bases) is not int or min_bases <= 0 or not 0 < fraction <= 1 or type(passes) is not int or not 1 <= passes <= 3:
        raise WorkflowError('Invalid finalist-region evidence or resource policy')
    by = {r['sample']:r for r in rows}
    if not by or len(by) != len(rows):
        raise WorkflowError('Empty or duplicate comparison samples')
    data, evidence_hashes, identity = {}, {}, None
    for name in sorted(by):
        saved, regions, snps = load_evidence(by[name])
        current = (saved['target_sha256'], saved['target_lengths'])
        if identity is not None and current != identity:
            raise WorkflowError('Candidates do not share one target coordinate system')
        identity = current
        positions = {}
        for contig,pos in snps:
            positions.setdefault(contig,[]).append(pos)
        # Retain compact numeric indexes, rather than millions of position tuples.
        indexed = {c:array('Q',sorted(values)) for c,values in positions.items()}
        data[name] = (regions,indexed)
        evidence_hashes[name] = by[name]['site_evidence']['sha256']
    del saved, snps, positions

    def size(mask):
        return sum(hi-lo+1 for spans in mask.values() for lo,hi in spans)

    def digest(mask):
        return hashlib.sha256(json.dumps(mask,sort_keys=True,separators=(',',':')).encode()).hexdigest()

    def count(name,mask):
        positions = data[name][1]
        return sum(bisect_right(positions[c],hi)-bisect_left(positions[c],lo)
                   for c,spans in mask.items() if c in positions for lo,hi in spans)

    names=sorted(by)
    anchor=names[0]
    discovery=[]
    for candidate in names[1:]:
        shared=intersect(data[anchor][0],data[candidate][0])
        bases=size(shared)
        ratio=bases/max(by[anchor]['target_aligned_bases'],by[candidate]['target_aligned_bases'])
        left,right=count(anchor,shared),count(candidate,shared)
        sufficient=bases>=min_bases and ratio>=fraction
        winner=(anchor if left<right else candidate if right<left else None) if sufficient else None
        discovery.append({'first':anchor,'second':candidate,'first_snps':left,'second_snps':right,
                          'shared_target_bases':bases,'shared_fraction_of_larger_aligned_span':ratio,
                          'relation':'FIRST' if winner==anchor else 'SECOND' if winner==candidate else 'TIE' if sufficient else 'INSUFFICIENT_DATA',
                          'winner':winner,'mask_sha256':digest(shared),'phase':'seed_discovery'})
        if winner==candidate or (winner is None and by[candidate]['target_aligned_bases']>by[anchor]['target_aligned_bases']):
            anchor=candidate
    panel={anchor}
    history=[]
    total_checks=len(discovery)
    certified=False
    supported=[]
    unknown=[]
    proposed=[]
    stop='verification_budget_exhausted'
    final={}
    for iteration in range(passes):
        ordered=sorted(panel)
        mask=data[ordered[0]][0]
        for name in ordered[1:]:mask=intersect(mask,data[name][0])
        bases=size(mask)
        ratio=bases/max(by[n]['target_aligned_bases'] for n in panel)
        final={'samples':ordered,'regions':mask,'mask_sha256':digest(mask),
               'shared_target_bases':bases,'shared_fraction_of_largest_aligned_span':ratio,
               'snp_counts':{},'minimum_snps':None,'nearest_samples':[],
               'status':'SUFFICIENT' if bases>=min_bases and ratio>=fraction else 'INSUFFICIENT_DATA'}
        if final['status']!='SUFFICIENT':
            stop='insufficient_common_finalist_regions'
            break
        scores={name:count(name,mask) for name in ordered}
        minimum=min(scores.values())
        leaders=sorted(n for n in scores if scores[n]==minimum)
        anchor=leaders[0]
        final.update(snp_counts=scores,minimum_snps=minimum,nearest_samples=leaders)
        proposed=[];unknown=[];bounds=[]
        for candidate in names:
            if candidate in panel:continue
            shared=intersect(mask,data[candidate][0])
            shared_bases=size(shared)
            lower=count(candidate,shared)
            overlap=shared_bases/max(bases,by[candidate]['target_aligned_bases'])
            # This is a bound on the SAME panel mask, not an estimate obtained by
            # treating missing bases as matches or comparing different denominators.
            if shared_bases>=min_bases and lower>minimum:
                disposition='RULED_OUT_BY_SNP_LOWER_BOUND'
            elif shared_bases>=min_bases and overlap>=fraction:
                disposition='ADD_TO_FINALIST_PANEL';proposed.append(candidate)
            else:
                disposition='UNRESOLVED_OVERLAP';unknown.append(candidate)
            bounds.append({'sample':candidate,'shared_target_bases':shared_bases,
                           'shared_fraction_of_larger_span':overlap,'observed_snps_lower_bound':lower,
                           'panel_minimum_snps':minimum,'disposition':disposition})
        total_checks+=len(bounds)
        history.append({'panel_samples':ordered,'mask_sha256':final['mask_sha256'],
                        'shared_target_bases':bases,'panel_minimum_snps':minimum,'nearest_samples':leaders,
                        'anchor':anchor,'proposed_additions':proposed,'unresolved':unknown,
                        'potential_samples':sorted(set(leaders+proposed+unknown)),
                        'candidate_bounds':bounds})
        if not proposed:
            certified=True;supported=leaders;stop='verified_common_finalist_regions';break
        if iteration+1<passes:panel.update(proposed)
    sensitivity={'status':'NOT_EVALUATED','alternatives':[],'comparisons':[]}
    if certified and policy.get('check_region_sensitivity',True):
        sensitivity['status']='CHECKED'
        for candidate in names:
            if candidate==anchor:continue
            shared=intersect(data[anchor][0],data[candidate][0])
            shared_bases=size(shared)
            overlap=shared_bases/max(by[anchor]['target_aligned_bases'],by[candidate]['target_aligned_bases'])
            left,right=count(anchor,shared),count(candidate,shared)
            sufficient=shared_bases>=min_bases and overlap>=fraction
            contradictory=sufficient and right<=left and candidate not in supported
            if contradictory:sensitivity['alternatives'].append(candidate)
            sensitivity['comparisons'].append({'sample':candidate,'anchor_snps':left,'candidate_snps':right,
                                               'shared_target_bases':shared_bases,'shared_fraction':overlap,
                                               'mask_sha256':digest(shared),'sufficient':sufficient,
                                               'contradicts_unique_panel_preference':contradictory})
        total_checks+=len(sensitivity['comparisons'])
        if sensitivity['alternatives']:sensitivity['status']='REGION_SENSITIVE'
    uncertain=sorted(set(unknown) | set(sensitivity['alternatives']) |
                     (set(panel)|set(proposed) if not certified else set()))
    return {'status':'VERIFIED' if certified else 'UNRESOLVED','basis':'common_finalist_regions',
            'anchor':anchor,'supported_samples':supported,'unresolved_samples':uncertain,
            'stop_reason':stop,'verification_passes':history,'comparisons':discovery,
            'finalist_panel':final,'region_sensitivity':sensitivity,'comparison_count':total_checks,
            'comparison_upper_bound':(passes+2)*max(0,len(names)-1),
            'evidence_sha256':evidence_hashes,'policy':policy,'target_sha256':identity[0],
            'interpretation':'Finalists use one identical target mask. Every outsider is rechecked after mask changes. An observed SNP lower bound greater than the panel minimum rules an outsider out; other poorly overlapping alternatives remain unresolved. No candidate-to-candidate alignments or all-pairs matrix are constructed.'}
