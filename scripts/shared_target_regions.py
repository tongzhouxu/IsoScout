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


def save_evidence(path, metrics, target_hash, candidate_hash, lengths):
    evidence = {"schema": "target-sites-v1", "coordinates": "one_based_inclusive",
                "target_sha256": target_hash, "candidate_sha256": candidate_hash,
                "target_lengths": lengths,
                "regions": metrics["target_unique_alignment_intervals"],
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
    if data.get("schema") != "target-sites-v1" or data.get("coordinates") != "one_based_inclusive":
        raise WorkflowError("Unsupported per-position evidence format")
    if data.get("target_sha256") != row.get("target_sha256") or data.get("candidate_sha256") != row.get("candidate_sha256"):
        raise WorkflowError("Shared-region evidence belongs to different inputs")
    lengths = data["target_lengths"]
    if any(type(n) is not int or n <= 0 for n in lengths.values()) or sum(lengths.values()) != row["target_bases"]:
        raise WorkflowError("Shared-region target lengths differ")
    regions = {name: merge(spans) for name, spans in data["regions"].items()}
    for name, spans in regions.items():
        if name not in lengths or any(hi > lengths[name] for lo, hi in spans):
            raise WorkflowError("Shared-region coordinates outside target")
    if sum(hi-lo+1 for spans in regions.values() for lo, hi in spans) != row["target_aligned_bases"]:
        raise WorkflowError("Shared-region coverage differs from pair metrics")
    snps = [tuple(site) for site in data["snps"]]
    if len(set(snps)) != len(snps) or len(snps) != row["snp_distance"]:
        raise WorkflowError("Shared-region SNP count differs from pair metrics")
    if any(type(pos) is not int or not within(regions, name, pos) for name, pos in snps):
        raise WorkflowError("SNP outside recorded target regions")
    return data, regions, snps


def challenge_candidates(rows, policy):
    """Find and verify a candidate with no observed challenger that beats it.

    Discovery is one linear scan. At most three verification passes reuse saved
    target calls, never align references or build an all-pairs distance matrix.
    Each challenge has one identical position mask for its two target distances.
    A verification certificate, not scan order, determines the conclusion.
    """
    min_bases = policy['min_shared_bases']
    fraction = policy['min_shared_fraction_of_larger_aligned_span']
    passes = policy['max_verification_passes']
    if type(min_bases) is not int or min_bases <= 0 or not 0 < fraction <= 1 or type(passes) is not int or not 1 <= passes <= 3:
        raise WorkflowError('Invalid candidate-challenge evidence or resource policy')
    by = {r['sample']: r for r in rows}
    if len(by) != len(rows):
        raise WorkflowError('Duplicate comparison sample')
    data, evidence_hashes, identity = {}, {}, None
    for name in sorted(by):
        saved, regions, snps = load_evidence(by[name])
        current = (saved['target_sha256'], saved['target_lengths'])
        if identity is not None and current != identity:
            raise WorkflowError('Candidates do not share one target coordinate system')
        identity = current
        data[name] = (regions, snps)
        evidence_hashes[name] = by[name]['site_evidence']['sha256']
    audited = {}

    def compare(left, right):
        first, second = sorted((left, right))
        key = (first, second)
        if key not in audited:
            shared = intersect(data[first][0], data[second][0])
            bases = sum(hi-lo+1 for spans in shared.values() for lo,hi in spans)
            ratio = bases / max(by[first]['target_aligned_bases'], by[second]['target_aligned_bases'])
            counts = [sum(within(shared, contig, pos) for contig,pos in data[name][1]) for name in key]
            sufficient = bases >= min_bases and ratio >= fraction
            winner = (first if counts[0] < counts[1] else second if counts[1] < counts[0] else None) if sufficient else None
            relation = 'FIRST' if winner == first else 'SECOND' if winner == second else 'TIE' if sufficient else 'INSUFFICIENT_DATA'
            audited[key] = {'first': first, 'second': second, 'shared_target_bases': bases,
                            'shared_fraction_of_larger_aligned_span': ratio,
                            'first_snps': counts[0], 'second_snps': counts[1],
                            'relation': relation, 'winner': winner,
                            'mask_sha256': hashlib.sha256(json.dumps(shared,sort_keys=True,separators=(',',':')).encode()).hexdigest()}
        return audited[key]

    def verdict(left, right):
        measured = compare(left, right)
        if measured['relation'] in ('TIE', 'INSUFFICIENT_DATA'):
            return measured['relation']
        return 'WIN' if measured['winner'] == left else 'LOSE'

    names = sorted(by)
    if not names:
        raise WorkflowError('No candidate position evidence')
    anchor = names[0]
    for candidate in names[1:]:
        relation = verdict(anchor, candidate)
        if relation == 'LOSE' or (relation in ('TIE','INSUFFICIENT_DATA') and
                                  by[candidate]['target_aligned_bases'] > by[anchor]['target_aligned_bases']):
            anchor = candidate
    history, visited = [], set()
    certified = False
    stop = 'verification_budget_exhausted'
    ties, unknown, defeats = [], [], []
    for iteration in range(passes):
        visited.add(anchor)
        ties, unknown, defeats = [], [], []
        for candidate in names:
            if candidate == anchor:
                continue
            relation = verdict(anchor, candidate)
            if relation == 'LOSE': defeats.append(candidate)
            elif relation == 'TIE': ties.append(candidate)
            elif relation == 'INSUFFICIENT_DATA': unknown.append(candidate)
        history.append({'anchor':anchor,'tied':ties,'unresolved':unknown,'defeaters':defeats})
        if not defeats:
            certified = True
            stop = 'verified_no_observed_defeater'
            break
        next_anchor = defeats[0]
        if next_anchor in visited:
            stop = 'comparison_cycle'
            break
        if iteration + 1 < passes:
            anchor = next_anchor
    # A single candidate has no opponent from which to obtain a support check.
    if len(names) == 1 and by[anchor]['target_aligned_bases'] < min_bases:
        certified = False
        stop = 'insufficient_sequence_evidence'
    supported = sorted({anchor, *ties}) if certified else []
    uncertain = sorted(set(unknown) | set(defeats) | (visited if not certified else set()))
    return {'status':'VERIFIED' if certified else 'UNRESOLVED',
            'basis':'candidate_anchor_shared_regions', 'anchor':anchor,
            'supported_samples':supported, 'unresolved_samples':uncertain,
            'stop_reason':stop, 'verification_passes':history,
            'comparisons':list(audited.values()), 'comparison_count':len(audited),
            'comparison_upper_bound':(passes+1)*max(0,len(names)-1),
            'evidence_sha256':evidence_hashes, 'policy':policy,
            'target_sha256':identity[0],
            'interpretation':'Each challenger is evaluated on its shared target positions with the anchor. There is no universal SNP count, total distance ranking or candidate-to-candidate alignment. Low overlap remains unresolved evidence.'}

