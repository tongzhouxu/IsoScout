"""Compare saved target-relative SNP calls on one identical region mask.

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


def rank_on_shared_regions(rows, minimum_fraction):
    """Use every qualifying reference; never prune a reference to improve a mask."""
    mask = None
    sites = {}
    identity = None
    evidence_hashes = {}
    for row in sorted(rows, key=lambda r: r["sample"]):
        data, regions, snps = load_evidence(row)
        current = (data["target_sha256"], data["target_lengths"])
        if identity is not None and identity != current:
            raise WorkflowError("Candidates do not share one target coordinate system")
        identity = current
        mask = regions if mask is None else intersect(mask, regions)
        sites[row["sample"]] = snps
        evidence_hashes[row["sample"]] = row["site_evidence"]["sha256"]
    mask = {name: spans for name, spans in (mask or {}).items() if spans}
    bases = sum(hi-lo+1 for spans in mask.values() for lo, hi in spans)
    target_bases = sum(identity[1].values()) if identity else 0
    fraction = bases / target_bases if target_bases else 0
    sufficient = bool(bases and fraction >= minimum_fraction)
    digest = hashlib.sha256(json.dumps(mask, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    ranking = []
    for row in rows:
        count = sum(within(mask, name, pos) for name, pos in sites[row["sample"]])
        ranking.append({**row, "ranking_snp_distance": count,
                        "ranking_snps_per_target_mb": count*1e6/bases if bases else None})
    ranking.sort(key=lambda r: (r["ranking_snp_distance"], r["sample"]))
    audit = {"status": "SUFFICIENT" if sufficient else "INSUFFICIENT_DATA",
             "basis": "intersection_of_all_qualifying_target_regions",
             "coordinate_system": "one_based_inclusive", "target_sha256": identity[0] if identity else None,
             "target_bases": target_bases, "shared_target_bases": bases,
             "shared_target_fraction": fraction, "minimum_target_fraction": minimum_fraction,
             "mask_sha256": digest, "regions": mask, "evidence_sha256": evidence_hashes,
             "samples": sorted(sites),
             "interpretation": "Identical aligned target positions for every ranked reference; not an exact callable-site mask, strain cutoff or NCBI SNP distance."}
    return ranking, audit
