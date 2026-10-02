#!/usr/bin/env python3
"""Compare adaptive selection with exhaustive SKA2 in a supplied held-out pool.

This measures search recall within that finite pool, not biological accuracy.
Reference files and query data remain local. No database or genome is downloaded.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, WorkflowError, load_json, sha256_file, write_json
from interpret_snp_resolution import interpret
from run_ska import run_ska
from select_snp_targets import select_targets


def benchmark(manifest_path: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise WorkflowError("Benchmark output must be a new directory.")
    manifest = load_json(manifest_path)
    root = manifest_path.resolve().parent
    def local(value):
        return str((root / value).resolve())
    references = manifest["references"]
    accessions = [row["accession"] for row in references]
    if len(set(accessions)) != len(accessions) or "QUERY" in accessions:
        raise WorkflowError("Benchmark reference IDs must be unique and cannot be QUERY.")
    if not manifest.get("held_out_query"):
        raise WorkflowError("The manifest must explicitly assert held_out_query after removing self/duplicate records.")
    query = manifest["query"]
    query_paths = [local(query)] if isinstance(query, str) else [local(path) for path in query]
    query_input = query_paths[0] if isinstance(query, str) else query_paths
    query_hashes = {sha256_file(Path(path)) for path in query_paths}
    genomes: dict[str, Any] = {"QUERY": query_input}
    checksums = []
    for row in references:
        path = local(row["path"])
        checksum = sha256_file(Path(path))
        if checksum in query_hashes:
            raise WorkflowError("A byte-identical query file appears in the reference pool.")
        genomes[row["accession"]] = path
        checksums.append({"accession": row["accession"], "path": path, "sha256": checksum})
    workflow = load_json(CONFIG_DIR / "workflow.json")
    selection = select_targets(Path(local(manifest["mashpit_output_dir"])), load_json(CONFIG_DIR / "snp-resolution-policy.json"))
    selected = {row["accession"] for row in selection["targets"]}
    if not selected.issubset(accessions):
        raise WorkflowError("The exhaustive pool must contain every adaptively selected reference.")
    output.mkdir(parents=True)
    ska = run_ska(genomes, output / "exhaustive", workflow["snp_resolution"]["kmer_size"])
    full = interpret(ska["distances"], references, None)
    subset_rows = [row for row in ska["distances"] if row["sample1"] in selected | {"QUERY"} and row["sample2"] in selected | {"QUERY"}]
    subset = interpret(subset_rows, [row for row in references if row["accession"] in selected], None)
    nearest = set(full.get("nearest_samples", []))
    expected = set(manifest.get("expected_nearest", []))
    if expected and not expected.issubset(accessions):
        raise WorkflowError("Expected nearest accessions must belong to the supplied pool.")
    result = {
        "schema_version": "1.0.0", "scope": "Exhaustive SKA2 comparison within a supplied finite pool; not biological validation.",
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_file(manifest_path)},
        "label_source": manifest.get("label_source"),
        "query_checksums": [{"path": path, "sha256": sha256_file(Path(path))} for path in query_paths],
        "reference_checksums": checksums, "selection": selection,
        "pool_size": len(references), "selected_count": len(selected),
        "exhaustive": full, "adaptive": subset,
        "nearest_neighbor_recall": len(nearest & selected) / len(nearest) if nearest else None,
        "all_exhaustive_nearest_retained": nearest.issubset(selected) if nearest else None,
        "expected_nearest": sorted(expected),
        "expected_set_recovered": set(subset.get("nearest_samples", [])) == expected if expected else None,
        "commands": ska["commands"],
    }
    write_json(output / "benchmark.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    try:
        result = benchmark(Path(args.manifest), Path(args.output_dir))
    except (WorkflowError, OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f"Benchmark failed: {error}\n")
    print(f"Nearest-neighbor recall in supplied pool: {result['nearest_neighbor_recall']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
