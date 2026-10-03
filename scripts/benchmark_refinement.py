#!/usr/bin/env python3
"""Compare frozen, revised, and all-returned selection with SKA2 in a finite pool.

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
from frozen_select_snp_targets_v2 import select_targets as select_frozen


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
    mashpit_output = Path(local(manifest["mashpit_output_dir"]))
    current_policy = load_json(CONFIG_DIR / "snp-resolution-policy.json")
    benchmark_budget = manifest.get("max_reference_genomes", current_policy["max_total_genomes"])
    if not isinstance(benchmark_budget, int) or isinstance(benchmark_budget, bool) or benchmark_budget <= 0:
        raise WorkflowError("Benchmark max_reference_genomes must be a positive integer.")
    if len(references) > benchmark_budget:
        raise WorkflowError("Finite reference pool exceeds the recorded benchmark reference budget.")
    selection = select_targets(mashpit_output, current_policy, query_accession=manifest.get("query_accession"))
    frozen = select_frozen(mashpit_output, load_json(CONFIG_DIR / "frozen-snp-resolution-policy-v2.json"))
    all_returned = select_targets(mashpit_output, current_policy, mode="all_returned", query_accession=manifest.get("query_accession"))
    if all_returned["status"] != "SELECTED":
        raise WorkflowError("All-returned baseline lacks sufficient recorded selection budget.")
    selected = {row["accession"] for row in selection["targets"]}
    selection_sets = {"frozen_v2": {row["accession"] for row in frozen["targets"]},
                      "revised_adaptive": selected,
                      "all_returned": {row["accession"] for row in all_returned["targets"]}}
    if not all(group.issubset(accessions) for group in selection_sets.values()):
        raise WorkflowError("The finite pool must contain every selected returned reference from each policy.")
    expected = set(manifest.get("expected_nearest", []))
    if expected and not expected.issubset(accessions):
        raise WorkflowError("Expected nearest accessions must belong to the supplied pool.")
    output.mkdir(parents=True)
    ska = run_ska(genomes, output / "exhaustive", workflow["snp_resolution"]["kmer_size"])
    full = interpret(ska["distances"], references, None)
    subset_rows = [row for row in ska["distances"] if row["sample1"] in selected | {"QUERY"} and row["sample2"] in selected | {"QUERY"}]
    subset = interpret(subset_rows, [row for row in references if row["accession"] in selected], None)
    nearest = set(full.get("nearest_samples", []))
    comparisons = {}
    for name, chosen in selection_sets.items():
        rows = [row for row in ska["distances"] if row["sample1"] in chosen | {"QUERY"} and row["sample2"] in chosen | {"QUERY"}]
        compared = interpret(rows, [row for row in references if row["accession"] in chosen], None)
        comparisons[name] = {
            "selected_count": len(chosen), "selected_accessions": sorted(chosen),
            "nearest_neighbor_recall": len(nearest & chosen) / len(nearest) if nearest else None,
            "all_exhaustive_nearest_retained": nearest.issubset(chosen) if nearest else None,
            "expected_set_recovered": set(compared.get("nearest_samples", [])) == expected if expected else None,
            "interpretation": compared,
        }
    result = {
        "schema_version": "2.0.0", "scope": "Three selection policies versus SKA2 in a supplied finite pool; not biological validation.",
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_file(manifest_path)},
        "implementations": {
            "frozen_v2": {"selector_sha256": sha256_file(Path(__file__).with_name("frozen_select_snp_targets_v2.py")),
                          "policy_sha256": sha256_file(CONFIG_DIR / "frozen-snp-resolution-policy-v2.json")},
            "revised_and_all_returned": {"selector_sha256": sha256_file(Path(__file__).with_name("select_snp_targets.py")),
                                         "policy_sha256": sha256_file(CONFIG_DIR / "snp-resolution-policy.json")},
            "refinement_policy_sha256": sha256_file(CONFIG_DIR / "refinement-policy.json"),
        },
        "label_source": manifest.get("label_source"),
        "query_checksums": [{"path": path, "sha256": sha256_file(Path(path))} for path in query_paths],
        "reference_checksums": checksums, "selection": selection,
        "pool_size": len(references), "benchmark_reference_budget": benchmark_budget,
        "selected_count": len(selected),
        "exhaustive": full, "adaptive": subset,
        "policy_comparison": comparisons,
        "selection_audits": {"frozen_v2": frozen, "revised_adaptive": selection, "all_returned": all_returned},
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
