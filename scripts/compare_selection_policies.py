#!/usr/bin/env python3
"""Reproduce old/revised/all-returned selection on saved Mashpit runs only.

Expected clusters are optional development-data audit labels; no selector sees them.
This script does not run SNP comparisons or claim accuracy.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, WorkflowError, load_json, sha256_file, write_json
from frozen_select_snp_targets_v2 import select_targets as select_frozen
from select_snp_targets import select_targets


def compare(manifest_path: Path, output_path: Path) -> dict[str, Any]:
    if output_path.exists():
        raise WorkflowError("Comparison output already exists; preserve prior results.")
    manifest = load_json(manifest_path)
    old_policy_path = CONFIG_DIR / "frozen-snp-resolution-policy-v2.json"
    new_policy_path = CONFIG_DIR / "snp-resolution-policy.json"
    old_policy, new_policy = load_json(old_policy_path), load_json(new_policy_path)
    cases = []
    for case in manifest["cases"]:
        source = (manifest_path.resolve().parent / case["mashpit_output_dir"]).resolve()
        accession = case.get("query_accession")
        old = select_frozen(source, old_policy)
        revised = select_targets(source, new_policy, mode="adaptive", query_accession=accession)
        baseline = select_targets(source, new_policy, mode="all_returned", query_accession=accession)
        expected = case.get("expected_cluster")
        def summarize(selection: dict[str, Any]) -> dict[str, Any]:
            targets = selection["targets"]
            return {
                "status": selection["status"], "selected_count": len(targets),
                "selected_accessions": [row["accession"] for row in targets],
                "selected_clusters": sorted({row["cluster"] for row in targets}),
                "expected_cluster_references_selected": sum(row["cluster"] == expected for row in targets) if expected else None,
                "expected_cluster_references_returned": sum(row["cluster"] == expected and row["accession"] != accession for row in selection["decisions"]) if expected else None,
                "selection_audit": selection,
            }
        cases.append({"case_id": case["case_id"], "mashpit_output_dir": str(source),
                      "expected_cluster_label_is_audit_only": expected,
                      "frozen_v2": summarize(old), "revised_adaptive": summarize(revised),
                      "all_returned": summarize(baseline)})
    result = {
        "schema_version": "1.0.0", "scope": "Selection coverage only; not SNP or strain-label accuracy.",
        "data_role": manifest.get("data_role", "development"),
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_file(manifest_path)},
        "implementations": {"frozen_v2": {"path": "scripts/frozen_select_snp_targets_v2.py", "source_sha256": sha256_file(Path(__file__).with_name("frozen_select_snp_targets_v2.py")), "policy_sha256": sha256_file(old_policy_path)},
                            "revised_adaptive": {"path": "scripts/select_snp_targets.py", "source_sha256": sha256_file(Path(__file__).with_name("select_snp_targets.py")), "policy_sha256": sha256_file(new_policy_path)},
                            "all_returned": {"path": "scripts/select_snp_targets.py", "source_sha256": sha256_file(Path(__file__).with_name("select_snp_targets.py")), "policy_sha256": sha256_file(new_policy_path)}},
        "cases": cases,
    }
    write_json(output_path, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        result = compare(Path(args.manifest), Path(args.output))
    except (WorkflowError, OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f"Comparison failed: {error}\n")
    print(f"Compared selection for {len(result['cases'])} saved runs; no SNP accuracy inferred.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
