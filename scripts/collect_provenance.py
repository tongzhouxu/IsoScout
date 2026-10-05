#!/usr/bin/env python3
"""Collect checksums and pinned/runtime versions for a screening run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, executable_version, load_json, sha256_file, utc_now, write_json


def collect(
    inputs: list[Path], assembly: Path | None, database_metadata: dict[str, Any] | None,
    commands: list[list[str]], extra: dict[str, Any] | None = None,
    refinement_policy_path: Path | None = None,
) -> dict[str, Any]:
    workflow = load_json(CONFIG_DIR / "workflow.json")
    qc_policy = load_json(CONFIG_DIR / "qc-policy.json")
    selection_policy = load_json(CONFIG_DIR / "snp-resolution-policy.json")
    refinement_policy_path = refinement_policy_path or CONFIG_DIR / "refinement-policy.json"
    refinement_policy = load_json(refinement_policy_path)
    provenance = {
        "timestamp_utc": utc_now(),
        "skill_version": workflow["skill_version"],
        "workflow_version": workflow["workflow_version"],
        "qc_policy_version": qc_policy["policy_version"],
        "snp_selection_policy": selection_policy,
        "refinement_policy": refinement_policy,
        "snp_selection_policy_sha256": sha256_file(CONFIG_DIR / "snp-resolution-policy.json"),
        "snp_selector_sha256": sha256_file(Path(__file__).with_name("select_snp_targets.py")),
        "refinement_scheduler_sha256": sha256_file(Path(__file__).with_name("refinement_expansion.py")),
        "refinement_policy_sha256": sha256_file(refinement_policy_path),
        "refinement_policy_path": str(refinement_policy_path.resolve()),
        "container": {
            "image": workflow["container_image"],
            "digest": workflow["container_digest"],
        },
        "pinned_tool_versions": workflow["tools"],
        "observed_tool_versions": {
            name: executable_version("quast.py" if name == "quast" else name)
            for name in workflow["tools"]
        },
        "inputs": [
            {"path": str(path.resolve()), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for path in inputs if path.is_file()
        ],
        "assembly": (
            {"path": str(assembly.resolve()), "sha256": sha256_file(assembly)}
            if assembly and assembly.is_file() else None
        ),
        "database": database_metadata,
        "mashpit_profile": workflow["mashpit"],
        "commands": commands,
    }
    if extra:
        provenance["run_state"] = extra
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--assembly")
    parser.add_argument("--database-metadata")
    parser.add_argument("--commands-json")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    database = load_json(Path(args.database_metadata)) if args.database_metadata else None
    commands = json.loads(Path(args.commands_json).read_text()) if args.commands_json else []
    result = collect(
        [Path(value) for value in args.input],
        Path(args.assembly) if args.assembly else None,
        database,
        commands,
    )
    write_json(Path(args.output), result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
