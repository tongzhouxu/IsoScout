#!/usr/bin/env python3
"""Run the versioned fastp -> SKESA -> QUAST Illumina workflow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import CONFIG_DIR, WorkflowError, load_json, require_executable, run_logged, write_json
from validate_fastq import validate_pair


def run_workflow(r1: Path, r2: Path, output_dir: Path) -> dict:
    validation = validate_pair(r1, r2)
    if validation["status"] != "PASS":
        raise WorkflowError("Paired FASTQ validation failed: " + "; ".join(validation["errors"]))
    workflow = load_json(CONFIG_DIR / "workflow.json")
    settings = workflow["assembly"]
    fastp = require_executable("fastp")
    skesa = require_executable("skesa")
    quast = require_executable("quast.py")
    output_dir.mkdir(parents=True, exist_ok=True)
    trimmed_r1 = output_dir / "trimmed_R1.fastq.gz"
    trimmed_r2 = output_dir / "trimmed_R2.fastq.gz"
    fastp_json = output_dir / "fastp.json"
    commands: list[list[str]] = []
    fastp_command = [
        fastp,
        "--in1", str(r1), "--in2", str(r2),
        "--out1", str(trimmed_r1), "--out2", str(trimmed_r2),
        "--json", str(fastp_json), "--html", str(output_dir / "fastp.html"),
        *settings["fastp_args"],
    ]
    commands.append(fastp_command)
    if run_logged(fastp_command, output_dir, output_dir / "fastp.stdout.log", output_dir / "fastp.stderr.log"):
        raise WorkflowError("fastp failed; inspect assembly/fastp.stderr.log.")
    assembly = output_dir / "assembly.fasta"
    skesa_command = [
        skesa,
        "--fastq", f"{trimmed_r1},{trimmed_r2}",
        "--contigs_out", str(assembly),
        *settings["skesa_args"],
    ]
    commands.append(skesa_command)
    if run_logged(skesa_command, output_dir, output_dir / "skesa.stdout.log", output_dir / "skesa.stderr.log"):
        raise WorkflowError("SKESA failed; inspect assembly/skesa.stderr.log.")
    if not assembly.is_file() or assembly.stat().st_size == 0:
        raise WorkflowError("SKESA completed without producing a non-empty assembly.")
    quast_dir = output_dir / "quast"
    quast_command = [quast, "-o", str(quast_dir), *settings["quast_args"], str(assembly)]
    commands.append(quast_command)
    if run_logged(quast_command, output_dir, output_dir / "quast.stdout.log", output_dir / "quast.stderr.log"):
        raise WorkflowError("QUAST failed; inspect assembly/quast.stderr.log.")
    fastp_data = json.loads(fastp_json.read_text(encoding="utf-8"))
    summary = fastp_data.get("summary", {}).get("after_filtering", {})
    q30_rate = summary.get("q30_rate")
    result = {
        "status": "PASS",
        "assembly_path": str(assembly),
        "read_paths": [str(trimmed_r1), str(trimmed_r2)],
        "read_qc": {
            "status": "WARN" if q30_rate is not None and q30_rate < 0.80 else "PASS",
            "q30_fraction": q30_rate,
            "total_bases": summary.get("total_bases"),
            "read_pairs_validated": validation["read_pairs"],
        },
        "commands": commands,
        "workflow_version": workflow["workflow_version"],
    }
    write_json(output_dir / "assembly_workflow.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("r1")
    parser.add_argument("r2")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    try:
        result = run_workflow(Path(args.r1), Path(args.r2), Path(args.output_dir))
    except WorkflowError as error:
        result = {"status": "FAIL", "error": str(error)}
        write_json(Path(args.output_dir) / "assembly_workflow.json", result)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
