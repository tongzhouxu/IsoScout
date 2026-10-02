#!/usr/bin/env python3
"""Build a merged ska2 split-kmer file and compute pairwise SNP distances.

Pinned to ska2 0.5.1's ``ska build``/``ska distance`` contract: a tab
separated ``name<TAB>path`` file list for ``build`` (see the ska2 module
docs), and the fixed long-form ``ska distance`` header
``Sample1  Sample2  Distance  Mismatches (proportion)  Match count  Mismatch count``
(see ``generic_modes.rs::distance`` in ska.rust). Re-verify both if the
pinned ska version changes.
"""

from __future__ import annotations

import argparse
import csv
import math
from itertools import combinations
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, WorkflowError, load_json, require_executable, run_logged, write_json


DISTANCE_HEADER = (
    "Sample1", "Sample2", "Distance", "Mismatches (proportion)", "Match count", "Mismatch count"
)


def build_file_list(genomes: dict[str, str | list[str]], path: Path) -> None:
    lines = []
    for sample_id, value in sorted(genomes.items()):
        paths = [value] if isinstance(value, str) else list(value)
        if (not isinstance(value, str) and len(paths) != 2) or not paths:
            raise WorkflowError("Read inputs must contain exactly two paired FASTQ paths.")
        fields = [sample_id, *paths]
        if any(not isinstance(field, str) or not field or any(c in field for c in "\t\r\n") for field in fields):
            raise WorkflowError("SKA file-list fields must be nonempty strings without tabs/newlines.")
        if isinstance(value, str) and value.lower().endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
            raise WorkflowError("SKA FASTQ input requires an R1/R2 pair, not a single path.")
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_distance_table(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or tuple(header) != DISTANCE_HEADER:
            raise WorkflowError(f"Unexpected ska distance header in {path}: {header}")
        rows = []
        for parts in reader:
            if len(parts) != len(DISTANCE_HEADER):
                raise WorkflowError(f"Malformed ska distance row in {path}: {parts}")
            sample1, sample2, distance, mismatch_prop, match_count, mismatch_count = parts
            rows.append({
                "sample1": sample1,
                "sample2": sample2,
                "snp_distance": float(distance),
                "mismatch_proportion": float(mismatch_prop),
                "match_count": int(match_count),
                "mismatch_count": int(mismatch_count),
            })
    pairs = set()
    for row in rows:
        pair = tuple(sorted((row["sample1"], row["sample2"])))
        if pair in pairs or pair[0] == pair[1]:
            raise WorkflowError("Duplicate or self comparison in SKA distance table.")
        if not math.isfinite(row["snp_distance"]) or row["snp_distance"] < 0 or not math.isfinite(row["mismatch_proportion"]) or not 0 <= row["mismatch_proportion"] <= 1 or row["match_count"] < 0 or row["mismatch_count"] < 0:
            raise WorkflowError("Nonfinite or invalid values in SKA distance table.")
        pairs.add(pair)
    return rows


def run_ska(genomes: dict[str, str | list[str]], output_dir: Path, kmer_size: int) -> dict[str, Any]:
    if len(genomes) < 2:
        raise WorkflowError(
            "ska2 SNP resolution requires at least two genomes (query plus one representative)."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    ska = require_executable("ska")
    file_list = output_dir / "ska_input.tsv"
    build_file_list(genomes, file_list)
    skf_file = output_dir / "merged.skf"
    build_command = [
        ska, "build", "-f", str(file_list), "-o", str(output_dir / "merged"), "-k", str(kmer_size),
    ]
    if any(not isinstance(value, str) for value in genomes.values()):
        filters = load_json(CONFIG_DIR / "workflow.json")["snp_resolution"]["read_filters"]
        build_command.extend(["--min-count", str(filters["min_count"]), "--min-qual", str(filters["min_qual"]), "--qual-filter", filters["qual_filter"]])
    returncode = run_logged(
        build_command, output_dir, output_dir / "ska_build.stdout.log", output_dir / "ska_build.stderr.log"
    )
    if returncode != 0 or not skf_file.is_file():
        raise WorkflowError(f"ska build failed; inspect {output_dir / 'ska_build.stderr.log'}.")
    distance_path = output_dir / "distances.tsv"
    distance_command = [ska, "distance", "-o", str(distance_path), str(skf_file)]
    returncode = run_logged(
        distance_command, output_dir,
        output_dir / "ska_distance.stdout.log", output_dir / "ska_distance.stderr.log",
    )
    if returncode != 0 or not distance_path.is_file():
        raise WorkflowError(f"ska distance failed; inspect {output_dir / 'ska_distance.stderr.log'}.")
    rows = parse_distance_table(distance_path)
    observed = {tuple(sorted((row["sample1"], row["sample2"]))) for row in rows}
    if observed != set(combinations(sorted(genomes), 2)):
        raise WorkflowError("SKA distance matrix does not contain exactly the expected sample pairs.")
    result = {
        "status": "PASS",
        "commands": [build_command, distance_command],
        "distance_table": str(distance_path),
        "distances": rows,
    }
    write_json(output_dir / "run_ska.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genome", action="append", default=[], metavar="SAMPLE_ID=PATH",
        help="Repeatable sample_id=path pair; one entry must use the sample id QUERY.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--kmer-size", type=int, default=31)
    args = parser.parse_args()
    genomes = dict(entry.split("=", 1) for entry in args.genome)
    try:
        run_ska(genomes, Path(args.output_dir), args.kmer_size)
    except WorkflowError as error:
        write_json(Path(args.output_dir) / "run_ska.json", {"status": "FAIL", "error": str(error)})
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
