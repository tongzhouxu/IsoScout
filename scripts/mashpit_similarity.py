"""Shared validated representative data and Mashpit sketch-tolerance arithmetic."""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from common import WorkflowError, load_json
from parse_mashpit_results import CLUSTER_FIELDS, SCORE_FIELDS, first_value

ACCESSION_FIELDS = ("asm_acc", "assembly_accession", "assembly_acc")
BIOSAMPLE_FIELDS = ("biosample_acc", "biosample_accession")


def within_tolerance(gap: float, tolerance: float) -> bool:
    # Account only for floating-point roundoff at an inclusive boundary.
    return gap <= tolerance or math.isclose(gap, tolerance, rel_tol=1e-12, abs_tol=1e-15)


def read_scores(path: Path) -> tuple[list[dict[str, Any]], int]:
    """Fail visibly on unusable rows; never silently distort a distribution."""
    by_accession: dict[str, dict[str, Any]] = {}
    duplicate_rows = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        for aliases in (ACCESSION_FIELDS, CLUSTER_FIELDS, SCORE_FIELDS):
            if not any(field in fields for field in aliases):
                raise WorkflowError(f"Missing required representative column in {path}: {aliases}")
        for line, row in enumerate(reader, 2):
            accession = first_value(row, ACCESSION_FIELDS)
            cluster = first_value(row, CLUSTER_FIELDS)
            try:
                score = float(first_value(row, SCORE_FIELDS))
            except (TypeError, ValueError) as error:
                raise WorkflowError(f"Invalid similarity score at {path}:{line}") from error
            if not accession or not cluster or not math.isfinite(score) or not 0 <= score <= 1:
                raise WorkflowError(f"Invalid representative at {path}:{line}")
            item = {"accession": str(accession), "cluster": str(cluster), "score": score,
                    "biosample": first_value(row, BIOSAMPLE_FIELDS)}
            if accession in by_accession:
                if {key: value for key, value in by_accession[accession].items() if key != "source_rows"} != item:
                    raise WorkflowError(f"Conflicting scores or clusters for accession {accession}")
                duplicate_rows += 1
                by_accession[accession]["source_rows"].append(line)
            else:
                item["source_rows"] = [line]
                by_accession[accession] = item
    return sorted(by_accession.values(), key=lambda item: (-item["score"], item["accession"])), duplicate_rows


def recorded_query_context(output_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Use the historical run's settings, never guessed sketch sizes or flags."""
    run = load_json(output_dir / "mashpit_run.json")
    try:
        command = run["command"]
        profile = {key: convert(command[command.index(flag) + 1]) for key, flag, convert in (
            ("number", "--number", int), ("threshold", "--threshold", float),
            ("tie_tolerance_hashes", "--tie-tolerance-hashes", float),
        )}
        return run["database"], profile
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise WorkflowError("Recorded Mashpit query settings are missing or malformed.") from error


def score_tolerance(database: dict[str, Any], profile: dict[str, Any]) -> float:
    size = database.get("mashpit_database_settings", {}).get("hash_number")
    hashes = profile.get("tie_tolerance_hashes")
    if not isinstance(size, int) or isinstance(size, bool) or size <= 0:
        raise WorkflowError("A positive database hash_number is required; no sketch size is assumed.")
    if not isinstance(hashes, (int, float)) or isinstance(hashes, bool) or not math.isfinite(hashes) or hashes < 0:
        raise WorkflowError("Mashpit hash tolerance must be finite and nonnegative.")
    return hashes / size
