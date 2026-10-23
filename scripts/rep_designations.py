#!/usr/bin/env python3
"""Look up CDC REP designations for called NCBI SNP clusters.

Database packages from the `databases-rep-2026-10` release add a REP_CLUSTER
table to the Mashpit .db (PDS_acc primary key, rep_designation, organism,
anchor_pnusa, anchor_pdt, evidence_status, source_url). Older packages
(databases-v2) have no such table. Matching is by exact PDS accession only:
no version stripping, prefix matching, or fuzzy lookup.

A missing table or a missing row means "no REP label in this database", not
"not a REP strain": REP tables list only the clusters CDC has publicly linked
to a REP strain at the time the release was built.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

from common import WorkflowError


REP_TABLE = "REP_CLUSTER"
REP_FIELDS = ("rep_designation", "organism", "anchor_pnusa", "anchor_pdt", "evidence_status", "source_url")
NO_LABEL = "No REP label in this database"
MEMBERSHIP_CAVEAT = (
    "A REP designation means CDC has linked this NCBI SNP cluster to a REP strain in this database "
    "release. It does not confirm that this isolate is a member of the outbreak or REP strain; that "
    "requires validated typing and public health investigation."
)
ABSENCE_CAVEAT = (
    "A cluster without a REP label may still be linked to a REP strain; the table lists only "
    "clusters publicly linked to a REP strain when the release was built."
)


def database_file(database_metadata: dict[str, Any]) -> Path:
    directory = Path(database_metadata["database_directory"])
    matches = list(directory.glob("*.db"))
    if len(matches) != 1:
        raise WorkflowError("Mashpit database must contain exactly one .db file.")
    return matches[0]


def read_rep_table(path: Path, clusters: list[str]) -> dict[str, dict[str, Any]] | None:
    """Return exact-match REP rows keyed by PDS accession, or None when the table is absent."""
    connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        present = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?", (REP_TABLE,)
        ).fetchone()
        if not present:
            return None
        rows: dict[str, dict[str, Any]] = {}
        for cluster in dict.fromkeys(clusters):
            row = connection.execute(
                f"SELECT {', '.join(REP_FIELDS)} FROM {REP_TABLE} WHERE PDS_acc = ?", (cluster,)
            ).fetchone()
            if row is not None:
                rows[cluster] = dict(zip(REP_FIELDS, row))
        return rows
    finally:
        connection.close()


def final_called_clusters(result: dict[str, Any]) -> tuple[str | None, list[str]]:
    """The clusters the screen finally calls, and the evidence they come from."""
    snp = result.get("snp_resolution") or {}
    if snp.get("status") in {"COMPARED", "AMBIGUOUS", "RESOLVED"} and snp.get("nearest_clusters"):
        if snp.get("nearest_cluster"):
            return "snp_resolution", [snp["nearest_cluster"]]
        return "snp_resolution_unresolved", list(snp.get("cluster_candidates") or snp["nearest_clusters"])
    match = result.get("mashpit_result") or {}
    if match.get("best_candidate") and not match.get("below_threshold"):
        return "mashpit_top_candidate", [match["best_candidate"]["cluster"]]
    return None, []


def _entry(cluster: str, rows: dict[str, dict[str, Any]] | None) -> dict[str, Any]:
    row = (rows or {}).get(cluster)
    if row is None:
        return {"cluster": cluster, "status": "NO_LABEL", "rep_designation": None, "source_url": None,
                "statement": NO_LABEL}
    return {"cluster": cluster, "status": "LABELLED", **row}


def annotate(result: dict[str, Any], database_metadata: dict[str, Any]) -> dict[str, Any]:
    """Build the result.json `rep_designations` block from the final call and Mashpit candidates."""
    basis, final_clusters = final_called_clusters(result)
    match = result.get("mashpit_result") or {}
    candidates = [row["cluster"] for row in
                  ([match["best_candidate"]] if match.get("best_candidate") else [])
                  + (match.get("alternative_candidates") or [])]
    annotation: dict[str, Any] = {
        "table": REP_TABLE, "matching": "exact_pds_accession",
        "final_call_basis": basis,
    }
    try:
        path = database_file(database_metadata)
        rows = read_rep_table(path, final_clusters + candidates)
    except (WorkflowError, sqlite3.Error, KeyError, OSError) as error:
        annotation.update({
            "status": "UNREADABLE", "table_present": None, "error": str(error),
            "final_clusters": [{"cluster": cluster, "status": "UNKNOWN"} for cluster in final_clusters],
            "candidate_clusters": [{"cluster": cluster, "status": "UNKNOWN"} for cluster in candidates],
            "statement": "REP labels could not be read from this database; no REP status is reported.",
        })
        return annotation
    annotation["database_file"] = str(path)
    annotation["table_present"] = rows is not None
    annotation["status"] = "AVAILABLE" if rows is not None else "TABLE_ABSENT"
    annotation["final_clusters"] = [_entry(cluster, rows) for cluster in final_clusters]
    annotation["candidate_clusters"] = [_entry(cluster, rows) for cluster in candidates]
    if rows is None:
        annotation["statement"] = (
            NO_LABEL + f": this database release has no {REP_TABLE} table. " + ABSENCE_CAVEAT
        )
    else:
        annotation["statement"] = MEMBERSHIP_CAVEAT + " " + ABSENCE_CAVEAT
    return annotation


def describe(entry: dict[str, Any]) -> str:
    if entry["status"] == "LABELLED":
        return f"{entry['cluster']}: {entry['rep_designation']} ({entry['source_url']})"
    if entry["status"] == "UNKNOWN":
        return f"{entry['cluster']}: REP label unreadable"
    return f"{entry['cluster']}: no REP label in this database"


def summary_line(annotation: dict[str, Any]) -> str:
    if annotation["status"] == "UNREADABLE":
        return "CDC REP designation: could not be read from this database."
    if annotation["status"] == "TABLE_ABSENT":
        return f"CDC REP designation: no REP label in this database (no {REP_TABLE} table)."
    finals = annotation["final_clusters"]
    if not finals:
        return "CDC REP designation: no final cluster was called."
    line = "CDC REP designation of called cluster(s): " + "; ".join(describe(entry) for entry in finals) + "."
    if any(entry["status"] == "LABELLED" for entry in finals):
        line += " This does not confirm outbreak-strain membership."
    return line
