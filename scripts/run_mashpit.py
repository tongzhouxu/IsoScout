#!/usr/bin/env python3
"""Validate a database and execute the fixed Mashpit query profile."""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

from common import (
    CONFIG_DIR,
    WorkflowError,
    load_json,
    require_executable,
    run_logged,
    sha256_file,
    write_json,
)


def validate_database(database_dir: Path, expected_name: str) -> dict:
    metadata_path = database_dir / "database.json"
    if not database_dir.is_dir() or not metadata_path.is_file():
        raise WorkflowError(
            f"Mashpit database '{expected_name}' is missing or has no database.json: {database_dir}"
        )
    metadata = load_json(metadata_path)
    required = ("name", "version", "build_date", "checksum", "source")
    missing = [key for key in required if not metadata.get(key)]
    if missing:
        raise WorkflowError("Database metadata is incomplete: " + ", ".join(missing))
    if metadata["name"] != expected_name:
        raise WorkflowError(
            f"Database name mismatch: routed '{expected_name}', metadata says '{metadata['name']}'."
        )
    artifact_name = metadata.get("checksum_file")
    if artifact_name:
        artifact = database_dir / artifact_name
        if not artifact.is_file():
            raise WorkflowError(f"Database checksum target is missing: {artifact}")
        observed = sha256_file(artifact)
        if observed != metadata["checksum"]:
            raise WorkflowError("Database checksum verification failed.")
    db_files = list(database_dir.glob("*.db"))
    sig_files = list(database_dir.glob("*.sig"))
    if len(db_files) != 1 or len(sig_files) != 1:
        raise WorkflowError("Mashpit database must contain exactly one .db and one .sig entry.")
    if db_files[0].stem != sig_files[0].stem:
        raise WorkflowError("Mashpit .db and .sig basenames do not match.")
    signature_checksum = metadata.get("signature_checksum")
    if signature_checksum and sha256_file(sig_files[0]) != signature_checksum:
        raise WorkflowError("Mashpit signature checksum verification failed.")
    try:
        with sqlite3.connect(db_files[0]) as connection:
            rows = dict(connection.execute("SELECT name, value FROM DESC").fetchall())
        if rows.get("Version") and rows["Version"] != metadata["version"]:
            raise WorkflowError("Database metadata release does not match its Mashpit DESC table.")
        metadata["mashpit_database_settings"] = {
            "type": rows.get("Type"),
            "species": rows.get("Species"),
            "hash_number": int(rows["Hash_number"]),
            "kmer_size": int(rows["Kmer_size"]),
        }
    except (sqlite3.Error, KeyError, TypeError, ValueError) as error:
        raise WorkflowError(f"Cannot read Mashpit database settings: {error}") from error
    metadata["database_directory"] = str(database_dir.resolve())
    if "membership_snapshot" in metadata:
        from expand_cluster_members import bundled_snapshot
        bundled_snapshot(metadata)
    return metadata


def run_mashpit(assembly: Path, database_dir: Path, expected_name: str, output_dir: Path) -> dict:
    if not assembly.is_file():
        raise WorkflowError(f"Assembly is missing: {assembly}")
    metadata = validate_database(database_dir, expected_name)
    workflow = load_json(CONFIG_DIR / "workflow.json")
    profile = workflow["mashpit"]
    mashpit = require_executable("mashpit")
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        mashpit,
        "query",
        str(assembly),
        str(database_dir),
        "--number",
        str(profile["number"]),
        "--threshold",
        str(profile["threshold"]),
        "--tie-tolerance-hashes",
        str(profile["tie_tolerance_hashes"]),
    ]
    returncode = run_logged(
        command, output_dir, output_dir / "mashpit.stdout.log", output_dir / "mashpit.stderr.log"
    )
    metadata_out = {
        "status": "PASS" if returncode == 0 else "FAIL",
        "returncode": returncode,
        "command": command,
        "profile": profile["profile"],
        "database": metadata,
        "output_directory": str(output_dir),
    }
    write_json(output_dir / "mashpit_run.json", metadata_out)
    if returncode:
        raise WorkflowError(f"Mashpit query failed; inspect {output_dir / 'mashpit.stderr.log'}.")
    return metadata_out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly")
    parser.add_argument("database")
    parser.add_argument("output_directory")
    parser.add_argument("--database-name", required=True)
    args = parser.parse_args()
    try:
        run_mashpit(
            Path(args.assembly), Path(args.database), args.database_name, Path(args.output_directory)
        )
    except WorkflowError as error:
        failure = {"status": "FAIL", "error": str(error)}
        write_json(Path(args.output_directory) / "mashpit_run.json", failure)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
