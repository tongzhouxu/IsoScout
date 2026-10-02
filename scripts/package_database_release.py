#!/usr/bin/env python3
"""Add exact-release membership tables to an existing Mashpit database package."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sqlite3
import tarfile
import tempfile
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, WorkflowError, load_json, sha256_file, write_json
from expand_cluster_members import bundled_snapshot, download_file, snapshot_urls
from run_mashpit import validate_database


def verify_representatives(db_path: Path, metadata_path: Path, isolates_path: Path) -> int:
    """Require the release tables to agree with every sketched representative."""
    with sqlite3.connect(db_path.as_uri() + "?mode=ro", uri=True) as connection:
        expected = {target: (cluster, assembly) for target, cluster, assembly in connection.execute(
            "SELECT PDT_acc, PDS_acc, asm_acc FROM REPRESENTATIVE"
        )}
    if not expected:
        raise WorkflowError("Database contains no representatives to verify against the release tables.")
    found_cluster: dict[str, str] = {}
    with isolates_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"target_acc", "PDS_acc"}.issubset(reader.fieldnames or []):
            raise WorkflowError("Membership table lacks target_acc or PDS_acc.")
        for row in reader:
            target = row["target_acc"]
            if target in expected:
                cluster = row["PDS_acc"]
                if target in found_cluster and found_cluster[target] != cluster:
                    raise WorkflowError(f"Conflicting memberships for representative {target}.")
                found_cluster[target] = cluster
    found_assembly: dict[str, str] = {}
    with metadata_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"target_acc", "asm_acc"}.issubset(reader.fieldnames or []):
            raise WorkflowError("Metadata table lacks target_acc or asm_acc.")
        for row in reader:
            target = row["target_acc"]
            if target in expected:
                assembly = row["asm_acc"]
                if target in found_assembly and found_assembly[target] != assembly:
                    raise WorkflowError(f"Conflicting assemblies for representative {target}.")
                found_assembly[target] = assembly
    for target, (cluster, assembly) in expected.items():
        if (found_cluster.get(target), found_assembly.get(target)) != (cluster, assembly):
            raise WorkflowError(f"Release tables disagree with database representative {target}.")
    return len(expected)


def compress_table(source: Path, destination: Path, url: str) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as original, destination.open("wb") as output:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=6, mtime=0, fileobj=output) as packed:
            shutil.copyfileobj(original, packed, length=1024 * 1024)
    return {"path": str(destination.parent.name + "/" + destination.name),
            "sha256": sha256_file(destination), "size_bytes": destination.stat().st_size,
            "source_url": url, "original_sha256": sha256_file(source),
            "original_size_bytes": source.stat().st_size}


def create_archive(package: Path, destination: Path) -> None:
    """Create a deterministic tar.gz asset from the validated package."""
    if destination.exists():
        raise WorkflowError(f"Archive already exists: {destination}")
    partial = destination.with_name(destination.name + ".partial")
    try:
        with partial.open("wb") as output:
            with gzip.GzipFile(filename="", mode="wb", compresslevel=6, mtime=0, fileobj=output) as zipped:
                with tarfile.open(fileobj=zipped, mode="w|") as archive:
                    for path in sorted(package.rglob("*")):
                        if not path.is_file():
                            continue
                        name = f"{package.name}/{path.relative_to(package).as_posix()}"
                        info = archive.gettarinfo(str(path), arcname=name)
                        info.uid = info.gid = 0
                        info.uname = info.gname = ""
                        info.mtime = 0
                        with path.open("rb") as handle:
                            archive.addfile(info, handle)
        partial.replace(destination)
    finally:
        partial.unlink(missing_ok=True)


def package_database(database_dir: Path, output_dir: Path, policy: dict[str, Any],
                     metadata_source: Path | None = None, isolates_source: Path | None = None,
                     archive: Path | None = None) -> dict[str, Any]:
    database_dir = database_dir.resolve()
    output_dir = output_dir.resolve()
    archive = archive.resolve() if archive else None
    if output_dir.exists() or (archive and archive.exists()):
        raise WorkflowError("Package directory and archive must not already exist.")
    if output_dir.is_relative_to(database_dir):
        raise WorkflowError("Package output must not be inside the source database.")
    if (metadata_source is None) != (isolates_source is None):
        raise WorkflowError("Supply both local release tables or allow both to be downloaded.")
    database = validate_database(database_dir, database_dir.name)
    if database.get("membership_snapshot"):
        raise WorkflowError("Source database already has a membership snapshot.")
    species = database["mashpit_database_settings"]["species"]
    release = database["version"]
    urls = snapshot_urls(species, release)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=output_dir.name + ".", dir=output_dir.parent) as temporary:
        root = Path(temporary)
        staging = root / output_dir.name
        shutil.copytree(database_dir, staging)
        sources: dict[str, Path] = {}
        for key, supplied in (("metadata", metadata_source), ("all_isolates", isolates_source)):
            if supplied is None:
                source = root / (f"{release}.metadata.tsv" if key == "metadata" else f"{release}.reference_target.all_isolates.tsv")
                download_file(urls[key], source, policy["max_metadata_file_bytes"])
            else:
                source = supplied.resolve()
                if not source.is_file() or source.stat().st_size > policy["max_metadata_file_bytes"]:
                    raise WorkflowError(f"Supplied {key} table is missing or exceeds the size ceiling.")
            sources[key] = source
        verified = verify_representatives(staging / f"{database['name']}.db", sources["metadata"], sources["all_isolates"])
        manifest = {"schema_version": "1.0.0", "release": release, "species": species,
                    "acquisition": "supplied_local_files" if metadata_source else "ncbi_exact_release",
                    "verified_representatives": verified}
        for key, source in sources.items():
            filename = f"{release}.metadata.tsv.gz" if key == "metadata" else f"{release}.reference_target.all_isolates.tsv.gz"
            manifest[key] = compress_table(source, staging / "membership" / filename, urls[key])
        persisted = load_json(staging / "database.json")
        persisted["signature_checksum"] = sha256_file(staging / f"{database['name']}.sig")
        persisted["membership_snapshot"] = manifest
        write_json(staging / "database.json", persisted)
        validate_database(staging, database["name"])
        bundled_snapshot({**persisted, "database_directory": str(staging), "mashpit_database_settings": database["mashpit_database_settings"]})
        staging.replace(output_dir)
    if archive:
        archive.parent.mkdir(parents=True, exist_ok=True)
        create_archive(output_dir, archive)
    return {"database_directory": str(output_dir), "archive": str(archive) if archive else None,
            "membership_snapshot": manifest}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--archive")
    parser.add_argument("--metadata-tsv")
    parser.add_argument("--all-isolates-tsv")
    args = parser.parse_args()
    try:
        result = package_database(
            Path(args.database_dir), Path(args.output_dir),
            load_json(CONFIG_DIR / "refinement-policy.json")["expansion"],
            Path(args.metadata_tsv) if args.metadata_tsv else None,
            Path(args.all_isolates_tsv) if args.all_isolates_tsv else None,
            Path(args.archive) if args.archive else None,
        )
    except (WorkflowError, OSError, ValueError, sqlite3.Error, csv.Error) as error:
        parser.exit(2, f"Database packaging failed: {error}\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
