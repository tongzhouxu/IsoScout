#!/usr/bin/env python3
"""Download representative assembly FASTA files from NCBI for SNP resolution.

Mirrors the dehydrated-download-then-rehydrate pattern Mashpit itself uses to
build taxon databases (see upstream ``mashpit build``), since a Mashpit
database only retains sourmash signatures, not the representative assemblies.
"""

from __future__ import annotations

import argparse
import csv
import re
import time
import zipfile
from pathlib import Path
from typing import Any

from common import WorkflowError, require_executable, run_logged, sha256_file, write_json


def read_local_manifest(path: Path) -> dict[str, dict[str, str]]:
    """Explicit accession/path/SHA-256 mapping; paths are relative to the manifest."""
    result = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"accession", "path", "sha256"} <= set(reader.fieldnames or []):
            raise WorkflowError("Assembly manifest requires accession, path and sha256 columns")
        for row in reader:
            accession = row["accession"]
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", accession) or accession == "QUERY" or accession in result:
                raise WorkflowError("Invalid, reserved or duplicate assembly manifest accession")
            if not re.fullmatch(r"[0-9a-f]{64}",row["sha256"]):
                raise WorkflowError("Invalid assembly manifest SHA-256")
            assembly = Path(row["path"]).expanduser()
            if not assembly.is_absolute(): assembly = path.resolve().parent/assembly
            result[accession] = {"path":str(assembly.resolve()),"sha256":row["sha256"]}
    return result


def locate_fasta(package_dir: Path, accession: str) -> Path | None:
    candidates = [
        path
        for path in (package_dir / "ncbi_dataset" / "data" / accession).glob("*_genomic.fna")
        if path.is_file() and path.stat().st_size > 0
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_size)


def download_batch(datasets_exe: str, accessions: list[str], batch_dir: Path) -> tuple[list[str], list[str]]:
    batch_dir.mkdir(parents=True, exist_ok=True)
    accession_file = batch_dir / "accessions.txt"
    accession_file.write_text("\n".join(accessions) + "\n", encoding="utf-8")
    package = batch_dir / "package.zip"
    download_command = [
        datasets_exe, "download", "genome", "accession",
        "--inputfile", str(accession_file),
        "--include", "genome",
        "--dehydrated",
        "--filename", str(package),
    ]
    returncode = run_logged(
        download_command, batch_dir,
        batch_dir / "datasets_download.stdout.log", batch_dir / "datasets_download.stderr.log",
    )
    if returncode != 0 or not package.is_file():
        raise WorkflowError(
            f"NCBI datasets download failed; inspect {batch_dir / 'datasets_download.stderr.log'}."
        )
    with zipfile.ZipFile(str(package), "r") as archive:
        archive.extractall(str(batch_dir))
    rehydrate_command = [datasets_exe, "rehydrate", "--directory", str(batch_dir)]
    returncode = run_logged(
        rehydrate_command, batch_dir,
        batch_dir / "datasets_rehydrate.stdout.log", batch_dir / "datasets_rehydrate.stderr.log",
    )
    if returncode != 0:
        raise WorkflowError(
            f"NCBI datasets rehydrate failed; inspect {batch_dir / 'datasets_rehydrate.stderr.log'}."
        )
    return download_command, rehydrate_command


def fetch_genomes(
    accessions: list[str], output_dir: Path, attempts: int, retry_delay_seconds: float,
    local_manifest: Path | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    if not accessions:
        result = {"status": "COMPLETE", "verified": {}, "unavailable": [], "errors": {}, "commands": []}
        write_json(output_dir / "fetch_reference_genomes.json", result)
        return result
    verified: dict[str, str] = {}
    errors: dict[str, str] = {}
    commands: list[list[str]] = []
    pending = sorted(set(accessions))
    local_records = read_local_manifest(local_manifest) if local_manifest else {}
    local_used = []
    for accession in pending:
        if accession not in local_records: continue
        row = local_records[accession]; path = Path(row["path"])
        if not path.is_file() or sha256_file(path) != row["sha256"]:
            raise WorkflowError(f"Local assembly missing or checksum mismatch: {accession}; no download substitution made")
        verified[accession] = str(path)
        local_used.append(accession)
    pending = [accession for accession in pending if accession not in verified]
    datasets_exe = require_executable("datasets") if pending else None
    for attempt in range(1, attempts + 1):
        if not pending:
            break
        attempt_dir = output_dir / f"attempt_{attempt:02d}"
        try:
            download_command, rehydrate_command = download_batch(datasets_exe, pending, attempt_dir)
            commands.extend([download_command, rehydrate_command])
        except WorkflowError as error:
            for accession in pending:
                errors[accession] = str(error)
            if attempt < attempts:
                time.sleep(retry_delay_seconds)
            continue
        still_pending = []
        for accession in pending:
            fasta = locate_fasta(attempt_dir, accession)
            if fasta:
                verified[accession] = str(fasta)
            else:
                still_pending.append(accession)
                errors[accession] = "genomic FASTA not found after download"
        pending = still_pending
        if pending and attempt < attempts:
            time.sleep(retry_delay_seconds)
    result = {
        "status": "COMPLETE" if not pending else ("PARTIAL" if verified else "FAILED"),
        "verified": verified,
        "unavailable": sorted(pending),
        "errors": errors,
        "commands": commands,
        "local_manifest": str(local_manifest) if local_manifest else None,
        "local_manifest_sha256": sha256_file(local_manifest) if local_manifest else None,
        "local_assemblies": local_used,
    }
    write_json(output_dir / "fetch_reference_genomes.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", action="append", default=[])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--attempts", type=int, default=3)
    parser.add_argument("--retry-delay-seconds", type=float, default=5.0)
    args = parser.parse_args()
    result = fetch_genomes(
        args.accession, Path(args.output_dir), args.attempts, args.retry_delay_seconds
    )
    return 0 if result["status"] != "FAILED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
