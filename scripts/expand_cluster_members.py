#!/usr/bin/env python3
"""Retrieve membership from the exact NCBI release used by a Mashpit database."""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import urllib.request
from pathlib import Path
from typing import Any

from common import CONFIG_DIR, WorkflowError, load_json, sha256_file, utc_now, write_json

ACCESSION = re.compile(r"GC[AF]_\d+\.\d+\Z")
RELEASE = re.compile(r"PDG\d+\.\d+\Z")


def download_file(url: str, destination: Path, max_bytes: int) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "IsoScout/1.5 release-membership"})
    size = 0
    with urllib.request.urlopen(request, timeout=60) as response:
        # Do not silently follow an archived-release URL to latest or another host.
        if response.geturl() != url:
            raise WorkflowError("NCBI release download redirected; refusing a different snapshot.")
        with partial.open("wb") as handle:
            while chunk := response.read(1024 * 1024):
                size += len(chunk)
                if size > max_bytes:
                    raise WorkflowError("Release metadata exceeds the configured download size ceiling.")
                handle.write(chunk)
    partial.replace(destination)
    return {"url": url, "path": str(destination), "sha256": sha256_file(destination), "size_bytes": size}


def _table(path: Path, required: set[str]):
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, mode="rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise WorkflowError(f"Missing required columns in {path}: {sorted(required)}")
        yield from reader


def snapshot_urls(species: str, release: str) -> dict[str, str]:
    if not isinstance(release, str) or not isinstance(species, str) or not RELEASE.fullmatch(release) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", species):
        raise WorkflowError("Pinned PDG release and NCBI species-directory metadata are required for expansion.")
    base = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{species}/{release}"
    return {
        "metadata": f"{base}/Metadata/{release}.metadata.tsv",
        "all_isolates": f"{base}/Clusters/{release}.reference_target.all_isolates.tsv",
    }


def bundled_snapshot(database: dict[str, Any]) -> tuple[Path, Path, list[dict[str, Any]]] | None:
    """Verify a packaged snapshot before use; never fall back after corruption."""
    snapshot = database.get("membership_snapshot")
    if snapshot is None:
        return None
    if not isinstance(snapshot, dict):
        raise WorkflowError("Bundled membership snapshot manifest is malformed.")
    release = database.get("version", "")
    species = database.get("mashpit_database_settings", {}).get("species") or database.get("species", "")
    urls = snapshot_urls(species, release)
    if snapshot.get("schema_version") != "1.0.0" or snapshot.get("release") != release or snapshot.get("species") != species:
        raise WorkflowError("Bundled membership snapshot does not match the database release and species.")
    directory = database.get("database_directory")
    if not directory:
        raise WorkflowError("Bundled membership snapshot needs its database directory.")
    root = Path(directory).resolve()
    paths = []
    sources = []
    for key, filename in (("metadata", f"{release}.metadata.tsv.gz"),
                          ("all_isolates", f"{release}.reference_target.all_isolates.tsv.gz")):
        record = snapshot.get(key)
        if not isinstance(record, dict) or record.get("path") != f"membership/{filename}" or record.get("source_url") != urls[key]:
            raise WorkflowError(f"Bundled {key} manifest is missing or points to the wrong release.")
        path = (root / record["path"]).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise WorkflowError(f"Bundled {key} table is missing or outside its database directory.")
        if path.stat().st_size != record.get("size_bytes") or sha256_file(path) != record.get("sha256"):
            raise WorkflowError(f"Bundled {key} table failed checksum or size verification.")
        paths.append(path)
        sources.append({"url": urls[key], "path": str(path), "sha256": record["sha256"],
                        "size_bytes": record["size_bytes"], "original_sha256": record.get("original_sha256"),
                        "origin": "bundled_database_release"})
    return paths[0], paths[1], sources


def parse_members(metadata: Path, isolates: Path, clusters: list[str]) -> dict[str, Any]:
    wanted = set(clusters)
    target_clusters: dict[str, str] = {}
    for row in _table(isolates, {"target_acc", "PDS_acc"}):
        target, cluster = row["target_acc"], row["PDS_acc"]
        if cluster not in wanted:
            continue
        if target in target_clusters and target_clusters[target] != cluster:
            raise WorkflowError(f"Conflicting cluster membership for target {target}.")
        target_clusters[target] = cluster
    members: dict[str, dict[str, str]] = {}
    seen_targets = set()
    without_assembly = set()
    for row in _table(metadata, {"target_acc", "asm_acc"}):
        target = row["target_acc"]
        if target not in target_clusters:
            continue
        seen_targets.add(target)
        accession = (row["asm_acc"] or "").strip()
        if not ACCESSION.fullmatch(accession):
            without_assembly.add(target)
            continue
        cluster = target_clusters[target]
        if accession in members and members[accession]["cluster"] != cluster:
            raise WorkflowError(f"Assembly {accession} is assigned to multiple requested clusters.")
        members[accession] = {"accession": accession, "cluster": cluster, "source": "release_cluster_member"}
    return {
        "members": sorted(members.values(), key=lambda item: (item["cluster"], item["accession"])),
        "requested_clusters": sorted(wanted), "isolate_records": len(target_clusters),
        "targets_without_metadata": sorted(set(target_clusters) - seen_targets),
        "targets_without_assembly": sorted(without_assembly),
        "clusters_without_members": sorted(wanted - {row["cluster"] for row in members.values()}),
    }


def discover_members(database: dict[str, Any], clusters: list[str], output_dir: Path, policy: dict[str, Any]) -> dict[str, Any]:
    release = database.get("version", "")
    species = database.get("species") or database.get("mashpit_database_settings", {}).get("species", "")
    result: dict[str, Any] = {"schema_version": "1.0.0", "release": release, "species": species, "checked_at": utc_now(), "sources": [], "members": []}
    try:
        bundled = bundled_snapshot(database)
        if bundled is not None:
            metadata, isolates, result["sources"] = bundled
            result["source_mode"] = "bundled_database_release"
        else:
            urls = snapshot_urls(species, release)
            metadata = output_dir / f"{release}.metadata.tsv"
            isolates = output_dir / f"{release}.reference_target.all_isolates.tsv"
            for key, path in (("metadata", metadata), ("all_isolates", isolates)):
                result["sources"].append(download_file(urls[key], path, policy["max_metadata_file_bytes"]))
            result["source_mode"] = "ncbi_exact_release"
            result["retrieved_at"] = utc_now()
        result.update(parse_members(metadata, isolates, clusters))
        result["status"] = "AVAILABLE"
    except (WorkflowError, OSError, ValueError, TypeError, csv.Error) as error:
        result.update(status="UNAVAILABLE", error=str(error))
    write_json(output_dir / "membership.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-metadata", required=True)
    parser.add_argument("--cluster", action="append", required=True)
    parser.add_argument("--policy", default=str(CONFIG_DIR / "refinement-policy.json"))
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    if output.exists():
        parser.error("Choose a new output directory to retain immutable provenance.")
    policy = load_json(Path(args.policy))
    database = load_json(Path(args.database_metadata))
    database.setdefault("database_directory", str(Path(args.database_metadata).resolve().parent))
    result = discover_members(database, args.cluster, output, policy.get("expansion", policy))
    return 0 if result["status"] == "AVAILABLE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
