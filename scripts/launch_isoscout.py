#!/usr/bin/env python3
"""Run IsoScout's container from the host with reproducible path mounts."""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from pathlib import Path


DEFAULT_IMAGE = "isoscout:local"
ORGANISMS = ("salmonella", "ecoli_shigella", "listeria", "campylobacter", "cronobacter")


def git_identity(root: Path) -> tuple[str, str]:
    # Installed wheels may live inside an unrelated Git checkout. Only report
    # a launcher commit when this file is the checkout's actual source script.
    if (root / "scripts" / "launch_isoscout.py").resolve() != Path(__file__).resolve():
        return "unknown", "unknown"
    try:
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        ).stdout.strip()
        dirty = "true" if subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            check=True, capture_output=True, text=True,
        ).stdout.strip() else "false"
        return commit, dirty
    except (OSError, subprocess.CalledProcessError):
        return "unknown", "unknown"


def image_identity(image: str) -> tuple[str, str]:
    try:
        value = subprocess.run(
            [
                "docker", "image", "inspect", image, "--format",
                '{{.Id}}\t{{index .Config.Labels "org.opencontainers.image.revision"}}',
            ],
            check=True, capture_output=True, text=True,
        ).stdout.strip()
        image_id, _, revision = value.partition("\t")
        return image_id, revision or "unknown"
    except FileNotFoundError as error:
        raise SystemExit("Docker was not found. Install and start Docker, then retry.") from error
    except subprocess.CalledProcessError as error:
        message = error.stderr.strip() or f"Docker image {image!r} is unavailable"
        raise SystemExit(message) from error


def add_same_path_mount(command: list[str], path: Path, mounted: set[tuple[str, str]]) -> Path:
    """Mount a config/manifest parent at its absolute host path in the container."""
    resolved = path.expanduser().resolve()
    parent = resolved if resolved.is_dir() else resolved.parent
    key = (str(parent), "ro")
    if key not in mounted:
        command.extend(["--volume", f"{parent}:{parent}:ro"])
        mounted.add(key)
    return resolved


def mount_manifest_assemblies(command: list[str], manifest: Path, mounted: set[tuple[str, str]]) -> None:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "path" not in (reader.fieldnames or []):
            raise SystemExit("Assembly manifest requires a path column.")
        for row in reader:
            assembly = Path(row["path"]).expanduser()
            if not assembly.is_absolute():
                assembly = manifest.parent / assembly
            assembly = assembly.resolve()
            if not assembly.is_file():
                raise SystemExit(f"Manifest assembly not found: {assembly}")
            add_same_path_mount(command, assembly, mounted)


def build_command(
    args: argparse.Namespace, root: Path, runtime_image_id: str, image_revision: str,
) -> list[str]:
    inputs = [Path(value).expanduser().resolve() for value in args.inputs]
    if len(inputs) not in (1, 2):
        raise SystemExit("Provide one assembly or two paired-end FASTQ files.")
    missing = [str(path) for path in inputs if not path.is_file()]
    if missing:
        raise SystemExit(f"Input file not found: {missing[0]}")

    output = Path(args.output).expanduser().resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    database_root = Path(args.database_root).expanduser().resolve()
    if not database_root.is_dir():
        raise SystemExit(f"Database root not found: {database_root}")

    launcher_commit, launcher_dirty = git_identity(root)
    command = [
        "docker", "run", "--rm", "--platform", "linux/amd64",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--env", f"ISOSCOUT_IMAGE_REFERENCE={args.image}",
        "--env", f"ISOSCOUT_IMAGE_ID={runtime_image_id}",
        "--env", f"ISOSCOUT_SOURCE_COMMIT={image_revision}",
        "--env", f"ISOSCOUT_LAUNCHER_COMMIT={launcher_commit}",
        "--env", f"ISOSCOUT_LAUNCHER_DIRTY={launcher_dirty}",
    ]
    mounted: set[tuple[str, str]] = set()
    for name in ("ISOSCOUT_CLI_VERSION", "ISOSCOUT_EXECUTOR"):
        if name in os.environ:
            command.extend(["--env", f"{name}={os.environ[name]}"])
    container_inputs = []
    for index, path in enumerate(inputs, start=1):
        suffix = "".join(path.suffixes)
        target = f"/inputs/input_{index}{suffix}"
        command.extend(["--volume", f"{path}:{target}:ro"])
        container_inputs.append(target)
    command.extend([
        "--volume", f"{database_root}:/databases:ro",
        "--volume", f"{output.parent}:/results",
    ])

    cache = None
    if args.comparison_cache:
        cache = Path(args.comparison_cache).expanduser().resolve()
        cache.mkdir(parents=True, exist_ok=True)
        command.extend(["--volume", f"{cache}:/comparison-cache"])

    for option in ("refinement_policy", "assembly_manifest", "assembly_comparison_policy"):
        value = getattr(args, option)
        if value:
            mounted_path = add_same_path_mount(command, Path(value), mounted)
            setattr(args, option, str(mounted_path))
            if option == "assembly_manifest":
                mount_manifest_assemblies(command, mounted_path, mounted)

    command.extend([args.image, *container_inputs, "--database-root", "/databases", "--output", f"/results/{output.name}"])
    if args.organism:
        command.extend(["--organism", args.organism])
    command.append("--screen-only" if args.screen_only else "--snp-resolve")
    if args.snp_expand:
        command.append("--snp-expand")
    for option in (
        "snp_selection_mode", "query_accession", "snp_backend", "comparison_workers",
        "refinement_policy", "assembly_manifest", "assembly_comparison_policy",
    ):
        value = getattr(args, option)
        if value is not None:
            command.extend([f"--{option.replace('_', '-')}", str(value)])
    if cache:
        command.extend(["--comparison-cache", "/comparison-cache"])
    return command


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=os.environ.get("ISOSCOUT_CLI_PROG"),
        description="Run IsoScout in its pinned Docker environment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("inputs", nargs="+", help="One assembly or an R1/R2 FASTQ pair")
    parser.add_argument("--output", required=True, help="New host output directory")
    parser.add_argument("--database-root", default=os.environ.get("ISOSCOUT_DATABASE_ROOT", "~/.isoscout/databases"))
    parser.add_argument("--image", default=DEFAULT_IMAGE)
    parser.add_argument("--organism", choices=ORGANISMS)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--screen-only", action="store_true", help="Run screening without SNP comparison.")
    mode.add_argument("--snp-resolve", action="store_true", help="Compatibility alias for the default full analysis.")
    parser.add_argument("--snp-expand", action="store_true")
    parser.add_argument("--snp-selection-mode", choices=("adaptive", "all_returned"))
    parser.add_argument("--query-accession")
    parser.add_argument("--refinement-policy")
    parser.add_argument("--snp-backend", choices=("minimap2", "mummer", "ska"))
    parser.add_argument("--assembly-manifest")
    parser.add_argument("--comparison-cache")
    parser.add_argument("--assembly-comparison-policy")
    parser.add_argument("--comparison-workers", type=int)
    parser.add_argument("--print-command", action="store_true", help="Print the resolved Docker command without running it")
    args = parser.parse_args(argv)
    if args.screen_only and any((args.snp_expand, args.snp_selection_mode, args.snp_backend,
                                args.refinement_policy, args.assembly_manifest, args.comparison_cache,
                                args.assembly_comparison_policy, args.comparison_workers is not None)):
        parser.error("--screen-only cannot be combined with SNP-comparison options")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    runtime_image_id, image_revision = image_identity(args.image)
    command = build_command(args, root, runtime_image_id, image_revision)
    if args.print_command:
        print(" ".join(shlex.quote(part) for part in command))
        return 0
    return subprocess.run(command).returncode


if __name__ == "__main__":
    raise SystemExit(main())
