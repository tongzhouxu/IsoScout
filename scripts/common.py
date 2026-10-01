#!/usr/bin/env python3
"""Shared, dependency-free helpers for the fixed Mashpit workflow."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SKILL_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = SKILL_ROOT / "config"


class WorkflowError(RuntimeError):
    """A conservative, user-reportable workflow failure."""


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise WorkflowError(f"Expected a JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def executable_version(executable: str, args: Iterable[str] = ("--version",)) -> dict[str, Any]:
    resolved = shutil.which(executable)
    if not resolved:
        return {"available": False, "version": None, "path": None}
    try:
        completed = subprocess.run(
            [resolved, *args], capture_output=True, text=True, timeout=30, check=False
        )
        text = (completed.stdout or completed.stderr).strip().splitlines()
        result = {
            "available": True,
            "version": text[0] if text else "unknown",
            "path": resolved,
            "returncode": completed.returncode,
        }
        # A nonzero exit means the tool itself failed rather than printing a
        # version; the first output line is then a traceback/error fragment,
        # not a version string, so don't record it as one (observed for real:
        # sourmash 4.6.1's own `--version` hits an unrelated upstream argparse
        # bug - duplicate `sbt_combine` subparser registration - even though
        # mashpit's actual use of sourmash as a library is unaffected).
        if completed.returncode != 0:
            result["version"] = None
            result["error"] = "\n".join(text[:3]) if text else "no output"
        return result
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"available": True, "version": None, "path": resolved, "error": str(error)}


def require_executable(executable: str) -> str:
    resolved = shutil.which(executable)
    if not resolved:
        raise WorkflowError(
            f"Required executable '{executable}' was not found. Run this skill in its pinned container."
        )
    return resolved


def run_logged(command: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> int:
    cwd.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, check=False)
    return completed.returncode


def effective_database_root(cli_value: str | None) -> Path:
    if cli_value:
        return Path(cli_value).expanduser().resolve()
    environment_value = os.environ.get("ISOSCOUT_DATABASE_ROOT")
    if environment_value:
        return Path(environment_value).expanduser().resolve()
    return Path.home() / ".isoscout" / "databases"
