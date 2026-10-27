"""Dispatch to the existing workflow in an isolated Python process.

Running the bundled scripts as scripts preserves their imports, policy paths,
and legacy entry points. No scientific decisions are duplicated in this CLI.
"""

from __future__ import annotations

import argparse
from importlib.resources import files
import os
import sys

from . import __version__


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="isoscout",
        description="Screen an isolate with the same fixed workflow used by the IsoScout skill.",
    )
    parser.add_argument("--version", action="version", version=f"IsoScout CLI {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    screen = commands.add_parser(
        "screen", add_help=False,
        help="Screen an assembly or paired Illumina reads; use screen --help for options.",
        description="Choose an executor, then pass the workflow's screening options.",
    )
    screen.add_argument(
        "--executor", choices=("docker", "native"), default="docker",
        help="docker (default): pinned container; native: tools already installed on PATH",
    )
    screen.add_argument("-h", "--help", action="store_true", dest="screen_help")
    args, workflow_args = parser.parse_known_args(argv)
    script_name = "launch_isoscout.py" if args.executor == "docker" else "screen_isolate.py"
    script = files("isoscout._workflow").joinpath(script_name)
    if args.screen_help:
        screen.print_help()
        print("\nOptions for the selected executor:")
        sys.stdout.flush()
        workflow_args = ["--help"]
    environment = os.environ.copy()
    environment["ISOSCOUT_CLI_VERSION"] = __version__
    environment["ISOSCOUT_EXECUTOR"] = args.executor
    environment["ISOSCOUT_CLI_PROG"] = f"isoscout screen --executor {args.executor}"
    # exec preserves the workflow's exit status and normal terminal signals.
    os.execve(sys.executable, [sys.executable, str(script), *workflow_args], environment)
    return 0  # Only reachable when execve is replaced by a test double.
