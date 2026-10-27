"""Wheel smoke tests: run with the Python environment containing the installed CLI.

These tests intentionally execute from temporary directories, without source on
PYTHONPATH. Docker is simulated; no biological tools, network, or genomes needed.
"""

from __future__ import annotations

from importlib import metadata, resources
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import isoscout

ROOT = Path(__file__).resolve().parents[1]
CLI = Path(sys.executable).parent / "isoscout"


class InstalledCliTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="isoscout-installed-")
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.bin = self.work / "bin"
        self.bin.mkdir()
        self.env = os.environ.copy()
        for key in list(self.env):
            if key.startswith("ISOSCOUT_") or key in {"PYTHONPATH", "PYTHONHOME"}:
                self.env.pop(key)
        self.env["PATH"] = str(self.bin)
        self.env["PYTHONPYCACHEPREFIX"] = str(self.work / "pycache")

    def run_cli(self, *args):
        return subprocess.run(
            [str(CLI), *args], cwd=self.work, env=self.env,
            capture_output=True, text=True, timeout=30,
        )

    def test_installed_package_contents_and_startup(self):
        self.assertNotIn(ROOT, Path(isoscout.__file__).resolve().parents)
        self.assertEqual(metadata.version("isoscout"), isoscout.__version__)
        policies = resources.files("isoscout.config")
        expected = sorted((ROOT / "config").glob("*.json"))
        self.assertTrue(expected)
        for source in expected:
            self.assertEqual(source.read_bytes(), policies.joinpath(source.name).read_bytes())
        result = self.run_cli("screen", "--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.run_cli("--version")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(isoscout.__version__, result.stdout)

    def test_native_qc_stop_matches_source_and_retains_provenance(self):
        # One QC stop exercises installed imports, policies and report generation.
        # Input-validation edge cases belong to the existing workflow unit tests.
        sample = self.work / "short.fna"
        sample.write_text(">contig\nACGT\n")
        installed = self.work / "installed"
        legacy = self.work / "source"
        result = self.run_cli(
            "screen", "--executor", "native", str(sample),
            "--organism", "salmonella", "--output", str(installed),
        )
        source = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "screen_isolate.py"),
             str(sample), "--organism", "salmonella", "--output", str(legacy)],
            cwd=self.work, env=self.env, capture_output=True, text=True, timeout=30,
        )
        self.assertEqual((result.returncode, source.returncode), (2, 2), result.stderr + source.stderr)
        self.assertEqual(json.loads((installed / "result.json").read_text()),
                         json.loads((legacy / "result.json").read_text()))
        self.assertEqual((installed / "report.md").read_bytes(), (legacy / "report.md").read_bytes())
        provenance = json.loads((installed / "provenance.json").read_text())
        self.assertEqual(provenance["interface"], {"cli_version": isoscout.__version__, "executor": "native"})

    def test_default_docker_forwards_arguments_provenance_and_exit_status(self):
        capture = self.work / "docker.json"
        docker = self.bin / "docker"
        docker.write_text(
            f"#!{sys.executable}\nimport json, pathlib, sys\n"
            "if sys.argv[1:3] == ['image', 'inspect']:\n"
            "    print('sha256:test-image\\timage-source')\n"
            "    sys.exit(0)\n"
            f"pathlib.Path({str(capture)!r}).write_text(json.dumps(sys.argv[1:]))\n"
            "sys.exit(17)\n"
        )
        docker.chmod(0o755)
        sample = self.work / "sample with spaces.fna"
        sample.write_text(">sample\nACGT\n")
        databases = self.work / "databases"
        databases.mkdir()
        result = self.run_cli(
            "screen", str(sample), "--output", str(self.work / "results"),
            "--database-root", str(databases), "--organism", "salmonella",
            "--image", "test:local", "--comparison-workers", "1",
        )
        self.assertEqual(result.returncode, 17, result.stderr)
        command = json.loads(capture.read_text())
        self.assertIn(f"{sample}:/inputs/input_1.fna:ro", command)
        self.assertIn(f"ISOSCOUT_CLI_VERSION={isoscout.__version__}", command)
        self.assertIn("ISOSCOUT_EXECUTOR=docker", command)
        self.assertIn("ISOSCOUT_LAUNCHER_COMMIT=unknown", command)
        self.assertIn("--snp-resolve", command)
        self.assertEqual(command[command.index("--comparison-workers") + 1], "1")


if __name__ == "__main__":
    unittest.main()
