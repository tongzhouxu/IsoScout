from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from launch_isoscout import build_command, parse_args


class LauncherTests(unittest.TestCase):
    @patch("launch_isoscout.git_identity", return_value=("abc123", "false"))
    def test_builds_container_command_with_host_provenance(self, _identity):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = root / "sample.fna"
            sample.write_text(">sample\nACGT\n", encoding="ascii")
            databases = root / "databases"
            databases.mkdir()
            candidate = root / "candidate.fna"
            candidate.write_text(">candidate\nACGT\n", encoding="ascii")
            manifest = root / "assemblies.tsv"
            manifest.write_text("accession\tpath\tsha256\nGCA_1\tcandidate.fna\t" + "0" * 64 + "\n")
            args = parse_args([
                str(sample), "--output", str(root / "result"),
                "--database-root", str(databases), "--organism", "salmonella",
                "--assembly-manifest", str(manifest),
            ])
            command = build_command(args, ROOT, "sha256:image", "image-commit")
        image_index = command.index("isoscout:local")
        self.assertTrue(all(option != "--volume" for option in command[image_index + 1:]))
        self.assertIn("--snp-resolve", command[image_index + 1:])
        self.assertIn("--user", command[:image_index])
        self.assertIn("ISOSCOUT_IMAGE_ID=sha256:image", command)
        self.assertIn("ISOSCOUT_SOURCE_COMMIT=image-commit", command)
        self.assertIn(str(manifest), command)
        self.assertEqual(command[command.index("--organism") + 1], "salmonella")


if __name__ == "__main__":
    unittest.main()
