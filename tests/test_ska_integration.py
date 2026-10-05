"""Opt-in offline integration against the pinned SKA executable, using random DNA."""
from __future__ import annotations

import os
import random
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from interpret_snp_resolution import interpret
from run_ska import run_ska


@unittest.skipUnless(os.environ.get("ISOSCOUT_TOOL_TESTS") == "1", "Set ISOSCOUT_TOOL_TESTS=1 in the pinned container")
class SkaIntegrationTests(unittest.TestCase):
    def test_assembly_and_reads_recover_known_differences(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rng = random.Random(1041)
            sequence = "".join(rng.choices("ACGT", k=150000))
            query = root / "query.fa"
            query.write_text(">query\n" + sequence + "\n")
            references = []
            genomes = {"QUERY": str(query)}
            for sample, positions in (("A", [10000]), ("B", [10000, 30000, 50000])):
                bases = list(sequence)
                for pos in positions:
                    bases[pos] = "A" if bases[pos] != "A" else "C"
                path = root / f"{sample}.fa"
                path.write_text(f">{sample}\n" + "".join(bases) + "\n")
                genomes[sample] = str(path)
                references.append({"accession": sample, "cluster": "C1", "path": str(path)})
            r1, r2 = root / "query_R1.fastq", root / "query_R2.fastq"
            reverse = str.maketrans("ACGT", "TGCA")
            with r1.open("w") as first, r2.open("w") as second:
                for i, start in enumerate(range(0, len(sequence) - 350, 30)):
                    first.write(f"@read{i}/1\n{sequence[start:start+150]}\n+\n{'I'*150}\n")
                    second.write(f"@read{i}/2\n{sequence[start+200:start+350].translate(reverse)[::-1]}\n+\n{'I'*150}\n")
            for mode, value in (("assembly", str(query)), ("reads", [str(r1), str(r2)])):
                distances = run_ska({**genomes, "QUERY": value}, root / mode, 31)
                result = interpret(distances["distances"], references, None)
                self.assertEqual(result["nearest_samples"], ["A"])
                self.assertEqual({row["sample"]: row["snp_distance"] for row in result["ranked"]}, {"A": 1., "B": 3.})
                self.assertFalse(result["excluded_comparisons"])
if __name__ == "__main__":
    unittest.main()
