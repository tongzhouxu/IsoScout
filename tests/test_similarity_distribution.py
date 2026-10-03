from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_similarity_distribution import analyze, run_diagnostics
from common import WorkflowError, load_json
from generate_report import generate_report
from screen_isolate import screen


class SimilarityDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "mashpit"
        self.source.mkdir()
        self.database = {"name": "listeria", "version": "test", "mashpit_database_settings": {"hash_number": 1000}}
        self.profile = {"number": 200, "threshold": .85, "tie_tolerance_hashes": 2}
        self.policy = {"policy_version": "3.0.0", "strategy": "global-ranked-cluster-coverage-v1", "default_mode": "adaptive", "coverage_strategy": "unrepresented_cluster_leaders_then_global_rank", "initial_target_genomes": 50, "max_returned_representatives": 100, "max_total_genomes": 100}
        self.workflow = {"mashpit": self.profile, "similarity_diagnostics": {"rank_checkpoints": [50, 100, 200]}}

    def fixture(self, scores, clusters=None):
        clusters = clusters or ["C1"] * len(scores)
        with (self.source / "sample_representative_matches.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["asm_acc", "PDS_acc", "similarity_score"])
            writer.writeheader()
            for index, (score, cluster) in enumerate(zip(scores, clusters)):
                writer.writerow({"asm_acc": f"GCA_{index:04d}", "PDS_acc": cluster, "similarity_score": score})
        with (self.source / "sample_cluster_candidates.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["PDS_acc", "best_similarity_score", "near_top"])
            writer.writeheader()
            for cluster in sorted(set(clusters)):
                score = max(score for score, group in zip(scores, clusters) if cluster == group)
                writer.writerow({"PDS_acc": cluster, "best_similarity_score": score, "near_top": max(scores) - score <= .002})

    def analyze(self):
        return analyze(self.source, self.database, self.profile, self.policy, [50, 100, 200])

    def test_dense_top_band_flags_both_caps_without_claiming_unseen_matches(self):
        self.fixture([1.] * 200)
        result = self.analyze()
        self.assertEqual(result["within_top_tolerance_count"], 200)
        self.assertEqual(result["selection_preview"]["omitted_within_top_tolerance_count"], 100)
        self.assertTrue(result["retrieval"]["may_split_near_tie"])
        self.assertIsNone(result["retrieval"]["additional_matches_exist"])
        self.assertTrue(result["rank_checkpoints"][0]["splits_near_tie"])
        self.assertIsNone(result["rank_checkpoints"][2]["splits_near_tie"])
        self.assertTrue(result["clusters"][0]["selection_splits_near_tie"])

    def test_clear_gap_after_rank_50_is_reported_not_chosen_as_a_new_cutoff(self):
        self.fixture([1.] * 50 + [.98] * 50 + [.90] * 20)
        result = self.analyze()
        self.assertFalse(result["rank_checkpoints"][0]["splits_near_tie"])
        self.assertAlmostEqual(result["rank_checkpoints"][0]["gap_to_next"], .02)
        self.assertEqual(result["within_top_tolerance_count"], 50)
        self.assertEqual(result["selection_preview"]["selected_genomes"], 100)
        self.assertFalse(result["retrieval"]["limit_reached"])

    def test_top_band_uses_gap_from_best_not_chaining_adjacent_gaps(self):
        self.fixture([1., .998, .996, .994])
        result = self.analyze()
        self.assertEqual(result["within_top_tolerance_count"], 2)
        self.assertTrue(result["retrieval"]["tail_is_near_tie"])
        self.assertFalse(result["retrieval"]["top_band_reaches_tail"])

    def test_low_scoring_tail_plateau_at_return_limit_is_also_flagged(self):
        self.fixture([1.] + [.91] * 199)
        result = self.analyze()
        self.assertFalse(result["retrieval"]["top_band_reaches_tail"])
        self.assertTrue(result["retrieval"]["may_split_near_tie"])

    def test_below_threshold_scores_do_not_imply_selection(self):
        self.fixture([.01, .01])
        result = self.analyze()
        self.assertFalse(result["selection_preview"]["eligible"])
        self.assertEqual(result["selection_preview"]["selected_genomes"], 0)
        self.assertIsNone(result["selection_preview"]["omitted_within_top_tolerance_count"])

    def test_empty_and_singleton_have_unknown_unobserved_boundaries(self):
        for scores in ([], [1.]):
            self.fixture(scores)
            result = self.analyze()
            self.assertEqual(result["returned_genomes"], len(scores))
            self.assertIsNone(result["retrieval"]["tail_adjacent_gap"])
            self.assertTrue(all(point["splits_near_tie"] is None for point in result["rank_checkpoints"]))

    def test_invalid_scores_fail_visibly(self):
        self.fixture([1.])
        path = self.source / "sample_representative_matches.csv"
        for score in ("nan", "inf", "-0.1", "1.1", "", "oops"):
            path.write_text(f"asm_acc,PDS_acc,similarity_score\nA,C1,{score}\n")
            with self.assertRaises(WorkflowError):
                self.analyze()

    def test_duplicate_accessions_are_counted_once_in_diagnostics(self):
        self.fixture([1.])
        path = self.source / "sample_representative_matches.csv"
        path.write_text("asm_acc,PDS_acc,similarity_score\nA,C1,1\nA,C1,1\n")
        result = self.analyze()
        self.assertEqual(result["returned_genomes"], 1)
        self.assertEqual(result["duplicate_rows"], 1)

    def test_no_default_hash_number_is_assumed(self):
        self.fixture([1.])
        self.database["mashpit_database_settings"] = {}
        with self.assertRaises(WorkflowError):
            self.analyze()

    def test_ranked_scores_and_gaps_are_sorted_numerically(self):
        self.fixture([.9, 1., .95])
        result = self.analyze()
        self.assertEqual([row["score"] for row in result["ranked"]], [1., .95, .9])
        self.assertAlmostEqual(result["ranked"][0]["gap_to_next"], .05)
        self.assertIsNone(result["ranked"][-1]["gap_to_next"])

    @patch("analyze_similarity_distribution.render_plot", return_value={"status": "FAIL", "error": "no plotting library"})
    def test_plot_failure_retains_json_and_report(self, _plot):
        self.fixture([1., .99])
        result = run_diagnostics(self.source, self.database, self.workflow, self.policy, self.root / "diagnostics")
        self.assertEqual(result["status"], "WARN")
        saved = load_json(self.root / "diagnostics" / "summary.json")
        self.assertEqual(saved["returned_genomes"], 2)
        report = generate_report({"sample": "test", "status": "COMPLETED_WITH_WARNINGS", "similarity_distribution": result})
        self.assertIn("policy preview", report)
        self.assertIn("numerical diagnostics are retained", report)

    @patch("analyze_similarity_distribution.render_plot", return_value={"status": "PASS", "path": "plot.png"})
    def test_screen_wires_diagnostics_into_result_report_and_provenance(self, _plot):
        self.fixture([1.] * 120)
        assembly = self.root / "query.fasta"
        assembly.write_text(">query\nACGT\n")
        fake_run = {"command": ["mashpit"], "database": self.database, "output_directory": str(self.source)}
        qc = {"status": "PASS", "metrics": {"total_length": 4}}
        with patch("screen_isolate.run_mashpit", return_value=fake_run), patch("screen_isolate.assess", return_value=qc), patch("screen_isolate.collect", return_value={}) as provenance:
            result = screen([assembly], self.root / "screen", self.root, organism="listeria")
        self.assertEqual(result["status"], "COMPLETED_WITH_WARNINGS")
        self.assertEqual(result["similarity_distribution"]["returned_genomes"], 120)
        self.assertNotIn("snp_resolution", result)
        self.assertIn("similarity_distribution", provenance.call_args.args[-1])
        self.assertTrue((self.root / "screen" / "similarity_distribution" / "summary.json").is_file())
        report = (self.root / "screen" / "report.md").read_text()
        self.assertIn("rank_similarity.png", report)
        self.assertIn("**20** returned representatives", report)

    def test_standalone_cli_uses_recorded_run_limit_not_current_config(self):
        self.fixture([1., 1.])
        run = {"database": self.database, "command": ["mashpit", "query", "q.fa", "db", "--number", "2", "--threshold", ".85", "--tie-tolerance-hashes", "1"]}
        (self.source / "mashpit_run.json").write_text(json.dumps(run))
        output = self.root / "standalone"
        completed = subprocess.run([sys.executable, str(ROOT / "scripts" / "analyze_similarity_distribution.py"), "--mashpit-output-dir", str(self.source), "--output-dir", str(output)], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = load_json(output / "summary.json")
        self.assertEqual(result["retrieval"]["requested_limit"], 2)
        self.assertEqual(result["tolerance"]["score_tolerance"], .001)
        self.assertTrue(result["retrieval"]["limit_reached"])

    def test_diagnostic_error_does_not_prevent_snp_refinement(self):
        self.fixture([1.])
        assembly = self.root / "query.fasta"
        assembly.write_text(">query\nACGT\n")
        self.database["mashpit_database_settings"] = {}
        fake_run = {"command": ["mashpit"], "database": self.database, "output_directory": str(self.source)}
        qc = {"status": "PASS", "metrics": {"total_length": 4}}
        snp = {"snp_resolution": {"status": "SKIPPED", "reason": "Test fixture"}, "commands": []}
        with patch("screen_isolate.run_mashpit", return_value=fake_run), patch("screen_isolate.assess", return_value=qc), patch("screen_isolate.collect", return_value={}), patch("screen_isolate.run_snp_resolution", return_value=snp) as refine:
            result = screen([assembly], self.root / "screen", self.root, organism="listeria", snp_resolve=True)
        self.assertEqual(result["status"], "COMPLETED_WITH_WARNINGS")
        self.assertEqual(result["similarity_distribution"]["status"], "ERROR")
        saved = load_json(self.root / "screen" / "similarity_distribution" / "summary.json")
        self.assertEqual(saved["status"], "ERROR")
        refine.assert_called_once()


if __name__ == "__main__":
    unittest.main()
