from __future__ import annotations

import csv
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import test_similarity_distribution as fixtures
from common import WorkflowError, load_json
from generate_report import generate_report
from screen_isolate import run_snp_resolution
from select_snp_targets import select_targets


class AdaptiveSelectionTests(unittest.TestCase):
    def setUp(self):
        self.case = fixtures.SimilarityDiagnosticsTests()
        self.case.setUp()
        self.addCleanup(self.case.temporary.cleanup)
        self.policy = load_json(ROOT / "config" / "snp-resolution-policy.json")

    def select(self, **kwargs):
        return select_targets(self.case.source, self.policy, self.case.database, self.case.profile, **kwargs)

    def test_rank_3_and_12_survive_gaps_above_old_tolerance(self):
        scores = [1., .999, .995] + [.995 - i * .0001 for i in range(1, 9)] + [.994] + [.994 - i * .0001 for i in range(1, 39)]
        clusters = ["C1"] * 50
        clusters[2], clusters[11] = "C2", "C3"
        self.case.fixture(scores, clusters)
        revised = self.select()
        selected = {row["accession"] for row in revised["targets"]}
        self.assertIn("GCA_0002", selected)
        self.assertIn("GCA_0011", selected)
        self.assertGreater(revised["decisions"][2]["gap_from_best"], .002)

    def test_global_top_fifty_guarantee_and_cluster_coverage_afterward(self):
        scores = [1.0 - i * .001 for i in range(120)]
        clusters = ["C1"] * 119 + ["C2"]
        self.case.fixture(scores, clusters)
        result = self.select()
        self.assertEqual([row["accession"] for row in result["targets"][:50]], [f"GCA_{i:04d}" for i in range(50)])
        self.assertEqual(result["selected_genomes"], 100)
        self.assertIn("GCA_0119", {row["accession"] for row in result["targets"]})
        self.assertEqual(result["decisions"][119]["reason"], "cluster_coverage_leader")

    def test_all_returned_includes_ranks_96_and_102(self):
        self.case.fixture([1.0 - i * .0003 for i in range(200)])
        result = self.select(mode="all_returned")
        self.assertEqual(result["selected_genomes"], 200)
        self.assertIn("GCA_0095", {row["accession"] for row in result["targets"]})
        self.assertIn("GCA_0101", {row["accession"] for row in result["targets"]})
        self.assertEqual(result["resource_budget"]["all_returned_required"], 200)

    def test_all_returned_rejects_partial_selection_when_budget_short(self):
        self.case.fixture([1.] * 120)
        self.policy["max_total_genomes"] = 100
        result = self.select(mode="all_returned")
        self.assertEqual(result["status"], "SKIPPED")
        self.assertEqual(result["selected_genomes"], 0)
        self.assertTrue(all(row["reason"] == "all_returned_budget_insufficient" for row in result["decisions"]))

    def test_rank_only_strategy_is_configurable_without_changing_first_fifty(self):
        self.case.fixture([1. - i * .001 for i in range(120)], ["C1"] * 119 + ["C2"])
        covered = self.select()
        self.policy["coverage_strategy"] = "global_rank_only"
        ranked = self.select()
        self.assertEqual(covered["targets"][:50], ranked["targets"][:50])
        self.assertIn("GCA_0119", {row["accession"] for row in covered["targets"]})
        self.assertNotIn("GCA_0119", {row["accession"] for row in ranked["targets"]})

    def test_duplicates_ties_query_self_and_determinism(self):
        self.case.fixture([1.] * 110)
        path = self.case.source / "sample_representative_matches.csv"
        with path.open("a") as handle:
            handle.write("GCA_0000,C1,1.0\n")
        first = self.select(query_accession="GCA_0000")
        self.assertEqual(first["duplicate_rows_removed"], 1)
        self.assertEqual(first["decisions"][0]["reason"], "query_self_excluded")
        self.assertNotIn("GCA_0000", {row["accession"] for row in first["targets"]})
        self.assertEqual([row["accession"] for row in first["targets"][:50]], [f"GCA_{i:04d}" for i in range(1, 51)])
        for csv_path in self.case.source.glob("*.csv"):
            with csv_path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                fields, rows = reader.fieldnames, list(reader)
            with csv_path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(reversed(rows))
        second = self.select(query_accession="GCA_0000")
        self.assertEqual([row["accession"] for row in first["targets"]], [row["accession"] for row in second["targets"]])

    def test_conflicting_duplicate_and_invalid_budget_fail(self):
        self.case.fixture([1.])
        path = self.case.source / "sample_representative_matches.csv"
        path.write_text("asm_acc,PDS_acc,similarity_score\nA,C1,1\nA,C2,1\n")
        with self.assertRaises(WorkflowError):
            self.select()
        self.case.fixture([1.])
        for change in ({"initial_target_genomes": 201}, {"max_returned_representatives": 0}, {"policy_version": "2.0.0"}):
            with self.assertRaises(WorkflowError):
                select_targets(self.case.source, {**self.policy, **change}, self.case.database, self.case.profile)

    def test_preview_matches_actual_download_request_and_records_losses(self):
        self.case.profile = load_json(ROOT / "config" / "workflow.json")["mashpit"]
        self.case.fixture([1.] * 150)
        self.case.policy = self.policy
        diagnostics = self.case.analyze()
        preview = diagnostics["selection_preview"]["audit"]
        accessions = [row["accession"] for row in preview["targets"]]
        self.assertEqual(set(accessions), {row["accession"] for row in diagnostics["ranked"] if row["selected_by_current_policy"]})
        fetch = {"commands": [], "verified": {acc: f"/tmp/{acc}.fa" for acc in accessions[:-1]}, "unavailable": accessions[-1:]}
        interpretation = {"status": "COMPARED", "ranked": [{"sample": accessions[0], "cluster": "C1", "snp_distance": 2.}], "excluded_comparisons": [], "warnings": []}
        with patch("screen_isolate.fetch_genomes", return_value=fetch) as download, patch("screen_isolate.run_ska", return_value={"commands": [], "distances": []}) as ska, patch("screen_isolate.interpret_snp_resolution", return_value=interpretation), patch("screen_isolate.render_from_interpretation", return_value={"status": "SKIPPED"}):
            actual = run_snp_resolution(self.case.source, self.case.root / "query.fa", "C1", self.case.root / "snp", self.case.database)
        self.assertEqual(download.call_args.args[0], accessions)
        self.assertEqual(set(ska.call_args.args[0]), {"QUERY", *accessions[:-1]})
        self.assertEqual(load_json(self.case.root / "snp" / "targets.json"), preview)
        coverage = actual["snp_resolution"]["coverage"]
        self.assertEqual(coverage["reference_downloads_unavailable"], accessions[-1:])
        self.assertEqual(len(coverage["returned_unexamined"]), 50)
        report = generate_report({"sample": "fixture", "status": "COMPLETED_WITH_WARNINGS",
                                  "mashpit_result": {"best_candidate": {"cluster": "C1", "score": 1.},
                                                     "alternative_candidates": [{"cluster": "C2", "score": .9}]},
                                  "snp_resolution": actual["snp_resolution"]})
        self.assertIn("Mashpit-displayed clusters versus SNP-tested references", report)
        self.assertIn("| C2 | 0 | 0 | 0 | 0 |", report)

    def test_all_returned_preview_matches_actual_two_hundred_attempts(self):
        self.case.profile = load_json(ROOT / "config" / "workflow.json")["mashpit"]
        self.case.fixture([1. - i * .0001 for i in range(200)])
        preview = self.select(mode="all_returned")
        wanted = [row["accession"] for row in preview["targets"]]
        fetch = {"commands": [], "verified": {acc: f"/tmp/{acc}.fa" for acc in wanted}, "unavailable": []}
        interpretation = {"status": "COMPARED", "ranked": [], "excluded_comparisons": [], "warnings": []}
        with patch("screen_isolate.fetch_genomes", return_value=fetch) as download, patch("screen_isolate.run_ska", return_value={"commands": [], "distances": []}), patch("screen_isolate.interpret_snp_resolution", return_value=interpretation), patch("screen_isolate.render_from_interpretation", return_value={"status": "SKIPPED"}):
            actual = run_snp_resolution(self.case.source, self.case.root / "query.fa", "C1", self.case.root / "snp", self.case.database, selection_mode="all_returned")
        self.assertEqual(download.call_args.args[0], wanted)
        self.assertEqual(len(wanted), 200)
        self.assertEqual(actual["snp_resolution"]["coverage"]["returned_unexamined"], [])

if __name__ == "__main__":
    unittest.main()
