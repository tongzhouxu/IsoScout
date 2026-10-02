from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import test_similarity_distribution as fixtures
from common import WorkflowError, load_json
from screen_isolate import run_snp_resolution
from select_snp_targets import select_targets


class AdaptiveSelectionTests(unittest.TestCase):
    def setUp(self):
        self.case = fixtures.SimilarityDiagnosticsTests()
        self.case.setUp()
        self.addCleanup(self.case.temporary.cleanup)
        self.policy = load_json(ROOT / "config" / "snp-resolution-policy.json")

    def select(self):
        return select_targets(self.case.source, self.policy, self.case.database, self.case.profile)

    def test_dense_band_expands_beyond_50_and_100(self):
        self.case.fixture([1.] * 50 + [.9995] * 100 + [.98] * 50)
        result = self.select()
        self.assertEqual(result["selected_genomes"], 150)
        self.assertEqual(result["boundary"]["stop_reason"], "score_gap")
        self.assertEqual(result["hard_ceiling_omitted"], 0)
        reasons = {row["accession"]: row["reasons"] for row in result["decisions"]}
        self.assertIn("boundary_near_tie_expansion", reasons["GCA_0149"])
        self.assertIn("beyond_boundary_score_gap", reasons["GCA_0150"])

    def test_clear_boundary_keeps_50(self):
        self.case.fixture([1.] * 50 + [.90] * 100)
        result = self.select()
        self.assertEqual(result["selected_genomes"], 50)
        self.assertTrue(result["desired_set_complete"])

    def test_adjacent_near_ties_extend_the_selection_boundary(self):
        self.case.fixture([1. - i * .001 for i in range(100)])
        result = self.select()
        self.assertEqual(result["selected_genomes"], 100)

    def test_hard_ceiling_records_each_omission(self):
        self.case.fixture([1.] * 250)
        self.case.profile["number"] = 300
        result = self.select()
        self.assertEqual(result["selected_genomes"], 200)
        self.assertEqual(result["desired_genomes"], 250)
        self.assertEqual(result["hard_ceiling_omitted"], 50)
        self.assertFalse(result["desired_set_complete"])
        self.assertEqual(sum(row["reasons"] == ["hard_resource_ceiling"] for row in result["decisions"]), 50)
        self.assertTrue(result["warnings"])

    def test_complete_returned_band_still_warns_at_retrieval_limit(self):
        self.case.fixture([1.] * 200)
        result = self.select()
        self.assertTrue(result["desired_set_complete"])
        self.assertTrue(result["retrieval_limit_reached"])
        self.assertTrue(any("unobserved" in warning for warning in result["warnings"]))

    def test_alternative_cluster_survives_dense_dominant_cluster(self):
        self.case.fixture([1.] * 240 + [.999] * 10, ["C1"] * 240 + ["C2"] * 10)
        result = self.select()
        selected = [row["cluster"] for row in result["targets"]]
        self.assertEqual(selected.count("C2"), 10)
        self.assertEqual(len(selected), 200)
        self.assertEqual(result["missing_plausible_clusters"], [])

    def test_impossible_cluster_coverage_is_explicit(self):
        self.case.fixture([1., 1., 1.], ["C1", "C2", "C3"])
        self.policy.update(initial_target_genomes=1, max_total_genomes=2)
        result = self.select()
        self.assertEqual(result["selected_genomes"], 2)
        self.assertEqual([row["accession"] for row in result["targets"]], ["GCA_0000", "GCA_0001"])
        self.assertEqual(result["missing_plausible_clusters"], ["C3"])
        self.assertFalse(result["desired_set_complete"])

    def test_duplicate_rows_do_not_consume_slots_or_collapse_distinct_accessions(self):
        self.case.fixture([1.] * 60)
        path = self.case.source / "sample_representative_matches.csv"
        with path.open("a") as handle:
            handle.write("GCA_0000,C1,1.0\n")
        result = self.select()
        self.assertEqual(result["selected_genomes"], 60)
        self.assertEqual(len(result["decisions"]), 60)
        self.assertEqual(result["duplicate_rows_removed"], 1)
        self.assertEqual(result["duplicate_decisions"][0]["reason"], "duplicate_accession")

    def test_conflicting_duplicate_is_rejected(self):
        self.case.fixture([1.])
        path = self.case.source / "sample_representative_matches.csv"
        for cluster, score in (("C1", 0.5), ("C2", 1.0)):
            with self.subTest(cluster=cluster, score=score):
                path.write_text(f"asm_acc,PDS_acc,similarity_score\nGCA_0000,C1,1.0\nGCA_0000,{cluster},{score}\n")
                with self.assertRaises(WorkflowError):
                    self.select()

    def test_selection_independent_of_input_row_order(self):
        self.case.fixture([1.] * 240, ["C2"] * 120 + ["C1"] * 120)
        original = self.select()
        for path in self.case.source.glob("*.csv"):
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                fields, rows = reader.fieldnames, list(reader)
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(reversed(rows))
        reversed_result = self.select()
        self.assertEqual([row["accession"] for row in original["targets"]], [row["accession"] for row in reversed_result["targets"]])
        self.assertEqual([(row["accession"], row["reasons"]) for row in original["decisions"]], [(row["accession"], row["reasons"]) for row in reversed_result["decisions"]])

    def test_skip_records_reasons_and_invalid_policy_fails(self):
        self.case.fixture([.1, .09])
        result = self.select()
        self.assertEqual(result["status"], "SKIPPED")
        self.assertTrue(all(row["reasons"] == ["no_candidate_above_screening_gate"] for row in result["decisions"]))
        for key, value in (("max_total_genomes", 0), ("initial_target_genomes", 201), ("policy_version", "1.0.0")):
            policy = {**self.policy, key: value}
            with self.assertRaises(WorkflowError):
                select_targets(self.case.source, policy, self.case.database, self.case.profile)

    def test_recorded_context_and_explicit_context_agree(self):
        self.case.fixture([1.] * 60)
        run = {"database": self.case.database, "command": ["mashpit", "query", "q.fa", "db", "--number", "200", "--threshold", ".85", "--tie-tolerance-hashes", "2"]}
        (self.case.source / "mashpit_run.json").write_text(json.dumps(run))
        self.assertEqual(self.select(), select_targets(self.case.source, self.policy))

    def test_actual_refinement_uses_preview_targets_and_retains_download_losses(self):
        self.case.profile = load_json(ROOT / "config" / "workflow.json")["mashpit"]
        self.case.fixture([1.] * 240 + [.999] * 10, ["C1"] * 240 + ["C2"] * 10)
        self.case.policy = self.policy
        diagnostics = self.case.analyze()
        preview = diagnostics["selection_preview"]["audit"]
        accessions = [row["accession"] for row in preview["targets"]]
        self.assertEqual(set(accessions), {row["accession"] for row in diagnostics["ranked"] if row["selected_by_current_policy"]})
        self.assertEqual([row["selected_genomes"] for row in diagnostics["clusters"]], [190, 10])
        fetch = {"commands": [], "verified": {acc: f"/tmp/{acc}.fa" for acc in accessions[:-1]}, "unavailable": accessions[-1:]}
        interpretation = {"status": "COMPARED", "warnings": []}
        with patch("screen_isolate.fetch_genomes", return_value=fetch) as download, patch("screen_isolate.run_ska", return_value={"commands": [], "distances": []}) as ska, patch("screen_isolate.interpret_snp_resolution", return_value=interpretation), patch("screen_isolate.render_from_interpretation", return_value={"status": "SKIPPED"}):
            result = run_snp_resolution(self.case.source, self.case.root / "query.fa", "C1", self.case.root / "snp", self.case.database)
        self.assertEqual(download.call_args.args[0], accessions)
        self.assertEqual(set(ska.call_args.args[0]), {"QUERY", *accessions[:-1]})
        saved = load_json(self.case.root / "snp" / "targets.json")
        self.assertEqual(saved, preview)
        self.assertEqual(result["snp_resolution"]["unavailable"], accessions[-1:])
        self.assertTrue(any("could not be downloaded" in warning for warning in result["snp_resolution"]["warnings"]))


if __name__ == "__main__":
    unittest.main()
