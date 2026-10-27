from __future__ import annotations

import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from generate_report import generate_report
from rep_designations import annotate, summary_line
from screen_isolate import screen
from site_fixture import MashpitFixture


REP_URL = "https://www.cdc.gov/salmonella/rep/REPJFX01.html"


def write_database(directory: Path, rep_rows: list[tuple] | None) -> dict:
    """A minimal Mashpit .db; rep_rows=None leaves out REP_CLUSTER (databases-v2 layout)."""
    directory.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(directory / "salmonella.db") as connection:
        connection.execute("CREATE TABLE DESC (name TEXT, value TEXT)")
        connection.executemany("INSERT INTO DESC VALUES (?, ?)",
                               [("Version", "test"), ("Hash_number", "1000"), ("Kmer_size", "21")])
        if rep_rows is not None:
            connection.execute(
                "CREATE TABLE REP_CLUSTER (PDS_acc TEXT PRIMARY KEY, rep_designation TEXT, organism TEXT, "
                "anchor_pnusa TEXT, anchor_pdt TEXT, evidence_status TEXT, source_url TEXT)"
            )
            connection.executemany("INSERT INTO REP_CLUSTER VALUES (?, ?, ?, ?, ?, ?, ?)", rep_rows)
    connection.close()
    return {"name": "salmonella", "version": "test", "database_directory": str(directory),
            "mashpit_database_settings": {"hash_number": 1000}}


REP_ROWS = [("PDS000001.3", "REPJFX01", "Salmonella enterica", "PNUSAS000001", "PDT000000001.1",
             "cdc_linked", REP_URL)]


def mashpit_result(best="PDS000001.3", alternatives=("PDS000002.1",), below=False):
    return {"best_candidate": {"cluster": best, "score": .99}, "below_threshold": below,
            "alternative_candidates": [{"cluster": cluster, "score": .9} for cluster in alternatives]}


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_labels_final_and_candidate_clusters_from_rep_table(self):
        database = write_database(self.root / "db", REP_ROWS)
        rep = annotate({"mashpit_result": mashpit_result()}, database)
        self.assertEqual(rep["status"], "AVAILABLE")
        self.assertEqual(rep["final_call_basis"], "mashpit_top_candidate")
        final, = rep["final_clusters"]
        self.assertEqual((final["cluster"], final["rep_designation"], final["source_url"]),
                         ("PDS000001.3", "REPJFX01", REP_URL))
        self.assertEqual([row["status"] for row in rep["candidate_clusters"]], ["LABELLED", "NO_LABEL"])
        self.assertEqual(rep["candidate_clusters"][1]["statement"], "No REP label in this database")
        self.assertIn("does not confirm", rep["statement"])

    def test_matching_requires_exact_pds_accession(self):
        database = write_database(self.root / "db", REP_ROWS)
        for cluster in ("PDS000001.2", "PDS000001", "pds000001.3", "PDS000001.3 "):
            with self.subTest(cluster=cluster):
                rep = annotate({"mashpit_result": mashpit_result(best=cluster, alternatives=())}, database)
                self.assertEqual(rep["final_clusters"][0]["status"], "NO_LABEL")

    def test_missing_table_is_no_label_not_negative_rep_status(self):
        database = write_database(self.root / "db", None)
        rep = annotate({"mashpit_result": mashpit_result()}, database)
        self.assertEqual(rep["status"], "TABLE_ABSENT")
        self.assertFalse(rep["table_present"])
        self.assertTrue(all(row["status"] == "NO_LABEL" for row in rep["final_clusters"] + rep["candidate_clusters"]))
        text = rep["statement"] + summary_line(rep)
        self.assertIn("no rep label in this database", text.lower())
        self.assertNotIn("not a rep strain", text.lower())

    def test_snp_resolved_cluster_is_the_final_call(self):
        database = write_database(self.root / "db", REP_ROWS)
        result = {"mashpit_result": mashpit_result(best="PDS000002.1", alternatives=("PDS000001.3",)),
                  "snp_resolution": {"status": "RESOLVED", "nearest_cluster": "PDS000001.3",
                                     "nearest_clusters": ["PDS000001.3"]}}
        rep = annotate(result, database)
        self.assertEqual(rep["final_call_basis"], "snp_resolution")
        self.assertEqual(rep["final_clusters"][0]["rep_designation"], "REPJFX01")
        unresolved = {**result, "snp_resolution": {"status": "AMBIGUOUS", "nearest_cluster": None,
                                                   "nearest_clusters": ["PDS000001.3", "PDS000002.1"]}}
        rep = annotate(unresolved, database)
        self.assertEqual(rep["final_call_basis"], "snp_resolution_unresolved")
        self.assertEqual([row["cluster"] for row in rep["final_clusters"]], ["PDS000001.3", "PDS000002.1"])

    def test_below_threshold_hit_has_no_final_call(self):
        database = write_database(self.root / "db", REP_ROWS)
        rep = annotate({"mashpit_result": mashpit_result(below=True)}, database)
        self.assertIsNone(rep["final_call_basis"])
        self.assertEqual(rep["final_clusters"], [])
        self.assertEqual(rep["candidate_clusters"][0]["rep_designation"], "REPJFX01")

    def test_malformed_table_is_unreadable_not_unlabelled(self):
        directory = self.root / "db"
        database = write_database(directory, None)
        with sqlite3.connect(directory / "salmonella.db") as connection:
            connection.execute("CREATE TABLE REP_CLUSTER (PDS_acc TEXT PRIMARY KEY)")
        connection.close()
        rep = annotate({"mashpit_result": mashpit_result()}, database)
        self.assertEqual(rep["status"], "UNREADABLE")
        self.assertEqual(rep["final_clusters"][0]["status"], "UNKNOWN")

    def test_report_gives_designation_url_and_membership_caveat(self):
        database = write_database(self.root / "db", REP_ROWS)
        result = {"sample": "q", "status": "COMPLETED", "mashpit_result": mashpit_result()}
        result["rep_designations"] = annotate(result, database)
        report = generate_report(result)
        self.assertIn("## CDC REP designations", report)
        self.assertIn("| PDS000001.3 | Final call | **REPJFX01** | cdc_linked | " + REP_URL + " |", report)
        self.assertIn("| PDS000002.1 | Mashpit candidate | No REP label in this database |", report)
        self.assertIn("does not confirm that this isolate is a member of the outbreak", report)

    def test_report_without_rep_field_is_unchanged(self):
        result = {"sample": "q", "status": "COMPLETED", "mashpit_result": mashpit_result()}
        self.assertNotIn("REP", generate_report(result))


class ScreenIntegrationTests(unittest.TestCase):
    def run_screen(self, root: Path, rep_rows):
        fixture = MashpitFixture(root)
        fixture.fixture([.99, .95], ["PDS000001.3", "PDS000002.1"])
        database = write_database(root / "db", rep_rows)
        assembly = root / "query.fasta"
        assembly.write_text(">query\nACGT\n")
        fake_run = {"command": ["mashpit"], "database": database, "output_directory": str(fixture.source)}
        qc = {"status": "PASS", "metrics": {"total_length": 4}}
        with patch("screen_isolate.run_mashpit", return_value=fake_run), \
                patch("screen_isolate.assess", return_value=qc), \
                patch("screen_isolate.collect", return_value={}), \
                patch("analyze_similarity_distribution.render_plot", return_value={"status": "SKIPPED"}):
            screen([assembly], root / "screen", root, organism="salmonella", snp_resolve=False)
        return (json.loads((root / "screen" / "result.json").read_text()),
                (root / "screen" / "report.md").read_text())

    def test_screen_reports_rep_labels_and_leaves_v2_databases_unchanged(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            labelled, labelled_report = self.run_screen(Path(first), REP_ROWS)
            plain, plain_report = self.run_screen(Path(second), None)
        self.assertEqual(labelled["rep_designations"]["final_clusters"][0]["rep_designation"], "REPJFX01")
        self.assertIn("CDC REP designation of called cluster(s): PDS000001.3: REPJFX01 (" + REP_URL + ")",
                      labelled["user_summary"])
        self.assertIn("**REPJFX01**", labelled_report)
        self.assertEqual(plain["rep_designations"]["status"], "TABLE_ABSENT")
        self.assertIn("No REP label in this database", plain_report)
        self.assertNotIn("REPJFX01", plain_report)
        # Status, warnings, and the Mashpit call do not depend on the REP table.
        for key in ("status", "warnings", "mashpit_result"):
            self.assertEqual(labelled[key] if key != "mashpit_result" else labelled[key]["best_candidate"],
                             plain[key] if key != "mashpit_result" else plain[key]["best_candidate"])


if __name__ == "__main__":
    unittest.main()
