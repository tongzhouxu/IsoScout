from __future__ import annotations

import json
import sqlite3
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import test_similarity_distribution as fixtures
from common import WorkflowError, load_json, sha256_file
from expand_cluster_members import parse_members, discover_members, download_file
from package_database_release import package_database
from run_mashpit import validate_database
from interpret_snp_resolution import interpret
from refinement_expansion import feedback, next_batch
from run_ska import build_file_list, parse_distance_table, DISTANCE_HEADER
from screen_isolate import run_snp_resolution, screen


def comparison(a, b, distance, matches=1000000, missing=.01):
    return {"sample1": a, "sample2": b, "snp_distance": distance, "match_count": matches, "mismatch_count": 100, "mismatch_proportion": missing}


class MembershipTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.meta, self.isolates = self.root / "meta.tsv", self.root / "isolates.tsv"
        self.meta.write_text("target_acc\tasm_acc\nT1\tGCA_000001.1\nT2\tGCA_000002.1\nT3\tNULL\nT4\tGCA_000004.1\n")
        self.isolates.write_text("target_acc\tPDS_acc\nT1\tC1\nT2\tC1\nT3\tC1\nT4\tC2\nT5\tC1\n")

    def test_membership_joins_targets_and_records_coverage_gaps(self):
        result = parse_members(self.meta, self.isolates, ["C1", "absent"])
        self.assertEqual([item["accession"] for item in result["members"]], ["GCA_000001.1", "GCA_000002.1"])
        self.assertEqual(result["targets_without_metadata"], ["T5"])
        self.assertEqual(result["targets_without_assembly"], ["T3"])
        self.assertEqual(result["clusters_without_members"], ["absent"])

    def test_conflicting_membership_fails(self):
        self.isolates.write_text("target_acc\tPDS_acc\nT1\tC1\nT1\tC2\n")
        with self.assertRaises(WorkflowError):
            parse_members(self.meta, self.isolates, ["C1", "C2"])

    def test_release_pinning_and_source_checksums(self):
        urls = []
        def fake_download(url, path, limit):
            urls.append(url)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((self.meta if "/Metadata/" in url else self.isolates).read_bytes())
            return {"url": url, "sha256": "fixture", "path": str(path)}
        with patch("expand_cluster_members.download_file", side_effect=fake_download):
            result = discover_members({"version": "PDG000000043.415", "species": "Cronobacter"}, ["C1"], self.root / "output", {"max_metadata_file_bytes": 100000})
        self.assertEqual(result["status"], "AVAILABLE")
        self.assertEqual(len(result["sources"]), 2)
        self.assertTrue(all("/PDG000000043.415/" in url and "latest" not in url for url in urls))
        self.assertTrue((self.root / "output" / "membership.json").exists())

    def test_unavailable_release_does_not_substitute_latest(self):
        with patch("expand_cluster_members.download_file", side_effect=OSError("archived release unavailable")) as download:
            result = discover_members({"version": "PDG000000043.415", "species": "Cronobacter"}, ["C1"], self.root / "output", {"max_metadata_file_bytes": 100000})
        self.assertEqual(result["status"], "UNAVAILABLE")
        self.assertEqual(download.call_count, 1)
        self.assertEqual(result["members"], [])

    def test_invalid_source_metadata_never_downloads(self):
        with patch("expand_cluster_members.download_file") as download:
            result = discover_members({"version": "latest", "species": "../other"}, ["C1"], self.root / "output", {})
        self.assertEqual(result["status"], "UNAVAILABLE")
        download.assert_not_called()

    def test_redirect_is_rejected(self):
        with patch("expand_cluster_members.urllib.request.urlopen") as open_url:
            open_url.return_value.__enter__.return_value.geturl.return_value = "https://ftp.ncbi.nlm.nih.gov/latest"
            with self.assertRaises(WorkflowError):
                download_file("https://ftp.ncbi.nlm.nih.gov/pinned", self.root / "file", 100)

    def source_database(self):
        directory = self.root / "cronobacter"
        directory.mkdir()
        db = directory / "cronobacter.db"
        with sqlite3.connect(db) as connection:
            connection.execute("CREATE TABLE DESC (name TEXT PRIMARY KEY, value TEXT)")
            connection.executemany("INSERT INTO DESC VALUES (?, ?)", [
                ("Species", "Cronobacter"), ("Hash_number", "1000"), ("Kmer_size", "31"),
            ])
            connection.execute("CREATE TABLE REPRESENTATIVE (PDT_acc TEXT, PDS_acc TEXT, asm_acc TEXT)")
            connection.execute("INSERT INTO REPRESENTATIVE VALUES ('T1', 'C1', 'GCA_000001.1')")
        (directory / "cronobacter.sig").write_text("fixture")
        (directory / "database.json").write_text(json.dumps({
            "name": "cronobacter", "version": "PDG000000043.415", "build_date": "test",
            "checksum_file": "cronobacter.db", "checksum": sha256_file(db), "source": "fixture",
        }))
        return directory

    def test_packaged_snapshot_supplies_members_offline_and_detects_corruption(self):
        source = self.source_database()
        package = self.root / "package" / "cronobacter"
        archive = self.root / "package" / "cronobacter.tar.gz"
        result = package_database(source, package, {"max_metadata_file_bytes": 100000}, self.meta, self.isolates, archive)
        self.assertEqual(result["membership_snapshot"]["verified_representatives"], 1)
        self.assertTrue(archive.is_file())
        with tarfile.open(archive, "r:gz") as contents:
            self.assertIn("cronobacter/database.json", contents.getnames())
            self.assertIn("cronobacter/membership/PDG000000043.415.metadata.tsv.gz", contents.getnames())
        database = validate_database(package, "cronobacter")
        with patch("expand_cluster_members.download_file") as download:
            found = discover_members(database, ["C1"], self.root / "analysis", {"max_metadata_file_bytes": 100000})
        download.assert_not_called()
        self.assertEqual(found["source_mode"], "bundled_database_release")
        self.assertEqual([item["accession"] for item in found["members"]], ["GCA_000001.1", "GCA_000002.1"])
        table = package / database["membership_snapshot"]["all_isolates"]["path"]
        with table.open("ab") as handle:
            handle.write(b"corrupt")
        with self.assertRaisesRegex(WorkflowError, "checksum or size"):
            validate_database(package, "cronobacter")
        with patch("expand_cluster_members.download_file") as download:
            unavailable = discover_members(database, ["C1"], self.root / "corrupt-analysis", {"max_metadata_file_bytes": 100000})
        self.assertEqual(unavailable["status"], "UNAVAILABLE")
        download.assert_not_called()

    def test_packaging_rejects_release_tables_that_disagree_with_representatives(self):
        source = self.source_database()
        self.isolates.write_text("target_acc\tPDS_acc\nT1\tC2\n")
        output = self.root / "package" / "cronobacter"
        with self.assertRaisesRegex(WorkflowError, "disagree"):
            package_database(source, output, {"max_metadata_file_bytes": 100000}, self.meta, self.isolates)
        self.assertFalse(output.exists())

    def test_snapshot_from_another_pdg_version_cannot_replace_bundled_membership(self):
        source = self.source_database()
        package = self.root / "package" / "cronobacter"
        package_database(source, package, {"max_metadata_file_bytes": 100000}, self.meta, self.isolates)
        manifest = load_json(package / "database.json")
        manifest["membership_snapshot"]["release"] = "PDG000000043.416"
        (package / "database.json").write_text(json.dumps(manifest))
        with self.assertRaisesRegex(WorkflowError, "does not match"):
            validate_database(package, "cronobacter")
        with patch("expand_cluster_members.download_file") as download:
            result = discover_members({**manifest, "database_directory": str(package)}, ["C1"], self.root / "analysis", {"max_metadata_file_bytes": 100000})
        self.assertEqual(result["status"], "UNAVAILABLE")
        download.assert_not_called()


class InterpretationTests(unittest.TestCase):
    def test_tiny_overlap_cannot_win_with_zero_snps(self):
        rows = [comparison("QUERY", "A", 0, 1, .999), comparison("QUERY", "B", 2), comparison("A", "B", 2, 1, .999)]
        result = interpret(rows, [{"accession": "A", "cluster": "C1"}, {"accession": "B", "cluster": "C2"}], "C1")
        self.assertEqual(result["nearest_sample"], "B")
        self.assertEqual(len(result["excluded_comparisons"]), 1)
        self.assertIn("insufficient_shared_split_kmers", result["excluded_comparisons"][0]["reasons"])

    def test_zero_ties_across_clusters_are_ambiguous(self):
        rows = [comparison("QUERY", "A", 0), comparison("QUERY", "B", 0), comparison("A", "B", 0)]
        result = interpret(rows, [{"accession": "A", "cluster": "C1"}, {"accession": "B", "cluster": "C2"}], "C1")
        self.assertEqual(result["status"], "AMBIGUOUS")
        self.assertEqual(result["nearest_samples"], ["A", "B"])
        self.assertIsNone(result["nearest_cluster"])
        self.assertIsNone(result["agrees_with_mash_top_candidate"])
        self.assertIn("no unique cluster", result["confidence"]["statement"])
        self.assertNotIn("exact match", result["confidence"]["statement"])

    def test_no_qualified_comparison_is_insufficient(self):
        result = interpret([comparison("QUERY", "A", 0, 1)], [{"accession": "A", "cluster": "C1"}], "C1")
        self.assertEqual(result["status"], "INSUFFICIENT_DATA")
        self.assertEqual(feedback(None, result)["action"], "STOP")

    def test_feedback_tracks_improvement_ties_and_stability(self):
        previous = {"ranked": [1], "nearest_snp_distance": 5, "nearest_samples": ["A"], "nearest_clusters": ["C1"]}
        changed = {**previous, "nearest_snp_distance": 2, "nearest_samples": ["B"]}
        self.assertIn("closer_neighbor_found", feedback(previous, changed)["reason"])
        self.assertTrue(feedback(previous, previous)["stable"])
        self.assertFalse(feedback(previous, {**previous, "nearest_snp_distance": 6})["stable"])
        self.assertFalse(feedback(previous, {**previous, "nearest_clusters": ["C1", "C2"]})["stable"])

    def test_batches_are_balanced_and_budgeted_without_repeats(self):
        members = [{"accession": a, "cluster": c} for a, c in [("A", "C1"), ("B", "C1"), ("C", "C2"), ("D", "C2")]]
        batch = next_batch(members, {"A"}, ["C1", "C2"], 25, 2)
        self.assertEqual([row["accession"] for row in batch], ["B", "C"])


class IterativeWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.case = fixtures.SimilarityDiagnosticsTests()
        self.case.setUp()
        self.addCleanup(self.case.temporary.cleanup)
        self.case.fixture([1.])
        self.membership = {"status": "AVAILABLE", "release": "fixture", "sources": [], "members": [{"accession": f"MEMBER_{i:03d}", "cluster": "C1"} for i in range(100)]}

    def execute(self, interpretations, membership=None, settings=None, ska_results=None, failed_members=False, query_accession=None):
        def fetch(accessions, *_):
            if failed_members and all(acc.startswith("MEMBER_") for acc in accessions):
                return {"commands": [], "verified": {}, "unavailable": accessions}
            return {"commands": [], "verified": {acc: "/tmp/" + acc + ".fa" for acc in accessions}, "unavailable": []}
        refinement = load_json(ROOT / "config" / "refinement-policy.json")
        refinement["expansion"].update(settings or {})
        real_load = load_json
        def config(path):
            return refinement if path.name == "refinement-policy.json" else real_load(path)
        with patch("screen_isolate.load_json", side_effect=config), patch("screen_isolate.fetch_genomes", side_effect=fetch) as download, patch("screen_isolate.run_ska", return_value={"commands": [], "distances": []}, side_effect=ska_results), patch("screen_isolate.interpret_snp_resolution", side_effect=interpretations), patch("screen_isolate.discover_members", return_value=membership or self.membership), patch("screen_isolate.render_from_interpretation", return_value={"status": "SKIPPED"}):
            result = run_snp_resolution(self.case.source, self.case.root / "query.fa", "C1", self.case.root / "snp", self.case.database, expand=True, query_accession=query_accession)
        return result["snp_resolution"], [call.args[0] for call in download.call_args_list]

    @staticmethod
    def comparison(distance=5, samples=None):
        return {"status": "COMPARED", "ranked": [{"sample": (samples or ["GCA_0000"])[0]}], "nearest_snp_distance": distance, "nearest_samples": samples or ["GCA_0000"], "nearest_clusters": ["C1"], "excluded_comparisons": [], "warnings": []}

    def test_stability_stop_is_explicitly_incomplete(self):
        result, downloads = self.execute([self.comparison() for _ in range(3)])
        self.assertEqual(result["expansion"]["stop_reason"], "stable_sampled_neighborhood")
        self.assertEqual([len(batch) for batch in downloads], [1, 25, 25])
        self.assertEqual(result["expansion"]["pending_available_members"], 50)
        self.assertEqual(len({acc for batch in downloads for acc in batch}), 51)
        self.assertTrue(any("incompletely examined" in warning for warning in result["warnings"]))

    def test_improvement_resets_stability_and_respects_round_limit(self):
        results = [self.comparison(), self.comparison(2, ["MEMBER_000"]), self.comparison(2, ["MEMBER_000"])]
        result, _ = self.execute(results, settings={"max_rounds": 3})
        self.assertEqual(result["expansion"]["stop_reason"], "round_limit")
        self.assertIn("closer_neighbor_found", result["expansion"]["rounds"][1]["feedback"]["reason"])

    def test_total_budget_is_never_exceeded(self):
        result, downloads = self.execute([self.comparison(), self.comparison()], settings={"max_total_genomes": 4})
        self.assertEqual([len(batch) for batch in downloads], [1, 3])
        self.assertEqual(result["expansion"]["stop_reason"], "hard_resource_ceiling")

    def test_separate_expansion_budget_limits_later_attempts(self):
        result, downloads = self.execute([self.comparison(), self.comparison()],
                                         settings={"max_additional_reference_attempts": 3})
        self.assertEqual([len(batch) for batch in downloads], [1, 3])
        self.assertEqual(result["expansion"]["additional_reference_attempts"], 3)
        self.assertEqual(result["expansion"]["stop_reason"], "additional_expansion_budget")

    def test_missing_snapshot_keeps_initial_comparison(self):
        result, downloads = self.execute([self.comparison()], membership={"status": "UNAVAILABLE", "error": "snapshot absent", "members": []})
        self.assertEqual(result["status"], "COMPARED")
        self.assertEqual(result["expansion"]["stop_reason"], "membership_unavailable")
        self.assertEqual(len(downloads), 1)

    def test_insufficient_overlap_stops_expansion(self):
        result, downloads = self.execute([{"status": "INSUFFICIENT_DATA", "ranked": [], "warnings": []}])
        self.assertEqual(result["expansion"]["stop_reason"], "insufficient_comparability")
        self.assertEqual(len(downloads), 1)

    def test_later_round_failure_preserves_successful_comparison(self):
        result, downloads = self.execute([self.comparison()], ska_results=[{"commands": [], "distances": []}, WorkflowError("round tool failure")])
        self.assertEqual(result["status"], "COMPARED")
        self.assertEqual(result["nearest_snp_distance"], 5)
        self.assertEqual(result["expansion"]["stop_reason"], "round_failed")
        self.assertEqual(result["expansion"]["rounds"][-1]["compared_total"], 1)
        self.assertEqual(len(downloads), 2)

    def test_first_round_failure_records_error(self):
        result, _ = self.execute([], ska_results=[WorkflowError("first round tool failure")])
        self.assertEqual(result["status"], "ERROR")
        self.assertEqual(result["error"], "first round tool failure")

    def test_membership_conflict_keeps_initial_comparison(self):
        membership = {"status": "AVAILABLE", "members": [{"accession": "GCA_0000", "cluster": "C2"}]}
        result, downloads = self.execute([self.comparison()], membership=membership)
        self.assertEqual(result["status"], "COMPARED")
        self.assertEqual(result["expansion"]["stop_reason"], "membership_conflict")
        self.assertIsNone(result["expansion"]["pending_available_members"])
        self.assertEqual(len(downloads), 1)

    def test_failed_downloads_do_not_count_as_stable_exploration(self):
        result, _ = self.execute([self.comparison() for _ in range(3)], settings={"max_rounds": 3}, failed_members=True)
        self.assertEqual(result["expansion"]["stop_reason"], "round_limit")
        self.assertEqual(result["expansion"]["rounds"][-1]["feedback"]["reason"], "no_new_comparisons")

    def test_exhausted_member_pool_has_no_unexamined_members(self):
        membership = {"status": "AVAILABLE", "members": [{"accession": "NEW", "cluster": "C1"}]}
        result, _ = self.execute([self.comparison(), self.comparison()], membership=membership)
        self.assertEqual(result["expansion"]["stop_reason"], "available_member_pool_exhausted")
        self.assertEqual(result["expansion"]["pending_available_members"], 0)

    def test_alternative_selected_cluster_members_join_expansion(self):
        self.case.fixture([1., .995], ["C1", "C2"])
        membership = {"status": "AVAILABLE", "members": [
            {"accession": "MEMBER_A", "cluster": "C1"},
            {"accession": "MEMBER_B", "cluster": "C2"},
        ]}
        nearest_c2 = {**self.comparison(), "nearest_clusters": ["C2"]}
        result, downloads = self.execute([nearest_c2, nearest_c2], membership=membership,
                                         settings={"batch_size": 2, "max_rounds": 2})
        self.assertEqual(downloads[0], ["GCA_0000", "GCA_0001"])
        self.assertEqual(downloads[1], ["MEMBER_B", "MEMBER_A"])
        self.assertEqual(result["coverage"]["members_unexamined"], [])

    def test_selected_rank_only_member_priority_is_configurable(self):
        self.case.fixture([1., .995], ["C1", "C2"])
        membership = {"status": "AVAILABLE", "members": [
            {"accession": "MEMBER_A", "cluster": "C1"},
            {"accession": "MEMBER_B", "cluster": "C2"},
        ]}
        nearest_c2 = {**self.comparison(), "nearest_clusters": ["C2"]}
        _, downloads = self.execute([nearest_c2, nearest_c2], membership=membership,
                                    settings={"batch_size": 2, "max_rounds": 2, "cluster_priority": "selected_rank_only"})
        self.assertEqual(downloads[1], ["MEMBER_A", "MEMBER_B"])

    def test_query_self_is_not_reintroduced_by_member_expansion(self):
        self.case.fixture([1., .995], ["C1", "C1"])
        membership = {"status": "AVAILABLE", "members": [
            {"accession": "GCA_0000", "cluster": "C1"},
            {"accession": "MEMBER_A", "cluster": "C1"},
        ]}
        result, downloads = self.execute([self.comparison(), self.comparison()], membership=membership,
                                         settings={"batch_size": 2, "max_rounds": 2}, query_accession="GCA_0000")
        self.assertEqual(downloads, [["GCA_0001"], ["MEMBER_A"]])
        self.assertTrue(result["coverage"]["query_self_excluded_from_members"])

    def test_read_screen_passes_cleaned_reads_to_refinement(self):
        root = self.case.root
        r1, r2 = root / "sample_R1.fastq", root / "sample_R2.fastq"
        for path in (r1, r2):
            path.write_text("@sample\nACGT\n+\nIIII\n")
        assembly = root / "assembly.fa"
        assembly.write_text(">query\nACGT\n")
        clean_reads = [str(root / "clean_R1.fastq.gz"), str(root / "clean_R2.fastq.gz")]
        built = {"commands": [], "assembly_path": str(assembly), "read_paths": clean_reads, "read_qc": {"status": "PASS", "total_bases": 400}}
        fake_run = {"command": ["mashpit"], "database": self.case.database, "output_directory": str(self.case.source)}
        with patch("screen_isolate.run_workflow", return_value=built), patch("screen_isolate.assess", return_value={"status": "PASS", "metrics": {"total_length": 4}}), patch("screen_isolate.run_mashpit", return_value=fake_run), patch("screen_isolate.run_snp_resolution", return_value={"snp_resolution": {"status": "SKIPPED"}, "commands": []}) as refine, patch("screen_isolate.collect", return_value={}), patch("analyze_similarity_distribution.render_plot", return_value={"status": "SKIPPED"}):
            screen([r1, r2], root / "screen", root, organism="listeria", snp_resolve=True)
        self.assertEqual(refine.call_args.args[-3], clean_reads)


class ReadAndDistanceContractTests(unittest.TestCase):
    def test_paired_read_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.tsv"
            build_file_list({"QUERY": ["r1.fastq.gz", "r2.fastq.gz"], "REF": "ref.fa"}, path)
            self.assertIn("QUERY\tr1.fastq.gz\tr2.fastq.gz", path.read_text())
            with self.assertRaises(WorkflowError):
                build_file_list({"QUERY": "reads.fastq"}, path)
            with self.assertRaises(WorkflowError):
                build_file_list({"QUERY": ["reads.fastq"]}, path)

    def test_malformed_and_nonfinite_distances_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "dist.tsv"
            for row in ("A\tB\t0", "A\tB\tnan\t0\t100\t0", "A\tB\t0\t2\t100\t0"):
                path.write_text("\t".join(DISTANCE_HEADER) + "\n" + row + "\n")
                with self.assertRaises(WorkflowError):
                    parse_distance_table(path)


if __name__ == "__main__":
    unittest.main()
