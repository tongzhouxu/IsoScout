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
from generate_report import generate_report
from package_database_release import package_database
from run_mashpit import validate_database
from interpret_snp_resolution import interpret
from refinement_expansion import assess_progress, member_coverage, schedule_focused_batch, validate_expansion_policy
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
        self.assertEqual(result["excluded_comparisons"][0]["sample"], "A")


class FocusedSchedulerTests(unittest.TestCase):
    def setUp(self):
        self.policy = load_json(ROOT / "config" / "refinement-policy.json")["expansion"]

    @staticmethod
    def interpretation(nearest=("C1",), qualifying=None, excluded=None):
        qualifying = qualifying or [("REP", nearest[0])]
        return {"ranked": [{"sample": sample, "cluster": cluster} for sample, cluster in qualifying],
                "excluded_comparisons": [{"sample": sample, "cluster": cluster} for sample, cluster in (excluded or [])],
                "nearest_clusters": list(nearest), "nearest_samples": [qualifying[0][0]], "nearest_snp_distance": 5}

    def test_focused_batch_gives_leader_meaningful_share_and_reserves_alternatives(self):
        clusters = ["C1"] + [f"C{i:02d}" for i in range(2, 71)]
        members = [{"accession": f"L{i:03d}", "cluster": "C1"} for i in range(100)]
        members += [{"accession": f"A{i:03d}", "cluster": cluster} for i, cluster in enumerate(clusters[1:])]
        batch, state, audit = schedule_focused_batch(members, set(), clusters,
            self.interpretation(), self.policy, 100, 7)
        self.assertEqual(len(batch), 25)
        self.assertGreaterEqual(sum(row["cluster"] == "C1" for row in batch), 19)
        self.assertGreaterEqual(sum(row["cluster"] != "C1" for row in batch), 2)
        self.assertIn("alternative_after", state)
        self.assertEqual(len({row["accession"] for row in batch}), 25)
        self.assertEqual(audit["leading_clusters"], ["C1"])

    def test_tied_leaders_both_receive_focus_and_alternative_remains(self):
        members = ([{"accession": f"A{i:03d}", "cluster": "C1"} for i in range(20)]
                   + [{"accession": f"B{i:03d}", "cluster": "C2"} for i in range(20)]
                   + [{"accession": f"C{i:03d}", "cluster": "C3"} for i in range(10)])
        batch, _, audit = schedule_focused_batch(members, set(), ["C1", "C2", "C3"],
            self.interpretation(("C1", "C2")), self.policy, 100, 7)
        counts = {cluster: sum(row["cluster"] == cluster for row in batch) for cluster in ("C1", "C2", "C3")}
        self.assertGreaterEqual(counts["C1"], 9)
        self.assertGreaterEqual(counts["C2"], 9)
        self.assertGreaterEqual(counts["C3"], 2)
        self.assertEqual(audit["leading_clusters"], ["C1", "C2"])

    def test_alternative_cursor_reaches_later_clusters_across_rounds(self):
        clusters = ["C1"] + [f"C{i:02d}" for i in range(2, 12)]
        members = ([{"accession": f"L{i:03d}", "cluster": "C1"} for i in range(100)]
                   + [{"accession": f"A{i:03d}", "cluster": cluster} for i, cluster in enumerate(clusters[1:])])
        attempted = set()
        state = {}
        seen_alternatives = []
        qualifying = [("REP", cluster) for cluster in clusters]
        for round_index in range(5):
            batch, state, _ = schedule_focused_batch(members, attempted, clusters,
                self.interpretation(qualifying=qualifying), self.policy, 100 - len(attempted), 7 - round_index, state)
            seen_alternatives.extend(row["cluster"] for row in batch if row["cluster"] != "C1")
            attempted.update(row["accession"] for row in batch)
        self.assertIn("C11", seen_alternatives)
        self.assertEqual(len(seen_alternatives), len(set(seen_alternatives)))

    def test_isolated_exclusion_does_not_block_but_unresolved_cluster_does(self):
        previous = self.interpretation(qualifying=[("REP", "C1")])
        current = self.interpretation(qualifying=[("REP", "C1"), ("NEW", "C1"), ("ALT", "C2")], excluded=[("BAD", "C1")])
        coverage = [{"cluster": "C1", "qualifying_additional": 20, "unexamined_additional": 10},
                    {"cluster": "C2", "qualifying_additional": 0, "unexamined_additional": 0}]
        result = assess_progress(previous, current, coverage, self.policy, False, 0)
        self.assertTrue(result["stable"])
        unresolved = coverage + [{"cluster": "C3", "qualifying_additional": 0, "unexamined_additional": 5}]
        result = assess_progress(previous, current, unresolved, self.policy, False, 0)
        self.assertFalse(result["stable"])
        self.assertIn("C3", result["unresolved_alternatives_with_pending_references"])

    def test_unresolved_challenges_block_stability_and_cluster_witness_can_support_it(self):
        previous = self.interpretation(qualifying=[("REP", "C1")])
        current = self.interpretation(qualifying=[("REP", "C1"), ("NEW", "C1")])
        coverage = [{"cluster": "C1", "qualifying_additional": 20, "unexamined_additional": 0}]
        for item in (previous, current):
            item.update(ranking_basis="candidate_anchor_shared_regions", cluster_status="AMBIGUOUS", genome_status="AMBIGUOUS")
        current["coverage_blockers"] = [{"sample": "NEW", "cluster": "C1"}]
        result = assess_progress(previous, current, coverage, self.policy, False, 1)
        self.assertFalse(result["stable"])
        self.assertEqual(result["reason"], "unresolved_candidate_challenges")
        current["cluster_status"] = "RESOLVED"
        self.assertEqual(assess_progress(previous,current,coverage,self.policy,False,1)["reason"], "nearest_result_changed")
        previous["cluster_status"] = "RESOLVED"
        result = assess_progress(previous,current,coverage,self.policy,False,0)
        self.assertTrue(result["stable"])
        self.assertEqual(result["stability_scope"],"cluster_support_only")

    def test_determinism_dedup_self_exclusion_and_strict_budget(self):
        members = ([{"accession": f"A{i}", "cluster": "C1"} for i in range(5)]
                   + [{"accession": "A1", "cluster": "C1"}, {"accession": "SELF", "cluster": "C2"},
                      {"accession": "OTHER", "cluster": "C2"}])
        args = (members, {"A0"}, ["C1", "C2"], self.interpretation(), self.policy, 3, 2)
        first, _, _ = schedule_focused_batch(*args, query_accession="SELF")
        second, _, _ = schedule_focused_batch(*args, query_accession="SELF")
        self.assertEqual(first, second)
        self.assertEqual(len(first), 3)
        self.assertEqual(len({row["accession"] for row in first}), 3)
        self.assertFalse({"A0", "SELF"} & {row["accession"] for row in first})
        with self.assertRaises(WorkflowError):
            schedule_focused_batch(members + [{"accession": "A1", "cluster": "C2"}], {"A0"},
                                   ["C1", "C2"], self.interpretation(), self.policy, 3, 2)

    def test_invalid_allocation_is_rejected(self):
        for change in ({"leading_share": 1}, {"alternative_min_slots": 25}, {"min_focused_qualifying_per_leader": 0}):
            with self.assertRaises(WorkflowError):
                validate_expansion_policy({**self.policy, **change})


class IterativeWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.case = fixtures.SimilarityDiagnosticsTests()
        self.case.setUp()
        self.addCleanup(self.case.temporary.cleanup)
        self.case.fixture([1.])
        self.membership = {"status": "AVAILABLE", "release": "fixture", "sources": [],
                           "members": [{"accession": f"MEMBER_{i:03d}", "cluster": "C1"} for i in range(100)]}

    def execute(self, membership=None, settings=None, distance_map=None, bad=None, query_accession=None, ska_fail_round=None, comparability=None):
        membership = membership or self.membership
        distance_map = distance_map or {}
        bad = bad or set()
        def fetch(accessions, *_):
            return {"commands": [], "verified": {acc: "/tmp/" + acc + ".fa" for acc in accessions}, "unavailable": []}
        def interpret_targets(_distances, targets, _best, _policy):
            self.assertEqual(_policy, refinement["comparability"])
            ranked = [{"sample": row["accession"], "cluster": row["cluster"],
                       "snp_distance": distance_map.get(row["accession"], 5)}
                      for row in targets if row["accession"] not in bad]
            ranked.sort(key=lambda row: (row["snp_distance"], row["sample"]))
            minimum = ranked[0]["snp_distance"] if ranked else None
            nearest = [row for row in ranked if row["snp_distance"] == minimum]
            return {"status": "COMPARED" if ranked else "INSUFFICIENT_DATA", "ranked": ranked,
                    "nearest_snp_distance": minimum, "nearest_samples": [row["sample"] for row in nearest],
                    "nearest_clusters": sorted({row["cluster"] for row in nearest}),
                    "excluded_comparisons": [{"sample": row["accession"], "cluster": row["cluster"],
                                               "reasons": ["insufficient_shared_split_kmers"]} for row in targets if row["accession"] in bad],
                    "warnings": []}
        refinement = load_json(ROOT / "config" / "refinement-policy.json")
        refinement["expansion"].update(settings or {})
        refinement["comparability"].update(comparability or {})
        real_load = load_json
        def config(path):
            return refinement if path.name == "refinement-policy.json" else real_load(path)
        ska_calls = [0]
        def ska(*_):
            ska_calls[0] += 1
            if ska_fail_round is not None and ska_calls[0] == ska_fail_round:
                raise WorkflowError("round tool failure")
            return {"commands": [], "distances": []}
        with patch("screen_isolate.load_json", side_effect=config), patch("screen_isolate.fetch_genomes", side_effect=fetch) as download, patch("screen_isolate.run_ska", side_effect=ska), patch("screen_isolate.interpret_snp_resolution", side_effect=interpret_targets), patch("screen_isolate.discover_members", return_value=membership), patch("screen_isolate.render_from_interpretation", return_value={"status": "SKIPPED"}):
            result = run_snp_resolution(self.case.source, self.case.root / "query.fa", "C1", self.case.root / "snp", self.case.database, snp_backend="ska",
                                        expand=True, query_accession=query_accession)
        return result["snp_resolution"], [call.args[0] for call in download.call_args_list]

    def test_closer_member_appears_after_focused_exploration(self):
        clusters = [f"C{i:02d}" for i in range(30)]
        self.case.fixture([1. - i * .001 for i in range(30)], clusters)
        membership = {"status": "AVAILABLE", "members": (
            [{"accession": f"LEAD_{i:03d}", "cluster": "C00"} for i in range(60)]
            + [{"accession": f"ALT_{i:03d}", "cluster": cluster} for i, cluster in enumerate(clusters[1:])])}
        distances = {f"GCA_{i:04d}": (5 if i == 0 else 20) for i in range(30)}
        distances.update({f"ALT_{i:03d}": 20 for i in range(29)})
        distances["LEAD_030"] = 1
        result, downloads = self.execute(membership=membership, distance_map=distances, settings={"max_rounds": 4})
        self.assertIn("LEAD_030", {acc for batch in downloads for acc in batch})
        self.assertEqual(result["nearest_snp_distance"], 1)
        self.assertGreaterEqual(sum(acc.startswith("LEAD_") for batch in downloads[1:] for acc in batch), 31)
        self.assertTrue(any(row.get("allocation", {}).get("leading_clusters") == ["C00"] for row in result["expansion"]["rounds"]))
        report = generate_report({"sample": "fixture", "status": "COMPLETED", "input_type": "assembly",
                                  "mashpit_result": {"best_candidate": {"cluster": "C00", "score": 1.0},
                                                     "screening_result": "CANDIDATE"},
                                  "snp_resolution": result})
        self.assertIn("SNP comparison coverage by returned cluster", report)
        self.assertIn("| Cluster | Attempted | Downloaded | Qualifying | Excluded |", report)
        self.assertIn("Expansion decisions after each comparison", report)
        self.assertIn("All equally nearest qualifying references: LEAD_030", report)

    def test_unchanged_early_rounds_cannot_stop_before_focused_minimum(self):
        membership = {"status": "AVAILABLE", "members": [{"accession": f"MEMBER_{i:03d}", "cluster": "C1"} for i in range(200)]}
        result, downloads = self.execute(membership=membership,
            settings={"min_focused_qualifying_per_leader": 60, "max_additional_reference_attempts": 100})
        self.assertGreaterEqual(len(downloads), 4)
        self.assertNotEqual(result["expansion"]["stop_reason"], "stable_sampled_neighborhood")
        self.assertEqual(result["expansion"]["additional_reference_attempts"], 100)
        self.assertTrue(result["expansion"]["remaining_uncertainty"]["budget_or_round_limit_prevented_required_exploration"])

    def test_alternative_exploration_corrects_misleading_initial_leader(self):
        self.case.fixture([1., .995], ["C1", "C2"])
        membership = {"status": "AVAILABLE", "members": ([{"accession": f"LEAD_{i:03d}", "cluster": "C1"} for i in range(50)]
                                                   + [{"accession": "ALT_CLOSE", "cluster": "C2"}])}
        distances = {"GCA_0000": 5, "GCA_0001": 6, "ALT_CLOSE": 1}
        result, downloads = self.execute(membership=membership, distance_map=distances,
                                         settings={"max_rounds": 3})
        self.assertIn("ALT_CLOSE", downloads[1])
        self.assertEqual(result["nearest_clusters"], ["C2"])
        self.assertEqual(result["nearest_snp_distance"], 1)

    def test_finite_leader_pool_completes_when_budget_allows(self):
        self.case.fixture([1., .995], ["C1", "C2"])
        membership = {"status": "AVAILABLE", "members": ([{"accession": f"LEAD_{i:03d}", "cluster": "C1"} for i in range(15)]
                                                   + [{"accession": f"ALT_{i:03d}", "cluster": "C2"} for i in range(3)])}
        distances = {"GCA_0000": 5, "GCA_0001": 6}
        result, downloads = self.execute(membership=membership, distance_map=distances,
                                         settings={"max_rounds": 3})
        self.assertTrue({f"LEAD_{i:03d}" for i in range(15)}.issubset({acc for batch in downloads for acc in batch}))
        self.assertEqual(result["expansion"]["pending_available_members"], 0)

    def test_insufficient_evidence_explores_broadly_then_reports_uncertainty(self):
        membership = {"status": "AVAILABLE", "members": [{"accession": "NEW", "cluster": "C1"}]}
        result, downloads = self.execute(membership=membership, bad={"GCA_0000", "NEW"})
        self.assertEqual(downloads, [["GCA_0000"], ["NEW"]])
        self.assertEqual(result["expansion"]["stop_reason"], "insufficient_evidence_pool_exhausted")
        self.assertEqual(result["expansion"]["remaining_uncertainty"]["unresolved_selected_clusters"], ["C1"])

    def test_resource_limits_are_strict_and_failed_round_keeps_prior_result(self):
        result, downloads = self.execute(settings={"max_total_genomes": 4})
        self.assertEqual([len(batch) for batch in downloads], [1, 3])
        self.assertEqual(result["expansion"]["stop_reason"], "hard_resource_ceiling")
        failed, failed_downloads = self.execute(ska_fail_round=2)
        self.assertEqual(failed["expansion"]["stop_reason"], "round_failed")
        self.assertEqual(len(failed_downloads), 2)
        self.assertEqual(failed["nearest_snp_distance"], 5)

    def test_custom_comparability_is_used_for_snp_interpretation(self):
        result, _ = self.execute(comparability={"min_shared_split_kmers": 42})
        self.assertEqual(result["status"], "COMPARED")

    def test_query_self_excluded_from_expansion_and_pinned_conflict_stops(self):
        self.case.fixture([1., .995], ["C1", "C1"])
        membership = {"status": "AVAILABLE", "members": [
            {"accession": "GCA_0000", "cluster": "C1"}, {"accession": "NEW", "cluster": "C1"}]}
        result, downloads = self.execute(membership=membership, query_accession="GCA_0000")
        self.assertEqual(downloads, [["GCA_0001"], ["NEW"]])
        self.assertTrue(result["coverage"]["query_self_excluded_from_members"])
        conflict = {"status": "AVAILABLE", "members": [{"accession": "GCA_0001", "cluster": "C2"}]}
        second_case = fixtures.SimilarityDiagnosticsTests()
        second_case.setUp()
        self.addCleanup(second_case.temporary.cleanup)
        # The conflict is checked before any expansion allocation.
        result2, _ = self.execute(membership=conflict, query_accession="GCA_0000")
        self.assertEqual(result2["expansion"]["stop_reason"], "membership_conflict")

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
        self.assertEqual(refine.call_args.args[6], clean_reads)


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
