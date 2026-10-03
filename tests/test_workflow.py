from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_snp_tree import neighbor_joining
from collect_provenance import collect
from classify_with_mlst import interpret_mlst, parse_mlst_csv, route_for_organism
from common import WorkflowError, executable_version
from fetch_reference_genomes import fetch_genomes
from inspect_input import inspect
from interpret_snp_resolution import interpret as interpret_snp_resolution
from parse_mashpit_results import interpret, load_candidates, locate_tree_image
from run_assembly_workflow import run_workflow
from run_mashpit import run_mashpit, validate_database
from render_snp_tree import render_from_interpretation
from run_ska import parse_distance_table, run_ska
from select_snp_targets import select_targets
from validate_assembly import assess
from validate_fastq import validate_pair


class InputAndQcTests(unittest.TestCase):
    def test_detects_assembly_and_fastq_pair(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly = root / "sample.fna"
            assembly.write_text(">a\nACGT\n", encoding="utf-8")
            self.assertEqual(inspect([assembly])["input_type"], "assembly")
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            record1 = "@read/1\nACGT\n+\nIIII\n"
            record2 = "@read/2\nTGCA\n+\nIIII\n"
            r1.write_text(record1, encoding="ascii")
            r2.write_text(record2, encoding="ascii")
            self.assertEqual(inspect([r2, r1])["input_type"], "illumina_paired_fastq")
            self.assertEqual(validate_pair(r1, r2)["status"], "PASS")

    def test_rejects_corrupt_fastq(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            r1.write_text("@read/1\nACGT\n+\nIII\n", encoding="ascii")
            r2.write_text("@read/2\nACGT\n+\nIIII\n", encoding="ascii")
            self.assertEqual(validate_pair(r1, r2)["status"], "FAIL")

    def test_supported_organism_length_policies(self):
        cases = {
            "salmonella": 4_700_000,
            "ecoli_shigella": 5_000_000,
            "listeria": 2_900_000,
            "campylobacter": 1_700_000,
            "cronobacter": 4_500_000,
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "assembly.fna"
            for organism, length in cases.items():
                with self.subTest(organism=organism):
                    path.write_text(">contig\n" + "A" * length + "\n", encoding="ascii")
                    self.assertEqual(assess(path, organism)["status"], "PASS")

    def test_rejects_corrupt_and_poor_assembly(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            corrupt = root / "corrupt.fna"
            corrupt.write_text("not-a-header\nACGTX\n", encoding="ascii")
            self.assertEqual(assess(corrupt)["status"], "FAIL")
            poor = root / "poor.fna"
            poor.write_text("".join(f">c{i}\n" + "A" * 2000 + "\n" for i in range(600)), encoding="ascii")
            result = assess(poor)
            self.assertEqual(result["status"], "FAIL")
            self.assertTrue(any("N50" in reason for reason in result["failures"]))


class MashpitParsingTests(unittest.TestCase):
    def test_locate_tree_image_finds_the_single_match(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "sample_tree.png").write_bytes(b"fake-png")
            self.assertEqual(locate_tree_image(root), root / "sample_tree.png")

    def test_locate_tree_image_returns_none_when_absent(self):
        # Not an error: Mashpit itself skips tree generation below --threshold
        # or when fewer than two candidates qualify.
        with tempfile.TemporaryDirectory() as temporary:
            self.assertIsNone(locate_tree_image(Path(temporary)))

    def write_candidates(self, path: Path, rows: list[dict]) -> None:
        fields = ["PDS_acc", "best_similarity_score", "near_top"]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def test_clear_mashpit_match(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample_cluster_candidates.csv"
            self.write_candidates(path, [
                {"PDS_acc": "PDS0001", "best_similarity_score": "0.991", "near_top": "True"},
                {"PDS_acc": "PDS0002", "best_similarity_score": "0.940", "near_top": "False"},
            ])
            result = interpret(load_candidates(path))
            self.assertEqual(result["status"], "MATCH")
            self.assertEqual(result["best_candidate"]["cluster"], "PDS0001")
            self.assertEqual(result["screening_result"], "Top-ranked candidate")
            self.assertFalse(result["below_threshold"])

    def test_below_threshold_top_hit_is_flagged_even_when_unambiguous(self):
        # Observed for real: a cross-genus query still produces a near_top-flagged
        # "candidate" at near-zero score, since --number/--tie-tolerance-hashes group
        # the top hits regardless of absolute score; --threshold only gates local tree
        # construction. The parser must catch this independently of ambiguity.
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample_cluster_candidates.csv"
            self.write_candidates(path, [
                {"PDS_acc": "PDS0001", "best_similarity_score": "0.007", "near_top": "True"},
            ])
            result = interpret(load_candidates(path), threshold=0.85)
            self.assertFalse(result["ambiguous"])
            self.assertTrue(result["below_threshold"])
            self.assertTrue(any("below its own query threshold" in warning for warning in result["warnings"]))

    def test_ambiguous_mashpit_match(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample_cluster_candidates.csv"
            self.write_candidates(path, [
                {"PDS_acc": "PDS0001", "best_similarity_score": "0.991", "near_top": "True"},
                {"PDS_acc": "PDS0002", "best_similarity_score": "0.990", "near_top": "True"},
            ])
            result = interpret(load_candidates(path))
            self.assertEqual(result["status"], "AMBIGUOUS")
            self.assertTrue(result["ambiguous"])


class MlstRoutingTests(unittest.TestCase):
    def mlst_row(self, scheme="senterica", status="PERFECT"):
        return {
            "FILE": "sample.fna", "SCHEME": scheme, "ST": "11",
            "STATUS": status, "SCORE": "100", "ALLELES": "aroC(1)",
        }

    def test_routes_supported_pubmlst_schemes(self):
        cases = {
            "senterica": "salmonella", "ecoli": "ecoli_shigella",
            "lmonocytogenes": "listeria", "campylobacter": "campylobacter",
            "cronobacter": "cronobacter",
        }
        for scheme, expected in cases.items():
            with self.subTest(scheme=scheme):
                result = interpret_mlst(self.mlst_row(scheme))
                self.assertEqual(result["status"], "SUPPORTED")
                self.assertEqual(result["database_name"], expected)

    def test_unsupported_and_uncertain_mlst_results(self):
        self.assertEqual(interpret_mlst(self.mlst_row("kpneumoniae"))["status"], "UNSUPPORTED")
        self.assertEqual(interpret_mlst(self.mlst_row("senterica", "MIXED"))["status"], "UNCERTAIN")
        warning = interpret_mlst(self.mlst_row("senterica", "MISSING"))
        self.assertEqual(warning["status"], "SUPPORTED")
        self.assertTrue(warning["warnings"])

    def test_parse_pinned_mlst_csv_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mlst.csv"
            path.write_text(
                "FILE,SCHEME,ST,STATUS,SCORE,ALLELES\n"
                "sample.fna,ecoli,10,PERFECT,100,adk(1);fumC(2)\n",
                encoding="utf-8",
            )
            self.assertEqual(parse_mlst_csv(path)["SCHEME"], "ecoli")

    def test_user_organism_selects_database_without_classification(self):
        for organism in ("salmonella", "ecoli_shigella", "listeria", "campylobacter", "cronobacter"):
            with self.subTest(organism=organism):
                route = route_for_organism(organism)
                self.assertEqual(route["source"], "user")
                self.assertEqual(route["database_name"], organism)


class FailureAndProvenanceTests(unittest.TestCase):
    def test_executable_version_records_real_version_on_success(self):
        result = executable_version("python3", args=("--version",))
        self.assertTrue(result["available"])
        self.assertEqual(result["returncode"], 0)
        self.assertIsNotNone(result["version"])

    @patch("common.subprocess.run")
    @patch("common.shutil.which", return_value="/opt/conda/bin/sourmash")
    def test_executable_version_does_not_report_traceback_as_a_version(self, _which, run):
        # Observed for real: sourmash 4.6.1's own `--version` hits an unrelated
        # upstream argparse bug and exits nonzero with a traceback on stderr.
        # A nonzero exit must not be recorded as if the first line were a version.
        run.return_value = type("Completed", (), {
            "stdout": "", "stderr": "Traceback (most recent call last):\nArgumentError: ...",
            "returncode": 1,
        })()
        result = executable_version("sourmash")
        self.assertTrue(result["available"])
        self.assertIsNone(result["version"])
        self.assertIn("Traceback", result["error"])

    def test_missing_database(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(WorkflowError):
                validate_database(Path(temporary) / "missing", "salmonella")

    @patch("run_mashpit.validate_database")
    @patch("run_mashpit.require_executable", return_value="/usr/bin/mashpit")
    @patch("run_mashpit.run_logged", return_value=1)
    def test_mashpit_runtime_failure(self, _run, _executable, database):
        database.return_value = {"name": "salmonella", "version": "test"}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly = root / "sample.fna"
            assembly.write_text(">a\nACGT\n", encoding="ascii")
            with self.assertRaises(WorkflowError):
                run_mashpit(assembly, root / "db", "salmonella", root / "output")

    @patch("run_mashpit.validate_database")
    @patch("run_mashpit.require_executable", return_value="/usr/bin/mashpit")
    @patch("run_mashpit.run_logged", return_value=0)
    def test_fixed_mashpit_invocation(self, logged, _executable, database):
        database.return_value = {
            "name": "salmonella",
            "version": "test",
            "mashpit_database_settings": {"hash_number": 1000, "kmer_size": 31},
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly = root / "sample.fna"
            assembly.write_text(">a\nACGT\n", encoding="ascii")
            result = run_mashpit(assembly, root / "db", "salmonella", root / "output")
            command = logged.call_args.args[0]
            self.assertEqual(command[:2], ["/usr/bin/mashpit", "query"])
            self.assertEqual(command[-6:], [
                "--number", "200", "--threshold", "0.85", "--tie-tolerance-hashes", "2"
            ])
            self.assertEqual(result["status"], "PASS")

    @patch("run_assembly_workflow.validate_pair", return_value={"status": "PASS", "read_pairs": 1})
    @patch("run_assembly_workflow.require_executable", side_effect=lambda name: name)
    @patch("run_assembly_workflow.run_logged", return_value=1)
    def test_failed_assembly_workflow(self, _run, _executable, _validate):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(WorkflowError):
                run_workflow(root / "r1.fastq", root / "r2.fastq", root / "assembly")

    @patch("collect_provenance.executable_version", return_value={"available": False})
    def test_provenance_contains_checksums_and_commit(self, _version):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample.fna"
            path.write_text(">a\nACGT\n", encoding="ascii")
            value = collect([path], path, {"name": "salmonella"}, [["mashpit", "query"]])
            self.assertEqual(len(value["inputs"][0]["sha256"]), 64)
            self.assertIn("538d342", value["pinned_tool_versions"]["mashpit"])


class SnpResolutionTests(unittest.TestCase):
    def write_mashpit_output(self, output_dir: Path, cluster_rows, representative_rows) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        with (output_dir / "sample_cluster_candidates.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["PDS_acc", "best_similarity_score", "near_top"])
            writer.writeheader()
            writer.writerows(cluster_rows)
        with (output_dir / "sample_representative_matches.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["asm_acc", "biosample_acc", "PDS_acc", "similarity_score"]
            )
            writer.writeheader()
            writer.writerows(representative_rows)

    def test_select_targets_unambiguous_uses_only_top_cluster(self):
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            self.write_mashpit_output(
                output_dir,
                [
                    {"PDS_acc": "PDS0001", "best_similarity_score": "0.99", "near_top": "False"},
                    {"PDS_acc": "PDS0002", "best_similarity_score": "0.80", "near_top": "False"},
                ],
                [
                    {"asm_acc": "GCA_1", "biosample_acc": "SAMN1", "PDS_acc": "PDS0001", "similarity_score": "0.99"},
                    {"asm_acc": "GCA_2", "biosample_acc": "SAMN2", "PDS_acc": "PDS0001", "similarity_score": "0.95"},
                    {"asm_acc": "GCA_3", "biosample_acc": "SAMN3", "PDS_acc": "PDS0002", "similarity_score": "0.80"},
                ],
            )
            policy = {"policy_version": "3.0.0", "strategy": "global-ranked-cluster-coverage-v1", "default_mode": "adaptive", "coverage_strategy": "unrepresented_cluster_leaders_then_global_rank", "initial_target_genomes": 5, "max_returned_representatives": 20, "max_total_genomes": 20}
            result = select_targets(output_dir, policy, {"mashpit_database_settings": {"hash_number": 1000}}, {"number": 200, "threshold": .85, "tie_tolerance_hashes": 2})
            self.assertEqual(result["status"], "SELECTED")
            self.assertEqual(result["selected_clusters"], ["PDS0001", "PDS0002"])
            accessions = {item["accession"] for item in result["targets"]}
            self.assertEqual(accessions, {"GCA_1", "GCA_2", "GCA_3"})

    def test_select_targets_ambiguous_includes_near_top_cluster(self):
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            self.write_mashpit_output(
                output_dir,
                [
                    {"PDS_acc": "PDS0001", "best_similarity_score": "0.99", "near_top": "True"},
                    {"PDS_acc": "PDS0002", "best_similarity_score": "0.985", "near_top": "True"},
                ],
                [
                    {"asm_acc": "GCA_1", "biosample_acc": "SAMN1", "PDS_acc": "PDS0001", "similarity_score": "0.99"},
                    {"asm_acc": "GCA_2", "biosample_acc": "SAMN2", "PDS_acc": "PDS0002", "similarity_score": "0.985"},
                ],
            )
            policy = {"policy_version": "3.0.0", "strategy": "global-ranked-cluster-coverage-v1", "default_mode": "adaptive", "coverage_strategy": "unrepresented_cluster_leaders_then_global_rank", "initial_target_genomes": 5, "max_returned_representatives": 20, "max_total_genomes": 20}
            result = select_targets(output_dir, policy, {"mashpit_database_settings": {"hash_number": 1000}}, {"number": 200, "threshold": .85, "tie_tolerance_hashes": 2})
            self.assertEqual(sorted(result["selected_clusters"]), ["PDS0001", "PDS0002"])
            accessions = {item["accession"] for item in result["targets"]}
            self.assertEqual(accessions, {"GCA_1", "GCA_2"})

    @patch("fetch_reference_genomes.download_batch")
    @patch("fetch_reference_genomes.require_executable", return_value="/usr/bin/datasets")
    def test_fetch_genomes_partial_success_and_retry(self, _executable, download_batch):
        def fake_download(datasets_exe, accessions, batch_dir):
            batch_dir.mkdir(parents=True, exist_ok=True)
            if "GOOD1" in accessions:
                target = batch_dir / "ncbi_dataset" / "data" / "GOOD1"
                target.mkdir(parents=True, exist_ok=True)
                (target / "GOOD1_genomic.fna").write_text(">contig\nACGT\n", encoding="ascii")
            return (["datasets", "download"], ["datasets", "rehydrate"])

        download_batch.side_effect = fake_download
        with tempfile.TemporaryDirectory() as temporary:
            result = fetch_genomes(["GOOD1", "BAD1"], Path(temporary) / "genomes", attempts=2, retry_delay_seconds=0)
            self.assertEqual(result["status"], "PARTIAL")
            self.assertIn("GOOD1", result["verified"])
            self.assertEqual(result["unavailable"], ["BAD1"])
            self.assertEqual(download_batch.call_count, 2)

    def test_parse_distance_table_rejects_unexpected_header(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "distances.tsv"
            path.write_text("Sample1\tSample2\tSNPs\n", encoding="utf-8")
            with self.assertRaises(WorkflowError):
                parse_distance_table(path)

    def test_run_ska_requires_at_least_two_genomes(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(WorkflowError):
                run_ska({"QUERY": "query.fasta"}, Path(temporary), 31)

    @patch("run_ska.require_executable", return_value="/usr/bin/ska")
    def test_run_ska_builds_file_list_and_parses_distances(self, _executable):
        def fake_run_logged(command, cwd, stdout_path, stderr_path):
            if command[1] == "build":
                (cwd / "merged.skf").write_text("fake-skf", encoding="utf-8")
            elif command[1] == "distance":
                Path(command[3]).write_text(
                    "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count\n"
                    "QUERY\tGCA_1\t4.00\t0.02\t99000\t80\n",
                    encoding="utf-8",
                )
            return 0

        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            with patch("run_ska.run_logged", side_effect=fake_run_logged) as logged:
                result = run_ska({"QUERY": "query.fasta", "GCA_1": "gca1.fasta"}, output_dir, 31)
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["distances"], [{
                "sample1": "QUERY", "sample2": "GCA_1", "snp_distance": 4.0,
                "mismatch_proportion": 0.02, "match_count": 99000, "mismatch_count": 80,
            }])
            build_command = logged.call_args_list[0].args[0]
            self.assertEqual(build_command[:2], ["/usr/bin/ska", "build"])
            self.assertIn("-k", build_command)
            file_list = (output_dir / "ska_input.tsv").read_text(encoding="utf-8")
            self.assertIn("QUERY\tquery.fasta", file_list)
            self.assertIn("GCA_1\tgca1.fasta", file_list)

    def test_interpret_snp_resolution_flags_disagreement_with_mash(self):
        rows = [
            {"sample1": "QUERY", "sample2": "GCA_1", "snp_distance": 5.0,
             "mismatch_proportion": 0.01, "match_count": 1000000, "mismatch_count": 1},
            {"sample1": "GCA_2", "sample2": "QUERY", "snp_distance": 2.0,
             "mismatch_proportion": 0.01, "match_count": 1000000, "mismatch_count": 1},
            {"sample1": "GCA_1", "sample2": "GCA_2", "snp_distance": 6.0,
             "mismatch_proportion": 0.01, "match_count": 1000000, "mismatch_count": 1},
        ]
        targets = [
            {"accession": "GCA_1", "cluster": "PDS0001"},
            {"accession": "GCA_2", "cluster": "PDS0002"},
        ]
        result = interpret_snp_resolution(rows, targets, "PDS0001")
        self.assertEqual(result["status"], "COMPARED")
        self.assertEqual(result["nearest_sample"], "GCA_2")
        self.assertEqual(result["nearest_cluster"], "PDS0002")
        self.assertFalse(result["agrees_with_mash_top_candidate"])
        self.assertEqual(result["confidence"]["nearest_cluster"], "PDS0002")
        self.assertEqual(result["confidence"]["next_cluster"], "PDS0001")
        self.assertAlmostEqual(result["confidence"]["separation_ratio"], 2.5)
        self.assertTrue(result["newick_tree"].startswith("("))
        # Tips are labeled with their cluster so a tree viewer shows more than an
        # opaque accession; QUERY itself has no cluster, so it stays unlabeled.
        self.assertIn("GCA_1_PDS0001", result["newick_tree"])
        self.assertIn("GCA_2_PDS0002", result["newick_tree"])
        self.assertIn("QUERY", result["newick_tree"])
        self.assertNotIn("QUERY_", result["newick_tree"])
        self.assertTrue(result["warnings"])

    def test_interpret_snp_resolution_agrees_with_mash(self):
        rows = [{"sample1": "QUERY", "sample2": "GCA_1", "snp_distance": 1.0,
                  "mismatch_proportion": 0.0, "match_count": 1000000, "mismatch_count": 0}]
        targets = [{"accession": "GCA_1", "cluster": "PDS0001"}]
        result = interpret_snp_resolution(rows, targets, "PDS0001")
        self.assertTrue(result["agrees_with_mash_top_candidate"])
        self.assertEqual(result["warnings"], [])

    def test_interpret_snp_resolution_insufficient_data_without_query_rows(self):
        rows = [{"sample1": "GCA_1", "sample2": "GCA_2", "snp_distance": 1.0,
                  "mismatch_proportion": 0.0, "match_count": 1000000, "mismatch_count": 0}]
        result = interpret_snp_resolution(rows, [], None)
        self.assertEqual(result["status"], "INSUFFICIENT_DATA")

    def test_neighbor_joining_recovers_exact_additive_tree(self):
        # Distances summed along a known tree ((A:2,B:3):1,(C:4,D:5):1); NJ is
        # guaranteed to recover an additive tree exactly, so this checks the
        # implementation bit-for-bit rather than just "it runs".
        rows = [
            {"sample1": "A", "sample2": "B", "snp_distance": 5},
            {"sample1": "A", "sample2": "C", "snp_distance": 8},
            {"sample1": "A", "sample2": "D", "snp_distance": 9},
            {"sample1": "B", "sample2": "C", "snp_distance": 9},
            {"sample1": "B", "sample2": "D", "snp_distance": 10},
            {"sample1": "C", "sample2": "D", "snp_distance": 9},
        ]
        newick = neighbor_joining(rows)
        self.assertEqual(newick, "((A:2.0000,B:3.0000):1.0000,(C:4.0000,D:5.0000):1.0000);")

    def test_render_from_interpretation_skips_without_tree(self):
        result = render_from_interpretation({"ranked": []}, Path("unused.png"))
        self.assertEqual(result["status"], "SKIPPED")

    @patch("render_snp_tree.render")
    def test_render_from_interpretation_reports_failure_without_raising(self, render):
        render.side_effect = RuntimeError("matplotlib backend unavailable")
        result = render_from_interpretation(
            {"newick_tree": "(A:1,B:1);", "ranked": [{"sample": "A"}]}, Path("unused.png")
        )
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("matplotlib", result["error"])

    @patch("render_snp_tree.render")
    def test_render_from_interpretation_passes_through_on_success(self, render):
        result = render_from_interpretation(
            {"newick_tree": "(A:1,B:1);", "ranked": [{"sample": "A"}]}, Path("tree.png")
        )
        self.assertEqual(result["status"], "PASS")
        render.assert_called_once()

    def test_neighbor_joining_two_taxa(self):
        # Labels are sorted, so GCA_1 sorts before QUERY.
        rows = [{"sample1": "QUERY", "sample2": "GCA_1", "snp_distance": 10}]
        self.assertEqual(neighbor_joining(rows), "(GCA_1:5.0000,QUERY:5.0000);")


if __name__ == "__main__":
    unittest.main()
