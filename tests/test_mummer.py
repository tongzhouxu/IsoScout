"""Software contracts for cached target-only assembly comparisons."""
from __future__ import annotations
import csv
import json
import os
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from common import WorkflowError,load_json,sha256_file
from run_mummer import run_mummer,union_length
from fetch_reference_genomes import fetch_genomes,read_local_manifest
from interpret_assembly_comparison import interpret
from generate_report import generate_report
from site_fixture import comparison

class PairTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.policy=load_json(ROOT/"config/mummer-comparison-policy.json")
        self.policy["ranking"]["min_shared_bases"]=100  # Scaled synthetic fixtures.
        self.tools={"version":"4.0.0","executables":{x:{"path":x,"sha256":"fixture"} for x in ("nucmer","delta-filter","show-snps","show-coords")}}
        self.genomes={}
        for i in range(126):
            name="QUERY" if i==0 else f"A{i:03d}"
            p=self.root/(name+".fna");p.write_text(">contig\n"+"A"*(300+i)+"\n");self.genomes[name]=str(p)
        self.launches=[]
    def execute(self,cmd,directory,output,timeout):
        tool=Path(cmd[0]).name
        if tool=="nucmer":
            self.launches.append(directory)
            (directory/"pair.delta").write_text("raw alignment")
            (directory/output).write_text("")
        elif tool=="delta-filter":(directory/output).write_text("unique alignment")
        elif tool=="show-coords":
            n=len((directory/"candidate.fna").read_text().splitlines()[1])
            (directory/output).write_text(f"1\t300\t1\t300\t300\t300\t100\t300\t{n}\t100\t90\tcontig\tcontig\n")
        elif tool=="show-snps":(directory/output).write_text("")
    def run_pairs(self,g,name="out",policy=None):
        with patch("run_mummer.tool_identity",return_value=self.tools),patch("run_mummer.command_run",side_effect=self.execute):
            return run_mummer(g,self.root/name,self.root/"cache",policy)
    def test_initial_100_expand_25_and_resume_only_new_pairs(self):
        first=dict(list(self.genomes.items())[:101]);a=self.run_pairs(first)
        self.assertEqual(a["alignment_jobs"],100);self.assertEqual(len(a["comparisons"]),100)
        b=self.run_pairs(self.genomes,"expanded")
        self.assertEqual(b["alignment_jobs"],25);self.assertEqual(b["cache_hits"],100)
        c=self.run_pairs(self.genomes,"resumed")
        self.assertEqual(c["alignment_jobs"],0);self.assertEqual(c["cache_hits"],125)
        self.assertEqual(len(self.launches),125)
    def test_input_and_alignment_changes_invalidate_but_coverage_does_not(self):
        g=dict(list(self.genomes.items())[:2]);self.run_pairs(g)
        self.policy["comparability"]["min_target_aligned_fraction"]=.9
        p=self.root/"policy.json";p.write_text(json.dumps(self.policy))
        self.assertEqual(self.run_pairs(g,"coverage",p)["cache_hits"],1)
        self.policy["alignment"]["min_match"]=21;p.write_text(json.dumps(self.policy))
        self.assertEqual(self.run_pairs(g,"alignment",p)["alignment_jobs"],1)
        Path(g["A001"]).write_text(">contig\n"+"C"*301+"\n")
        self.assertEqual(self.run_pairs(g,"changed",p)["alignment_jobs"],1)
    def test_failed_command_is_not_zero_or_cache_success(self):
        g=dict(list(self.genomes.items())[:2])
        with patch("run_mummer.tool_identity",return_value=self.tools),patch("run_mummer.command_run",side_effect=WorkflowError("failure")):
            result=run_mummer(g,self.root/"bad",self.root/"cache")
        self.assertEqual(result["comparisons"][0]["status"],"ERROR")
        self.assertIsNone(result["comparisons"][0]["snp_distance"])
        self.assertEqual(self.run_pairs(g,"retry")["alignment_jobs"],1)
    def test_corrupt_cached_artifact_recomputes(self):
        g=dict(list(self.genomes.items())[:2]);a=self.run_pairs(g)
        (Path(a["comparisons"][0]["artifact_directory"])/"snps.tsv").write_text("corrupt")
        self.assertEqual(self.run_pairs(g,"retry")["alignment_jobs"],1)
    def test_interval_union_does_not_double_count_overlap(self):
        self.assertEqual(union_length([(1,10),(8,20),(30,35)]),26)
    def test_local_manifest_checksum_failure_never_downloads(self):
        p=self.root/"manifest.tsv";g=Path(self.genomes["A001"])
        p.write_text("accession\tpath\tsha256\nA001\t"+str(g)+"\t"+sha256_file(g)+"\n")
        with patch("fetch_reference_genomes.require_executable") as exe:
            result=fetch_genomes(["A001"],self.root/"fetch",1,0,p)
            self.assertEqual(result["local_assemblies"],["A001"]);exe.assert_not_called()
            g.write_text(">changed\nACGT\n")
            with self.assertRaises(WorkflowError):fetch_genomes(["A001"],self.root/"badfetch",1,0,p)
            exe.assert_not_called()
    def test_ties_low_coverage_and_metric_conflict_remain_visible(self):
        def row(name,snps,covered):
            return comparison(self.root,name,snps,covered)
        targets=[{"accession":name,"cluster":"C1"} for name in "ABC"]
        result=interpret([row("A",0,500),row("B",3,1000),row("C",3,1000)],targets,"C1",self.policy)
        self.assertEqual(result["nearest_samples"],["A"]);self.assertEqual(result["excluded_comparisons"],[])
        self.assertEqual(result["nearest_sample"],"A");self.assertIsNone(result["newick_tree"])
        result=interpret([row("A",100,850),row("B",110,1000)],targets,"C1",self.policy)
        self.assertFalse(result["ranking_basis_conflict"]);self.assertEqual(result["nearest_sample"],"A")
        self.assertIsNone(result["nearest_snp_distance"])
        self.assertEqual(result["candidate_challenges"]["comparisons"][0]["shared_target_bases"],850)

@unittest.skipUnless(os.environ.get("ISOSCOUT_TOOL_TESTS")=="1","Set ISOSCOUT_TOOL_TESTS=1 with pinned tools")
class MummerIntegrationTests(unittest.TestCase):
    def test_known_substitutions_indels_ambiguity_and_reverse_complement(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);rng=random.Random(842)
            seq="".join(rng.choices("ACGT",k=150000))
            a=list(seq);a[10000]="A" if a[10000]!="A" else "C";a="".join(a)
            variants={"QUERY":seq,"ONE":a,"TIE":a,"RC":a.translate(str.maketrans("ACGT","TGCA"))[::-1],"INDEL":a[:30000]+"AGTCA"+a[30000:],"N":a[:50000]+"NNNNN"+a[50005:],"SHORT":a[:20000]}
            genomes={}
            for name,sequence in variants.items():
                p=root/(name+".fna");p.write_text(">contig\n"+sequence+"\n");genomes[name]=str(p)
            result=run_mummer(genomes,root/"run",root/"cache")
            self.assertEqual(result["status"],"PASS",result)
            by={row["sample"]:row for row in result["comparisons"]}
            for name in ("ONE","TIE","RC","INDEL","N"):
                self.assertEqual(by[name]["snp_distance"],1,(name,by[name]))
            self.assertEqual(by["INDEL"]["indel_bases"],5)
            self.assertEqual(by["N"]["ambiguous_difference_rows"],0)  # MUMmer suppresses N differences upstream.
            targets=[{"accession":name,"cluster":"C1"} for name in variants if name!="QUERY"]
            interpreted=interpret(result["comparisons"],targets,"C1",result["policy"])
            self.assertIn("SHORT",[x["sample"] for x in interpreted["coverage_blockers"]])
            self.assertEqual(interpreted["cluster_status"],"RESOLVED")
            self.assertIsNone(interpreted["nearest_sample"])
            again=run_mummer(genomes,root/"again",root/"cache")
            self.assertEqual(again["alignment_jobs"],0)
            self.assertEqual(again["cache_hits"],6)
