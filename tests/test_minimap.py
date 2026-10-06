"""Pinned minimap2/paftools integration and parser contracts."""
from __future__ import annotations
import os,random,sys,tempfile,unittest,json
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from run_minimap import run_minimap,parse_calls,prepare_paf
from interpret_assembly_comparison import interpret
from common import WorkflowError,load_json

@unittest.skipUnless(os.environ.get("ISOSCOUT_TOOL_TESTS")=="1","Set ISOSCOUT_TOOL_TESTS=1 with pinned tools")
class MinimapIntegrationTests(unittest.TestCase):
    def test_known_substitutions_indels_ambiguity_and_reverse_complement(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);rng=random.Random(842)
            seq="".join(rng.choices("ACGT",k=150000))
            a=list(seq);a[10000]="A" if a[10000]!="A" else "C";a="".join(a)
            variants={"QUERY":seq,"ONE":a,"TIE":a,"RC":a.translate(str.maketrans("ACGT","TGCA"))[::-1],"INDEL":a[:30000]+"AGTCA"+a[30000:],"N":a[:50000]+"NNNNN"+a[50005:],"SHORT":a[:20000]}
            genomes={}
            for name,sequence in variants.items():
                p=root/(name+".fna");p.write_text(">contig\n"+sequence+"\n");genomes[name]=str(p)
            result=run_minimap(genomes,root/"run",root/"cache")
            self.assertEqual(result["status"],"PASS",result)
            by={row["sample"]:row for row in result["comparisons"]}
            for name in ("ONE","TIE","RC","INDEL","N"):
                self.assertEqual(by[name]["snp_distance"],1,(name,by[name]))
            self.assertEqual(by["INDEL"]["indel_bases"],5)
            self.assertEqual(by["N"]["ambiguous_difference_rows"],0)  # paftools suppresses N differences upstream.
            targets=[{"accession":name,"cluster":"C1"} for name in variants if name!="QUERY"]
            interpreted=interpret(result["comparisons"],targets,"C1",result["policy"])
            self.assertIn("SHORT",[x["sample"] for x in interpreted["excluded_comparisons"]])
            self.assertIsNone(interpreted["nearest_sample"])
            again=run_minimap(genomes,root/"again",root/"cache")
            self.assertEqual(again["alignment_jobs"],0)
            self.assertEqual(again["cache_hits"],6)

class CallParserTests(unittest.TestCase):
    def test_overlapping_variants_are_excluded_and_unique_regions_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"calls.tsv"
            p.write_text("R\tref\t0\t100\nV\tref\t10\t11\t1\t60\ta\tc\tquery\t10\t11\t+\nV\tref\t20\t21\t2\t60\ta\tt\tquery\t20\t21\t+\n")
            row=parse_calls(p,{"ref":200})
            self.assertEqual(row["snp_distance"],1)
            self.assertEqual(row["target_aligned_bases"],100)
    def test_malformed_calls_fail_instead_of_returning_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"calls.tsv";p.write_text("invalid\n")
            with self.assertRaises(WorkflowError):parse_calls(p,{"ref":200})


class MinimapCacheTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.policy = load_json(ROOT / 'config/assembly-comparison-policy.json')
        self.launches = []
        self.calls = []
        for tool in ('minimap2', 'paftools.js', 'k8'):
            (self.root / tool).write_text('fixture tool')
        self.genomes = {}
        for i in range(126):
            name = 'QUERY' if i == 0 else f'A{i:03d}'
            path = self.root / (name + '.fna')
            path.write_text('>contig\n' + 'A' * (1000 + i) + '\n')
            self.genomes[name] = str(path)

    def execute(self, cmd, directory, output, timeout):
        if Path(cmd[0]).name == 'minimap2':
            self.launches.append(cmd)
            n = len(Path(cmd[-1]).read_text().splitlines()[1])
            (directory / output).write_text(
                f'contig\t{n}\t0\t1000\t+\tcontig\t1000\t0\t1000\t1000\t1000\t60\ttp:A:P\tcs:Z::1000\n')
        else:
            self.calls.append(cmd)
            (directory / output).write_text('R\tcontig\t0\t1000\n')

    def run_pairs(self, genomes, policy=None):
        with patch('run_minimap.require_executable', side_effect=lambda name: str(self.root / name)), \
             patch('run_minimap.subprocess.run', return_value=SimpleNamespace(stdout='2.31-r1302')), \
             patch('run_minimap.command_run', side_effect=self.execute):
            return run_minimap(genomes, self.root / 'out', self.root / 'cache', policy)

    def test_expansion_and_resume_align_only_new_target_pairs(self):
        first = dict(list(self.genomes.items())[:101])
        self.assertEqual(self.run_pairs(first)['alignment_jobs'], 100)
        expanded = self.run_pairs(self.genomes)
        self.assertEqual(expanded['alignment_jobs'], 25)
        self.assertEqual(expanded['cache_hits'], 100)
        resumed = self.run_pairs(self.genomes)
        self.assertEqual(resumed['alignment_jobs'], 0)
        self.assertEqual(resumed['cache_hits'], 125)
        self.assertEqual(len(self.launches), 125)
        self.assertTrue(all(cmd[-2] == self.genomes['QUERY'] for cmd in self.launches))

    def test_filters_recall_without_realignment_and_input_change_invalidates(self):
        genomes = dict(list(self.genomes.items())[:2])
        self.run_pairs(genomes)
        policy_path = self.root / 'policy.json'
        self.policy['comparability']['min_target_aligned_fraction'] = .9
        policy_path.write_text(json.dumps(self.policy))
        self.assertEqual(self.run_pairs(genomes, policy_path)['cache_hits'], 1)
        self.assertEqual(len(self.calls), 1)
        self.policy['variant_calling']['min_alignment_length'] = 600
        policy_path.write_text(json.dumps(self.policy))
        self.assertEqual(self.run_pairs(genomes, policy_path)['alignment_jobs'], 0)
        self.assertEqual(len(self.calls), 2)
        Path(genomes['A001']).write_text('>contig\n' + 'C' * 1001 + '\n')
        self.assertEqual(self.run_pairs(genomes, policy_path)['alignment_jobs'], 1)

    def test_corrupt_alignment_recomputes_and_failures_are_not_zero(self):
        genomes = dict(list(self.genomes.items())[:2])
        result = self.run_pairs(genomes)
        directory = Path(result['comparisons'][0]['artifact_directory'])
        (directory / 'alignment.paf').write_text('corrupt')
        original_execute = self.execute
        def fail(*args):
            raise WorkflowError('command failed')
        self.execute = fail
        result = self.run_pairs(genomes)
        self.assertEqual(result['status'], 'PARTIAL')
        self.assertIsNone(result['comparisons'][0]['snp_distance'])
        self.assertFalse((directory / 'alignment.json').exists())
        self.execute = original_execute
        self.assertEqual(self.run_pairs(genomes)['alignment_jobs'], 1)


class CoverageLookupTests(unittest.TestCase):
    def test_unsorted_overlapping_regions_preserve_boundary_membership(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / 'calls.tsv'
            regions = 'R\tref\t40\t60\nR\tref\t10\t90\nR\tref\t120\t150\n'
            variants = ''.join(f'V\tref\t{pos}\t{pos+1}\t1\t60\ta\tc\tquery\t{pos}\t{pos+1}\t+\n'
                               for pos in (9,10,40,89,90,119,120,149,150))
            p.write_text(regions + variants)
            row = parse_calls(p, {'ref':200})
            self.assertEqual(row['snp_distance'], 5)
            self.assertEqual(row['target_aligned_bases'], 110)
