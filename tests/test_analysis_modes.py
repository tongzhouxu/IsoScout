"""Mode and completion regressions; no network or biological tools required."""
from contextlib import ExitStack
from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import screen_isolate
from common import WorkflowError
from generate_report import _snp_section, generate_report
from launch_isoscout import parse_args


class AnalysisModeTests(unittest.TestCase):
    def test_default_requires_snp_and_screen_only_does_not(self):
        cases = [
            ({'status': 'COMPARED'}, None, True, 'COMPLETED'),
            ({'status': 'AMBIGUOUS'}, None, True, 'COMPLETED'),
            ({'status': 'INSUFFICIENT_DATA'}, None, True, 'INCOMPLETE'),
            ({'status': 'COMPARED', 'unavailable': ['missing']}, None, True, 'INCOMPLETE'),
            ({'status': 'AMBIGUOUS', 'query_comparisons': [{'status': 'ERROR'}]}, None, True, 'INCOMPLETE'),
            (None, WorkflowError('minimap2 not found'), True, 'INCOMPLETE'),
            (None, None, False, 'COMPLETED'),
        ]
        cases = [(*case, True) for case in cases] + [(None, None, True, 'INCOMPLETE', False)]
        for snp, error, full, expected, eligible in cases:
            with self.subTest(snp=snp, error=error, full=full), tempfile.TemporaryDirectory() as tmp, ExitStack() as stack:
                work = Path(tmp); sample = work/'sample.fna'; sample.write_text('>sample\nACGT\n')
                qc = {'status': 'PASS', 'metrics': {'total_length': 4}, 'warnings': []}
                data = {'name':'cronobacter', 'version':'test'}
                parsed = {'best_candidate': {'cluster':'PDS1','score':1.0} if eligible else None, 'screening_result':'candidate', 'warnings':[]}
                values = {
                    'inspect': {'input_type':'assembly','assembly_path':str(sample)}, 'assess':qc,
                    'run_mashpit': {'command':['mashpit'],'database':data,'output_directory':str(work/'mashpit')},
                    'locate_candidate_file':work/'candidates.csv', 'load_candidates':[], 'interpret':parsed,
                    'locate_tree_image':None, 'run_diagnostics':{'warnings':[]},
                    'annotate_rep_designations':{'status':'TABLE_ABSENT'}, 'collect':{},
                    'summary_text':'summary', 'generate_report':'report',
                }
                for name,value in values.items():stack.enter_context(patch.object(screen_isolate,name,return_value=value))
                run=stack.enter_context(patch.object(screen_isolate,'run_snp_resolution',
                    return_value={'snp_resolution':snp,'commands':[]},side_effect=error))
                kwargs = {} if full else {'snp_resolve':False}
                result=screen_isolate.screen([sample],work/'out',work/'databases',organism='cronobacter',**kwargs)
                self.assertEqual(result['status'],expected)
                self.assertEqual(run.call_count,int(full and eligible))
                self.assertEqual(result['analysis_mode'],'full' if full else 'screen_only')
                saved=json.loads((work/'out/result.json').read_text())
                self.assertEqual(saved['status'],expected)
                timings = json.loads((work/'out/stage_timings.json').read_text())
                self.assertEqual(timings, saved['stage_timings'])
                self.assertEqual(timings['mashpit']['calls'], 1)
                self.assertGreaterEqual(timings['mashpit']['elapsed_seconds'], 0)

    def test_stage_timer_accumulates_and_preserves_exceptions(self):
        timings = {}
        with patch.object(screen_isolate.time, 'monotonic', side_effect=[10., 12., 20., 23.]):
            with screen_isolate.measure_stage(timings, 'reference_retrieval'):
                pass
            with self.assertRaisesRegex(WorkflowError, 'download failed'):
                with screen_isolate.measure_stage(timings, 'reference_retrieval'):
                    raise WorkflowError('download failed')
        self.assertEqual(timings['reference_retrieval'],
                         {'elapsed_seconds': 5., 'calls': 2, 'failed_calls': 1})

    def test_unexecuted_report_does_not_claim_backend(self):
        # No synthetic comparison result should be manufactured for missing candidates.
        self.assertIn('Not run', '\n'.join(_snp_section({'analysis_mode':'full'})))
        screen_report='\n'.join(_snp_section({'analysis_mode':'screen_only'}))
        self.assertIn('screening-only',screen_report)
        self.assertNotIn('ska',screen_report)
        report = generate_report({'sample':'test', 'status':'INCOMPLETE',
                                  'analysis_mode':'full', 'stop_reason':'Missing tool',
                                  'snp_resolution':{'status':'ERROR','error':'Missing tool'}})
        self.assertIn('Step 3: SNP comparison', report)
        self.assertIn('Missing tool', report)

    def test_conflicting_screen_only_options_are_rejected(self):
        with patch('sys.stderr'):
            for options in [ ['--snp-resolve'], ['--snp-expand'], ['--snp-backend','minimap2'] ]:
                with self.subTest(options=options), self.assertRaises(SystemExit) as error:
                    parse_args(['query.fna','--output','out','--screen-only',*options])
                self.assertEqual(error.exception.code,2)

if __name__ == '__main__':
    unittest.main()
