"""Reads must reach SNP comparison only through their QC-checked assembly."""
from contextlib import ExitStack
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from common import load_json, sha256_file, WorkflowError
from comparison_fixture import comparison_result
from screen_isolate import screen
from site_fixture import MashpitFixture, comparison


class ReadRefinementTests(unittest.TestCase):
    def test_generated_assembly_reaches_backend_and_provenance(self):
        for backend, fail, screen_only, qc_fail in [
            (None, False, False, False), ('mummer', False, False, False),
            (None, True, False, False), (None, False, True, False),
            (None, False, False, True),
        ]:
            with self.subTest(backend=backend, fail=fail, screen_only=screen_only, qc_fail=qc_fail), tempfile.TemporaryDirectory() as tmp, ExitStack() as stack:
                root = Path(tmp)
                case = MashpitFixture(root)
                case.fixture([1.])
                reads = [root / 'sample_R1.fastq', root / 'sample_R2.fastq']
                for read in reads:
                    read.write_text('@sample\nACGT\n+\nIIII\n')
                assembly = root / 'assembled.fna'
                assembly.write_text('>query\nACGT\n')
                reference = root / 'reference.fna'
                reference.write_text('>ref\nACGT\n')
                chosen = backend or 'minimap2'
                policy = load_json(ROOT / 'config' / ('assembly-comparison-policy.json' if chosen == 'minimap2' else 'mummer-comparison-policy.json'))
                policy['ranking']['min_shared_bases'] = 100
                comparison_run = comparison_result(policy)
                comparison_run['comparisons'] = [comparison(root, 'GCA_0000', target_hash=sha256_file(assembly))]
                values = {
                    'run_workflow': {'commands': [], 'assembly_path': str(assembly),
                        'read_paths': [str(root / 'clean_R1.fastq'), str(root / 'clean_R2.fastq')],
                        'read_qc': {'status': 'PASS', 'total_bases': 400}},
                    'assess': {'status': 'FAIL' if qc_fail else 'PASS',
                        'metrics': {'total_length': 4}, 'failures': ['test QC stop'] if qc_fail else []},
                    'run_mashpit': {'command': ['mashpit'], 'database': case.database, 'output_directory': str(case.source)},
                    'fetch_genomes': {'commands': [], 'verified': {'GCA_0000': str(reference)}, 'unavailable': []},
                    'collect': {},
                }
                mocks = {name: stack.enter_context(patch('screen_isolate.' + name, return_value=value)) for name, value in values.items()}
                runner = stack.enter_context(patch('screen_isolate.run_' + ('minimap' if chosen == 'minimap2' else 'mummer'),
                    return_value=comparison_run, side_effect=WorkflowError('comparison failed') if fail else None))
                unused = stack.enter_context(patch('screen_isolate.run_' + ('mummer' if chosen == 'minimap2' else 'minimap')))
                stack.enter_context(patch('analyze_similarity_distribution.render_plot', return_value={'status': 'SKIPPED'}))
                stack.enter_context(patch('plot_assembly_comparison.render', return_value={'status': 'SKIPPED'}))
                result = screen(reads, root / 'output', root, organism='listeria', snp_resolve=not screen_only, snp_backend=backend)
                unused.assert_not_called()
                self.assertEqual(result['input_type'], 'illumina_paired_fastq')
                if screen_only or qc_fail:
                    runner.assert_not_called()
                    mocks['fetch_genomes'].assert_not_called()
                    self.assertEqual(result['status'], 'STOPPED' if qc_fail else 'COMPLETED_WITH_WARNINGS')
                    continue
                runner.assert_called_once()
                self.assertEqual(runner.call_args.args[0], {'QUERY': str(assembly), 'GCA_0000': str(reference)})
                snp = result['snp_resolution']
                self.assertEqual(snp['backend'], chosen)
                self.assertEqual(snp['query_input'], {'type': 'assembly', 'paths': [str(assembly)]})
                query = next(row for row in snp['file_checksums'] if row['sample'] == 'QUERY')
                self.assertEqual(query['sha256'], sha256_file(assembly))
                self.assertEqual(result['status'], 'INCOMPLETE' if fail else 'COMPLETED_WITH_WARNINGS')
                report = (root / 'output/report.md').read_text()
                self.assertIn('assembly generated from paired-end reads', report)
                self.assertNotIn('cleaned paired-end reads', report)
                self.assertNotIn('SNP-distance tree', report)
                if fail:
                    self.assertIn('comparison failed', report)
