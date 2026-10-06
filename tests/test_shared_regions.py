"""Behavioral tests for consistent regions, coverage exclusions and identity."""
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import WorkflowError,load_json
from interpret_assembly_comparison import interpret
from shared_target_regions import intersect,merge
from site_fixture import comparison
from screen_isolate import summary_text
from screen_isolate import run_snp_resolution
import test_similarity_distribution as fixtures

class SharedRegionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.policy=load_json(ROOT/'config/assembly-comparison-policy.json')
    def call(self,rows,clusters=None):
        clusters=clusters or ['C1']*len(rows)
        return interpret(rows,[{'accession':r['sample'],'cluster':c} for r,c in zip(rows,clusters)],'C1',self.policy)
    def test_different_regions_reverse_count_order_without_new_alignments(self):
        a=comparison(self.root,'A',20,900)
        positions=[['ref',i] for i in range(1,11)]+[['ref',i] for i in range(950,980)]
        b=comparison(self.root,'B',positions=positions)
        d=self.call([a,b],['C1','C2'])
        self.assertEqual(d['nearest_samples'],['B'])
        self.assertEqual(d['nearest_cluster'],'C2')
        self.assertEqual(d['nearest_snp_distance'],10)
        self.assertEqual(d['shared_regions']['shared_target_bases'],900)
        self.assertEqual(d['pair_specific_diagnostics']['minimum_count_samples'],['A'])
        self.assertEqual(b['snp_distance'],40)  # Original measurements stay intact.
    def test_low_coverage_alternative_blocks_wrong_cluster_resolution(self):
        d=self.call([comparison(self.root,'excluded',1,800),comparison(self.root,'eligible',40)],['C1','C2'])
        self.assertEqual(d['cluster_status'],'AMBIGUOUS')
        self.assertIsNone(d['nearest_cluster']);self.assertIsNone(d['nearest_sample'])
        self.assertEqual(d['cluster_candidates'],['C1','C2'])
        self.assertEqual(d['nearest_samples'],['eligible'])
        self.assertEqual(d['coverage_blockers'][0]['sample'],'excluded')
    def test_same_cluster_exclusion_blocks_only_genome_resolution(self):
        d=self.call([comparison(self.root,'excluded',1,800),comparison(self.root,'eligible',40)])
        self.assertEqual((d['cluster_status'],d['genome_status']),('RESOLVED','AMBIGUOUS'))
    def test_distant_excluded_reference_is_not_promoted(self):
        d=self.call([comparison(self.root,'excluded',200,800),comparison(self.root,'eligible',1)],['C1','C2'])
        self.assertEqual(d['nearest_cluster'],'C2');self.assertEqual(d['coverage_blockers'],[])
    def test_insufficient_intersection_does_not_fall_back_to_raw_counts(self):
        a=comparison(self.root,'A',0,900)
        b=comparison(self.root,'B',0,900,regions={'ref':[[101,1000]]})
        d=self.call([a,b],['C1','C2'])
        self.assertEqual(d['shared_regions']['shared_target_bases'],800)
        self.assertEqual(d['status'],'INSUFFICIENT_DATA')
        self.assertEqual(d['nearest_samples'],[]);self.assertIsNone(d['nearest_snp_distance'])
        text=summary_text({'sample':'fixture','status':'COMPLETED_WITH_WARNINGS','snp_resolution':d})
        self.assertIn('INSUFFICIENT_DATA',text)
    def test_missing_position_data_fails_closed(self):
        a=comparison(self.root,'A');a.pop('site_evidence')
        with self.assertRaises(WorkflowError):self.call([a])
    def test_corrupt_position_data_fails_closed(self):
        a=comparison(self.root,'A');Path(a['site_evidence']['path']).write_text('{}')
        with self.assertRaises(WorkflowError):self.call([a])
    def test_different_targets_cannot_share_a_mask(self):
        with self.assertRaises(WorkflowError):self.call([comparison(self.root,'A'),comparison(self.root,'B',target_hash='different')])
    def test_counts_must_match_position_evidence(self):
        a=comparison(self.root,'A',3);a['snp_distance']=4
        with self.assertRaises(WorkflowError):self.call([a])
    def test_one_cluster_with_exact_genome_ties(self):
        d=self.call([comparison(self.root,'A',3),comparison(self.root,'B',3)])
        self.assertEqual((d['cluster_status'],d['genome_status']),('RESOLVED','AMBIGUOUS'))
        self.assertEqual(d['nearest_by_aligned_snp_rate'],d['nearest_samples'])
    def test_cross_cluster_ties_are_not_broken_by_order(self):
        a=comparison(self.root,'A',3);b=comparison(self.root,'B',3)
        d=self.call([b,a],['C2','C1']);reverse=self.call([a,b],['C1','C2'])
        self.assertEqual(d['cluster_status'],'AMBIGUOUS');self.assertIsNone(d['nearest_cluster'])
        self.assertEqual(d['shared_regions']['mask_sha256'],reverse['shared_regions']['mask_sha256'])
        self.assertEqual(d['nearest_samples'],reverse['nearest_samples'])
    def test_no_qualifying_data_has_no_cluster_call(self):
        d=self.call([comparison(self.root,'A',0,500)])
        self.assertEqual(d['cluster_status'],'INSUFFICIENT_DATA');self.assertEqual(d['cluster_candidates'],[])
    def test_interval_boundaries_and_multiple_contigs(self):
        self.assertEqual(merge([[1,10],[10,20],[21,22],[40,50]]),[(1,22),(40,50)])
        self.assertEqual(intersect({'a':[(1,10)],'b':[(4,8)]},{'a':[(10,12)],'b':[(1,4)]}),{'a':[(10,10)],'b':[(4,4)]})
    def test_new_candidates_change_mask_and_counts_deterministically(self):
        a=comparison(self.root,'A',positions=[['ref',950]])
        first=self.call([a]);second=self.call([a,comparison(self.root,'B',0,900)])
        self.assertEqual(first['nearest_snp_distance'],1)
        self.assertEqual(second['nearest_snp_distance'],0)
        self.assertNotEqual(first['shared_regions']['mask_sha256'],second['shared_regions']['mask_sha256'])

    def test_workflow_stops_before_member_download_when_common_regions_are_insufficient(self):
        case=fixtures.SimilarityDiagnosticsTests();case.setUp()
        self.addCleanup(case.temporary.cleanup)
        case.fixture([1.,.999],['C1','C2'])
        rows=[comparison(self.root,'GCA_0000',0,900),
              comparison(self.root,'GCA_0001',0,900,regions={'ref':[[101,1000]]})]
        measurement={'comparisons':rows,'commands':[],'alignment_jobs':2,'cache_hits':0,
                     'elapsed_seconds':0,'status':'PASS','tools':{},'policy':self.policy,
                     'policy_sha256':'fixture','implementation_sha256':'fixture','cache_directory':str(self.root)}
        downloaded={'commands':[],'verified':{r['sample']:str(self.root/(r['sample']+'.fna')) for r in rows},'unavailable':[]}
        with patch('screen_isolate.fetch_genomes',return_value=downloaded) as fetch, \
             patch('screen_isolate.run_minimap',return_value=measurement), \
             patch('screen_isolate.discover_members') as members:
            out=run_snp_resolution(case.source,self.root/'query.fna','C1',self.root/'workflow',
                                   case.database,expand=True,snp_backend='minimap2')['snp_resolution']
        self.assertEqual(out['expansion']['stop_reason'],'insufficient_shared_target_regions')
        self.assertEqual(out['cluster_status'],'INSUFFICIENT_DATA')
        self.assertEqual(fetch.call_count,1);members.assert_not_called()
        self.assertFalse(out['search_scope']['global_nearest_established'])
if __name__=='__main__':unittest.main()
