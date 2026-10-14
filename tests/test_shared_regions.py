"""Candidate challenges: fair positions, retained uncertainty and bounded work."""
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import WorkflowError,load_json
from interpret_assembly_comparison import interpret
from shared_target_regions import intersect,merge,subtract
from site_fixture import MashpitFixture,comparison
from screen_isolate import run_snp_resolution

class SharedRegionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.policy=load_json(ROOT/'config/assembly-comparison-policy.json')
        self.policy['ranking']['min_shared_bases']=100  # Scale the evidence floor to the 1 kb software fixture.
    def call(self,rows,clusters=None):
        clusters=clusters or ['C1']*len(rows)
        return interpret(rows,[{'accession':r['sample'],'cluster':c} for r,c in zip(rows,clusters)],'C1',self.policy)
    def test_different_regions_reverse_count_order_without_new_alignments(self):
        a=comparison(self.root,'A',20,900)
        positions=[['ref',i] for i in range(1,11)]+[['ref',i] for i in range(950,980)]
        b=comparison(self.root,'B',positions=positions)
        d=self.call([a,b],['C1','C2'])
        self.assertEqual(d['nearest_samples'],['B']);self.assertEqual(d['nearest_cluster'],'C2')
        self.assertIsNone(d['nearest_snp_distance'])  # Different challenges cannot share one distance scalar.
        self.assertEqual(d['candidate_challenges']['comparisons'][0]['shared_target_bases'],900)
        self.assertEqual(b['snp_distance'],40)
    def test_coverage_below_reporting_reference_remains_and_can_win(self):
        d=self.call([comparison(self.root,'A',1,800),comparison(self.root,'B',40)],['C1','C2'])
        self.assertEqual(d['cluster_status'],'RESOLVED');self.assertEqual(d['nearest_cluster'],'C1')
        self.assertEqual(d['nearest_sample'],'A');self.assertEqual(d['excluded_comparisons'],[])
        self.assertTrue(next(r for r in d['ranked'] if r['sample']=='A')['coverage_flags'])
    def test_cycle_cannot_be_reported_as_resolved(self):
        rows=[comparison(self.root,'A',covered=200,regions={'ref':[[1,100],[201,300]]},positions=[['ref',201],['ref',202]]),
              comparison(self.root,'B',covered=200,regions={'ref':[[1,200]]},positions=[['ref',10]]),
              comparison(self.root,'C',covered=200,regions={'ref':[[101,300]]},positions=[['ref',110]])]
        d=self.call(rows,['C1','C2','C3'])
        self.assertEqual(d['genome_status'],'AMBIGUOUS');self.assertIsNone(d['nearest_cluster'])
        self.assertIsNone(d['nearest_sample'])
        self.assertLessEqual(d['candidate_challenges']['comparison_count'],d['candidate_challenges']['comparison_upper_bound'])
    def test_same_cluster_cycle_supports_cluster_but_never_unique_genome(self):
        rows=[comparison(self.root,'A',covered=200,regions={'ref':[[1,100],[201,300]]},positions=[['ref',201],['ref',202]]),
              comparison(self.root,'B',covered=200,regions={'ref':[[1,200]]},positions=[['ref',10]]),
              comparison(self.root,'C',covered=200,regions={'ref':[[101,300]]},positions=[['ref',110]]),
              comparison(self.root,'D',covered=300,positions=[['ref',i] for i in range(1,301)])]
        d=self.call(rows,['C1','C1','C1','C2'])
        self.assertEqual(d['genome_status'],'AMBIGUOUS')
        self.assertEqual(d['cluster_certificate']['status'],'VERIFIED')
        self.assertEqual(d['nearest_cluster'],'C1');self.assertIsNone(d['nearest_sample'])
        self.assertEqual(d['genome_status'],'AMBIGUOUS')
        self.assertTrue(d['cluster_certificate']['witnesses'])
        self.assertNotIn('D',d['decision_reference_samples'])
    def test_common_finalist_mask_resolves_a_cycle_with_shared_core(self):
        rows=[comparison(self.root,'A',covered=400,regions={'ref':[[1,300],[401,500]]},positions=[['ref',i] for i in range(210,220)]),
              comparison(self.root,'B',covered=400,regions={'ref':[[1,400]]},positions=[['ref',1]]+[['ref',i] for i in range(310,320)]),
              comparison(self.root,'C',covered=400,regions={'ref':[[1,200],[301,500]]},positions=[['ref',1],['ref',2]]+[['ref',i] for i in range(410,420)])]
        d=self.call(rows,['C1','C2','C3'])
        self.assertEqual(d['nearest_samples'],['A']);self.assertIsNone(d['nearest_cluster'])
        self.assertIsNone(d['nearest_sample'])
        self.assertEqual(d['candidate_challenges']['region_sensitivity']['status'],'REGION_SENSITIVE')
        self.assertIn('B',d['candidate_challenges']['region_sensitivity']['alternatives'])
        self.assertEqual(d['finalist_panel']['shared_target_bases'],200)
        self.assertEqual(d['finalist_panel']['snp_counts'],{'A':0,'B':1,'C':2})
        reverse=self.call(list(reversed(rows)),['C3','C2','C1'])
        self.assertEqual(d['candidate_challenges'],reverse['candidate_challenges'])
    def test_distant_partial_overlap_is_ruled_out_by_its_snp_lower_bound(self):
        d=self.call([comparison(self.root,'A',3),comparison(self.root,'B',200,400)],['C1','C2'])
        self.assertEqual(d['nearest_sample'],'A');self.assertEqual(d['nearest_cluster'],'C1')
        bound=d['candidate_challenges']['verification_passes'][-1]['candidate_bounds'][0]
        self.assertLess(bound['shared_fraction_of_larger_span'],.5)
        self.assertEqual(bound['disposition'],'RULED_OUT_BY_SNP_LOWER_BOUND')
        self.assertEqual(bound['observed_snps_lower_bound'],200)
    def test_close_partial_overlap_is_not_silently_ruled_out(self):
        d=self.call([comparison(self.root,'A',3),comparison(self.root,'B',0,400)],['C1','C2'])
        self.assertIsNone(d['nearest_sample']);self.assertIsNone(d['nearest_cluster'])
        self.assertEqual(d['decision_reference_samples'],['A','B'])
    def test_previously_ruled_out_reference_is_rechecked_after_mask_changes(self):
        rows=[comparison(self.root,'A',positions=[['ref',i] for i in range(850,860)]),
              comparison(self.root,'B',covered=900,positions=[['ref',i] for i in range(860,880)]),
              comparison(self.root,'C',covered=800,positions=[['ref',200]])]
        d=self.call(rows)
        steps=d['candidate_challenges']['verification_passes']
        first=next(r for r in steps[0]['candidate_bounds'] if r['sample']=='B')
        self.assertEqual(first['disposition'],'RULED_OUT_BY_SNP_LOWER_BOUND')
        self.assertIn('B',steps[1]['proposed_additions'])
        self.assertEqual(d['nearest_samples'],['A','B'])
        self.assertEqual(d['finalist_panel']['snp_counts'],{'A':0,'B':0,'C':1})
    def test_bounded_challenges_do_not_build_an_all_pairs_matrix(self):
        rows=[comparison(self.root,'A%03d'%i,i) for i in range(20)]
        d=self.call(rows)
        self.assertEqual(d['nearest_sample'],'A000')
        self.assertLessEqual(d['candidate_challenges']['comparison_count'],3*(len(rows)-1))
    def test_corrupt_position_data_fails_closed(self):
        a=comparison(self.root,'A');Path(a['site_evidence']['path']).write_text('{}')
        with self.assertRaises(WorkflowError):self.call([a])
    def test_one_cluster_with_exact_genome_ties(self):
        d=self.call([comparison(self.root,'A',3),comparison(self.root,'B',3)])
        self.assertEqual((d['cluster_status'],d['genome_status']),('RESOLVED','AMBIGUOUS'))
    def test_cross_cluster_ties_are_not_broken_by_order(self):
        a=comparison(self.root,'A',3);b=comparison(self.root,'B',3)
        d=self.call([b,a],['C2','C1']);reverse=self.call([a,b],['C1','C2'])
        self.assertEqual(d['cluster_status'],'AMBIGUOUS');self.assertIsNone(d['nearest_cluster'])
        self.assertEqual(d['nearest_samples'],reverse['nearest_samples'])
    def test_failed_candidate_remains_an_unresolved_alternative(self):
        error={'sample':'B','status':'ERROR','error':'fixture failure','snp_distance':None}
        d=self.call([comparison(self.root,'A',1),error],['C1','C2'])
        self.assertIsNone(d['nearest_cluster']);self.assertEqual(d['decision_reference_samples'],['A','B'])
    def test_interval_boundaries_and_multiple_contigs(self):
        self.assertEqual(merge([[1,10],[10,20],[21,22],[40,50]]),[(1,22),(40,50)])
        self.assertEqual(intersect({'a':[(1,10)],'b':[(4,8)]},{'a':[(10,12)],'b':[(1,4)]}),{'a':[(10,10)],'b':[(4,4)]})
        self.assertEqual(subtract({'a':[(1,10),(20,30)],'b':[(1,5)]},
                                  {'a':[(3,5),(8,22)],'b':[(1,10)]}),{'a':[(1,2),(6,7),(23,30)]})
    def test_workflow_keeps_a_low_coverage_winner(self):
        case=MashpitFixture(self.root/'workflow_case')
        case.fixture([1.,.999],['C1','C2'])
        rows=[comparison(self.root,'GCA_0000',1,800),comparison(self.root,'GCA_0001',40)]
        measurement={'comparisons':rows,'commands':[],'alignment_jobs':2,'cache_hits':0,'elapsed_seconds':0,'status':'PASS','tools':{},'policy':self.policy,'policy_sha256':'fixture','implementation_sha256':'fixture','cache_directory':str(self.root)}
        downloaded={'commands':[],'verified':{r['sample']:str(self.root/(r['sample']+'.fna')) for r in rows},'unavailable':[]}
        with patch('screen_isolate.fetch_genomes',return_value=downloaded),patch('screen_isolate.run_minimap',return_value=measurement):
            out=run_snp_resolution(case.source,self.root/'query.fna','C1',self.root/'workflow',case.database,expand=False,snp_backend='minimap2')['snp_resolution']
        self.assertEqual(out['nearest_sample'],'GCA_0000');self.assertEqual(out['cluster_status'],'RESOLVED')
        self.assertFalse(out['search_scope']['global_nearest_established'])
if __name__=='__main__':unittest.main()
