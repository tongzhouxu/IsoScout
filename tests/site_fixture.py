"""Small explicit target-coordinate examples for interpretation tests."""
from common import write_json,sha256_file

def comparison(root, name, count=0, covered=1000, regions=None, positions=None, target_hash='same-target'):
    regions = regions if regions is not None else {'ref': [[1, covered]]}
    positions = positions if positions is not None else [['ref', i] for i in range(1, count+1)]
    evidence = {'schema':'target-sites-v2','coordinates':'one_based_inclusive',
                'target_sha256':target_hash,'candidate_sha256':name,
                'target_lengths':{'ref':1000},'regions':regions,'snps':positions}
    path=root/(name+'.json');write_json(path,evidence)
    return {'sample':name,'status':'COMPARED','snp_distance':len(positions),
            'target_aligned_bases':covered,'candidate_aligned_bases':covered,
            'target_comparable_bases':covered,
            'target_bases':1000,'candidate_bases':1000,
            'target_aligned_fraction':covered/1000,'candidate_aligned_fraction':covered/1000,
            'snps_per_target_aligned_mb':len(positions)*1e6/covered,
            'target_sha256':target_hash,'candidate_sha256':name,'indel_bases':0,
            'site_evidence':{'path':str(path),'sha256':sha256_file(path),'schema':'target-sites-v2'}}
