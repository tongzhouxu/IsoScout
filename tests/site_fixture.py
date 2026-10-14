"""Small reusable fixtures shared by the workflow tests."""
import csv
from pathlib import Path

from common import write_json,sha256_file


class MashpitFixture:
    """Create the minimal Mashpit result files used by selection tests."""

    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.source = self.root / "mashpit"
        self.source.mkdir()
        self.database = {
            "name": "listeria", "version": "test",
            "mashpit_database_settings": {"hash_number": 1000},
        }
        self.profile = {"number": 200, "threshold": .85, "tie_tolerance_hashes": 2}
        self.policy = {
            "policy_version": "3.0.1", "strategy": "global-ranked-cluster-coverage-v1",
            "default_mode": "adaptive", "coverage_strategy": "unrepresented_cluster_leaders_then_global_rank",
            "initial_target_genomes": 50, "max_returned_representatives": 100,
            "max_total_genomes": 100,
        }
        self.workflow = {
            "mashpit": self.profile,
            "similarity_diagnostics": {"rank_checkpoints": [50, 100, 200]},
        }

    def fixture(self, scores, clusters=None):
        clusters = clusters or ["C1"] * len(scores)
        with (self.source / "sample_representative_matches.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["asm_acc", "PDS_acc", "similarity_score"])
            writer.writeheader()
            for index, (score, cluster) in enumerate(zip(scores, clusters)):
                writer.writerow({"asm_acc": f"GCA_{index:04d}", "PDS_acc": cluster, "similarity_score": score})
        with (self.source / "sample_cluster_candidates.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["PDS_acc", "best_similarity_score", "near_top"])
            writer.writeheader()
            for cluster in sorted(set(clusters)):
                score = max(score for score, group in zip(scores, clusters) if cluster == group)
                writer.writerow({"PDS_acc": cluster, "best_similarity_score": score, "near_top": max(scores) - score <= .002})

    def analyze(self):
        from analyze_similarity_distribution import analyze
        return analyze(self.source, self.database, self.profile, self.policy, [50, 100, 200])

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
