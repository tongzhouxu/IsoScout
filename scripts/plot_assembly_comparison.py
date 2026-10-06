"""Plot target-to-candidate distances and alignment coverage."""
from pathlib import Path

def render(result: dict, output: Path) -> dict:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        rows=result.get("ranked",[])
        if not rows or result.get("nearest_snp_distance") is None: return {"status":"SKIPPED","reason":"No qualifying candidate comparisons"}
        fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True)
        x=list(range(1,len(rows)+1))
        axes[0].scatter(x,[r["ranking_snp_distance"] for r in rows],s=14)
        axes[0].set_ylabel("SNPs on shared target regions")
        axes[0].set_title("Target-to-candidate comparisons on identical target positions")
        axes[1].plot(x,[100*r["target_aligned_fraction"] for r in rows],label="Target aligned")
        axes[1].plot(x,[100*r["candidate_aligned_fraction"] for r in rows],label="Candidate aligned")
        axes[1].set_ylabel("Aligned (%)"); axes[1].set_xlabel("Candidate rank (ties retained)")
        axes[1].legend(); fig.tight_layout(); fig.savefig(output,dpi=150); plt.close(fig)
        return {"status":"PASS","path":str(output)}
    except (ImportError,ValueError,OSError) as error:
        return {"status":"UNAVAILABLE","reason":str(error)}
