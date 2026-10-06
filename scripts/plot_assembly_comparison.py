"""Display common-panel counts and separately labeled partial-evidence bounds."""
from pathlib import Path


def render(result: dict, output: Path) -> dict:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        audit = result.get("candidate_challenges") or {}
        panel = audit.get("finalist_panel") or {}
        scores = panel.get("snp_counts") or {}
        if not scores:
            return {"status": "SKIPPED", "reason": "No sufficient common finalist region"}
        names = sorted(scores)
        steps = audit.get("verification_passes") or []
        bounds = steps[-1].get("candidate_bounds", []) if steps else []
        sensitive = set((audit.get("region_sensitivity") or {}).get("alternatives", []))
        fig, axes = plt.subplots(2, 1, figsize=(10, 7), gridspec_kw={"height_ratios": [1, 1]})
        axes[0].bar(range(len(names)), [scores[n] for n in names],
                    color=["#d97706" if n in sensitive else "#3478a3" for n in names])
        axes[0].set_xticks(range(len(names)))
        if len(names) <= 15:
            axes[0].set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        else:
            axes[0].set_xticklabels([])
            axes[0].set_xlabel("Finalist index in accession order; labels in the audit")
        axes[0].set_ylabel("SNPs on identical positions")
        axes[0].set_title(f"{panel['shared_target_bases']:,} common aligned target bases\n"
                          f"Cluster: {result.get('cluster_status')}; genome: {result.get('genome_status')}")
        colors = ["#3478a3" if r["disposition"] == "RULED_OUT_BY_SNP_LOWER_BOUND" and r["sample"] not in sensitive
                  else "#d97706" for r in bounds]
        axes[1].scatter(range(1, len(bounds)+1), [r["observed_snps_lower_bound"] for r in bounds], c=colors, s=15)
        axes[1].axhline(panel["minimum_snps"], color="black", ls="--", lw=1, label="Finalist minimum")
        axes[1].set_yscale("symlog", linthresh=10)
        axes[1].set_ylabel("Observed SNP lower bound\n(partial regions; not a distance)")
        axes[1].set_xlabel("Outside-reference index; orange = retained or contradictory evidence")
        axes[1].legend(loc="best")
        if not bounds:
            axes[1].text(.5, .5, "All examined references are in the finalist panel", transform=axes[1].transAxes, ha="center")
        fig.tight_layout()
        fig.savefig(output, dpi=150)
        plt.close(fig)
        return {"status": "PASS", "path": str(output)}
    except (ImportError, ValueError, OSError) as error:
        return {"status": "UNAVAILABLE", "reason": str(error)}
