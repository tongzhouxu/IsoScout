"""Display final-anchor challenges; counts are comparable within each pair only."""
from pathlib import Path


def render(result: dict, output: Path) -> dict:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        audit = result.get("candidate_challenges") or {}
        anchor = audit.get("anchor")
        rows = [r for r in audit.get("comparisons", []) if anchor in (r["first"], r["second"])]
        if not rows:
            return {"status": "SKIPPED", "reason": "No anchor challenges"}
        x = list(range(1, len(rows) + 1))
        anchor_counts = [r["first_snps"] if r["first"] == anchor else r["second_snps"] for r in rows]
        other_counts = [r["second_snps"] if r["first"] == anchor else r["first_snps"] for r in rows]
        colors = ["#d97706" if r["relation"] == "INSUFFICIENT_DATA" else "#3478a3" for r in rows]
        fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
        axes[0].scatter(x, [b-a for a, b in zip(anchor_counts, other_counts)], c=colors, s=14)
        axes[0].axhline(0, color="gray", lw=1)
        axes[0].set_ylabel("Challenger SNPs − anchor SNPs")
        axes[0].set_yscale("symlog", linthresh=10)
        axes[0].set_title(f"Anchor {anchor}: {audit.get('status')}\nPositive: fewer anchor SNPs on that pair's shared positions")
        axes[1].scatter(x, [r["shared_target_bases"]/1e6 for r in rows], c=colors, s=14)
        axes[1].set_ylabel("Shared target sequence (Mb)")
        axes[1].set_xlabel("Challenge index (not distance rank); orange = insufficient evidence")
        fig.tight_layout()
        fig.savefig(output, dpi=150)
        plt.close(fig)
        return {"status": "PASS", "path": str(output)}
    except (ImportError, ValueError, OSError) as error:
        return {"status": "UNAVAILABLE", "reason": str(error)}
