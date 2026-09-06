"""
Regenerate Fig_family_dispersion_threshold_sweep_rev.png using the shared seeded
engine (family_lib), so the figure matches Table_family_threshold_sweep_rev.tsv
exactly. Identity sweep extended to 20 and 40 (7 panels).

The x-axis is the standard deviation of log10(half-life); an SD of a log is
invariant to the time unit, so there is no hours/days distinction on this axis
(unlike the half-life-axis figures).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import family_lib as fl

IDENTITIES = [20, 30, 40, 50, 60, 70, 80]
OUTDIR = os.environ.get("FIG_OUTDIR", "/mnt/workspace/output")
BINS = np.linspace(0.0, 0.7, 36)


def main():
    hl, pairs = fl.load()
    pool = np.log10(np.array(list(hl.values())))
    fig, axes = plt.subplots(2, 4, figsize=(18, 8.5))
    axes = axes.ravel()
    fig.suptitle("Within-family half-life dispersion vs permuted null (threshold sweep)",
                 fontsize=14, fontweight="bold")

    for ax, thr in zip(axes, IDENTITIES):
        r = fl.analyse(thr, hl, pairs, pool)
        # Two histograms per panel: random null (grey) and within-family (green).
        # Series legend and the fam-med/rand-med/p annotation are omitted per request;
        # colours are defined in the figure caption.
        ax.hist(r["pool_sds"], bins=BINS, density=True, color="0.6", alpha=0.65)
        ax.hist(r["within"], bins=BINS, density=True, color="tab:green", alpha=0.7)
        ax.set_title("identity ≥ %d%%, coverage ≥ 50%%" % thr, fontsize=10)
        ax.set_xlabel("Std. of log₁₀(half-life) within set")
        ax.set_ylabel("Density")

    for ax in axes[len(IDENTITIES):]:
        ax.axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(OUTDIR, "Fig_family_dispersion_threshold_sweep_rev.png")
    fig.savefig(out, dpi=190)
    print("wrote", out)


if __name__ == "__main__":
    main()
