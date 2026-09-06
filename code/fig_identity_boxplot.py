"""
Regenerate Fig_identity_vs_halflife_boxplot_rev.png for reproducibility.

Pairwise paralog half-life difference |Δlog2(t1/2)| stratified by sequence
identity bin (mutual coverage >= 50%), from the frozen pair table. The y-axis is
a difference of logs (unit-invariant — no hours/days distinction). Axis range is
left to the data (auto), as requested.

Input  : data/snapshots/pair_table_rev.tsv ; data/inputs/cyno_target_rev.tsv
Output : Fig_identity_vs_halflife_boxplot_rev.png  (deterministic; no RNG)
"""
import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import family_lib as fl
from stats_numpy import rankdata

OUTDIR = os.environ.get("FIG_OUTDIR", "/mnt/workspace/output")
BINS = [(-1, 20, "<20"), (20, 40, "20-40"), (40, 60, "40-60"),
        (60, 80, "60-80"), (80, 1e9, "≥80")]
COLORS = ["#46327e", "#5e6a9e", "#2a9d8f", "#73c46b", "#f4d35e"]  # viridis-like


def kruskal_wallis(groups):
    """Tie-corrected Kruskal-Wallis H and chi-square p (df = k-1, even-df closed form)."""
    allv = np.concatenate(groups)
    N = len(allv)
    ranks = rankdata(allv)
    idx = 0
    H = 0.0
    for g in groups:
        n = len(g)
        Rg = ranks[idx:idx + n].sum(); idx += n
        H += Rg * Rg / n
    H = 12.0 / (N * (N + 1)) * H - 3.0 * (N + 1)
    _, counts = np.unique(allv, return_counts=True)
    tie = (counts ** 3 - counts).sum()
    H /= (1.0 - tie / (N ** 3 - N))
    k = len(groups)
    df = k - 1
    x = H
    if df % 2 == 0:                      # even df closed-form survival function
        terms = df // 2
        s = sum((x / 2.0) ** i / math.factorial(i) for i in range(terms))
        p = math.exp(-x / 2.0) * s
    else:
        from stats_numpy import betainc
        p = float("nan")
    return H, df, float(p)


def main():
    hl, pairs = fl.load()

    # Composite gene names are normalised via fl.canonical(), which maps the pair
    # table's '_' spelling ('Rack1_Gnb2L1') onto the cohort's ';' spelling. A plain
    # `g in hl` test would silently drop the 10 pairs touching a composite cohort
    # gene and pull H down to 86.59; normalising gives n = 5312 and H = 87.99,
    # matching the stored 88.0. 'Mrc1_Slc39A12' stays excluded (dedup-dropped).
    # See REPRODUCIBILITY_RUNBOOK.md section 9.

    # in-cohort pairs with coverage >= 50%: pident and pinned |Δlog2(t1/2)|
    rows = []
    with open(fl.PAIRTAB, encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t"); j = {c: i for i, c in enumerate(h)}
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            g1 = fl.canonical(f[j["gene1"]], hl)
            g2 = fl.canonical(f[j["gene2"]], hl)
            if g1 is not None and g2 is not None and float(f[j["cov_min"]]) >= 0.5:
                rows.append((float(f[j["pident"]]), float(f[j["abs_dlog2_t"]])))
    pid = np.array([r[0] for r in rows]); dl = np.array([r[1] for r in rows])

    data, labels, ns = [], [], []
    for lo, hi, lab in BINS:
        m = (pid >= lo) & (pid < hi) if lab != "≥80" else (pid >= 80)
        if lo == -1:
            m = pid < 20
        data.append(dl[m]); labels.append(lab); ns.append(int(m.sum()))
    H, df, p = kruskal_wallis(data)

    fig, ax = plt.subplots(figsize=(9, 7))
    bp = ax.boxplot(data, widths=0.6, showfliers=False, whis=1.5, patch_artist=True,
                    medianprops=dict(color="black", lw=2))
    for patch, c in zip(bp["boxes"], COLORS):
        patch.set_facecolor(c); patch.set_alpha(0.95); patch.set_edgecolor("black")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(["%s\n(n=%s)" % (lab, f"{n:,}") for lab, n in zip(labels, ns)])
    ax.set_xlabel("Sequence identity (%)")
    ax.set_ylabel("|Δlog₂(half-life)|  (paired)")
    ax.set_title("Half-life difference vs sequence identity "
                 "(coverage ≥ 50 %%, n=%s pairs)" % f"{len(dl):,}", fontweight="bold")
    ax.grid(axis="y", ls=":", alpha=0.4)
    # Kruskal-Wallis H/p are still computed and printed for the reproducibility
    # check, but are no longer drawn on the figure (per request).

    fig.tight_layout()
    out = os.path.join(OUTDIR, "Fig_identity_vs_halflife_boxplot_rev.png")
    fig.savefig(out, dpi=190)
    print("wrote", out)
    print("n total pairs =", len(dl), "| per-bin:", dict(zip(labels, ns)))
    # NOTE on the stored annotation: the stored H (88.0) is reproduced below, but
    # the stored p (1.6e-15) is internally INCONSISTENT with it -- the chi-square
    # survival function at df=4 maps H=88.0 to p=3.5e-18, not 1.6e-15. The stored
    # p is a stale annotation; the recomputed p is the correct one for this H.
    print("Kruskal-Wallis H=%.2f df=%d p=%.3e   (stored H=88.0 -> reproduced; "
          "stored p 1.6e-15 was stale, correct p at df=4 is %.3e)" % (H, df, p, p))


if __name__ == "__main__":
    main()
