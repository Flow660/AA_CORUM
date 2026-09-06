"""
Regenerate Fig_seq_features_rev.png from the frozen snapshot, with every
half-life axis in DAYS (values unchanged; only the unit label differs from the
original 'h' rendering). Cohort n corrected to 3,514 in the title.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from stats_numpy import ols
SNAP = os.path.join(PACK, "data", "snapshots", "seq_features_rev.tsv")
OUTDIR = os.environ.get("FIG_OUTDIR", "/mnt/workspace/output")
YLAB = "Protein half-life (days, log)"


def load():
    rows = []
    with open(SNAP, encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        j = {c: i for i, c in enumerate(h)}
        for ln in fh:
            if not ln.strip():
                continue
            rows.append(ln.rstrip("\n").split("\t"))
    col = lambda n, t=float: np.array([t(r[j[n]]) for r in rows])
    nend = [r[j["nend_class"]] for r in rows]
    return (col("length"), col("pI"), col("disorder_frac"), col("pest_max"),
            col("half_life_h"), nend)


def deciles(x, y, nb=10):
    edges = np.quantile(x, np.linspace(0, 1, nb + 1))
    edges[-1] += 1e-9
    idx = np.digitize(x, edges[1:-1])
    cx, my, lo, hi = [], [], [], []
    for b in range(nb):
        m = idx == b
        if m.sum() < 3:
            continue
        cx.append(np.median(x[m])); my.append(np.median(y[m]))
        lo.append(np.percentile(y[m], 25)); hi.append(np.percentile(y[m], 75))
    cx, my, lo, hi = map(np.array, (cx, my, lo, hi))
    return cx, my, my - lo, hi - my


def scatter_panel(ax, x, y, rho, xlabel, title, logx=False):
    ax.scatter(x, y, s=4, color="0.7", alpha=0.35, lw=0, zorder=1)
    cx, my, el, eu = deciles(x, y)
    ax.errorbar(cx, my, yerr=[el, eu], fmt="s", color="tab:red", ms=5, lw=1.2,
                capsize=2, zorder=3)
    if logx:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel(YLAB); ax.set_title(title, loc="left")
    ax.text(0.03, 0.97, f"Spearman ρ = {rho:+.3f}\nn = {len(x)}",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))


def main():
    length, pI, dis, pest, hl, nend = load()
    n = len(hl)
    fig, ax = plt.subplots(2, 3, figsize=(15, 9))

    scatter_panel(ax[0, 0], length, hl, -0.063, "Protein length (aa, log)", "(a) Length", logx=True)
    scatter_panel(ax[0, 1], pI, hl, +0.108, "Isoelectric point (pI)", "(b) pI")
    scatter_panel(ax[0, 2], dis, hl, -0.051, "Disorder-promoting AA fraction", "(c) Intrinsic disorder")
    scatter_panel(ax[1, 0], pest, hl, -0.101, "PEST max score", "(d) PEST motif strength")

    # (e) N-end rule violins
    axe = ax[1, 1]
    classes = ["STAB", "TYPE2", "TYPE1", "NDEG"]
    labels = ["Stable\n(M/A/V/G/P/S/T/C)", "Type-2\n(L/F/W/Y/I)",
              "Type-1\n(R/K/H)", "Asn/Gln/\nAsp/Glu"]
    cols = ["#9e9e9e", "tab:orange", "tab:blue", "tab:red"]
    nend = np.array(nend)
    data = [np.log10(hl[nend == c]) for c in classes]
    parts = axe.violinplot(data, showextrema=False)
    for pc, c in zip(parts["bodies"], cols):
        pc.set_facecolor(c); pc.set_alpha(0.7)
    for i, d in enumerate(data, 1):
        axe.hlines(np.median(d), i - 0.3, i + 0.3, color="black", lw=1.6)
    axe.set_xticks(range(1, 5))
    axe.set_xticklabels([f"{l}\nn={int((nend==c).sum())}" for l, c in zip(labels, classes)], fontsize=7.5)
    axe.set_ylabel("Protein half-life (days, log10)")
    axe.set_title("(e) N-end rule class", loc="left")

    # (f) Multivariate OLS forest
    axf = ax[1, 2]
    y = np.log10(hl); ll = np.log10(length)
    z = lambda v: (v - v.mean()) / v.std(ddof=0)
    t1 = (nend == "TYPE1").astype(float); t2 = (nend == "TYPE2").astype(float); nd = (nend == "NDEG").astype(float)
    fit = ols(np.column_stack([z(ll), z(pI), z(dis), z(pest), t1, t2, nd]), y)
    names = ["log10(length)", "pI", "Disorder", "PEST max", "N-end TYPE1", "N-end TYPE2", "N-end NDEG"]
    beta, se, p = fit["beta"][1:], fit["se"][1:], fit["p"][1:]
    yp = np.arange(len(names))[::-1]
    for i, (b, s, pv, yy) in enumerate(zip(beta, se, p, yp)):
        c = "tab:red" if pv < 0.05 else "0.6"
        axf.errorbar([b], [yy], xerr=[[1.96 * s], [1.96 * s]], fmt="o", color=c, capsize=3)
        axf.text(b, yy + 0.18, ("p=%.0e" % pv) if pv < 0.05 else "p=%.2f" % pv,
                 ha="center", fontsize=7, color=c)
    axf.axvline(0, color="black", lw=0.8)
    axf.set_yticks(yp); axf.set_yticklabels(names)
    axf.set_xlabel("Standardised β (Δlog₁₀ t½ per +1 SD)")
    axf.set_title("(f) Multivariate OLS (R² = %.3f)" % fit["r2"], loc="left")

    fig.suptitle(f"Sequence features vs protein half-life — cyno proteome (n = {n})", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(OUTDIR, "Fig_seq_features_rev.png")
    fig.savefig(out, dpi=200)
    print("wrote", out, "| n=%d  R2=%.3f" % (n, fit["r2"]))


if __name__ == "__main__":
    main()
