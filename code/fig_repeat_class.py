"""
Regenerate Fig_repeat_class_halflife_rev.png from the frozen snapshot, with the
y-axis and in-figure annotations in DAYS (values unchanged; only the unit label
differs from the original 'h' rendering).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
SNAP = os.path.join(PACK, "data", "snapshots", "Repeat_class_halflife_rev.tsv")
OUTDIR = os.environ.get("FIG_OUTDIR", "/mnt/workspace/output")

ORDER = ["OTHER", "TPR", "ANK", "WD40", "LRR", "HEAT_ARM", "KELCH_FBX", "SPECTRIN"]
DISP = {"OTHER": "Other", "TPR": "TPR", "ANK": "ANK", "WD40": "WD40", "LRR": "LRR",
        "HEAT_ARM": "HEAT/ARM", "KELCH_FBX": "Kelch/F-box", "SPECTRIN": "Spectrin"}
COL = {"OTHER": "#9e9e9e", "TPR": "tab:blue", "ANK": "tab:orange", "WD40": "tab:green",
       "LRR": "tab:red", "HEAT_ARM": "tab:purple", "KELCH_FBX": "tab:brown",
       "SPECTRIN": "tab:pink"}


def load():
    g, hl, cls = [], [], []
    with open(SNAP, encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        j = {c: i for i, c in enumerate(h)}
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            g.append(f[j["alias"]]); hl.append(float(f[j["half_life_h"]])); cls.append(f[j["repeat_class"]])
    return np.array(g), np.array(hl), np.array(cls)


def main():
    g, hl, cls = load()
    other_med = float(np.median(hl[cls == "OTHER"]))
    allrep = hl[cls != "OTHER"]
    thr = float(np.percentile(hl, 90))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 5.6),
                                   gridspec_kw=dict(width_ratios=[1, 1.25]))

    # ---- Panel A: boxplot by class ----
    data = [hl[cls == c] for c in ORDER]
    bp = axL.boxplot(data, widths=0.62, showfliers=False, whis=1.5,
                     patch_artist=True, medianprops=dict(color="black", lw=1.6))
    for patch, c in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COL[c]); patch.set_alpha(0.85); patch.set_edgecolor("black")
    axL.set_yscale("log")
    axL.set_ylabel("Protein half-life (days, log scale)")
    axL.set_xlabel("Repeat-domain class")
    axL.set_xticks(range(1, len(ORDER) + 1))
    axL.set_xticklabels([f"{DISP[c]} (n={int((cls==c).sum()):,})" for c in ORDER],
                        rotation=30, ha="right", fontsize=9)

    # ---- Panel B: long-lived tail ----
    tail_mask = (hl >= thr) & (cls != "OTHER")
    pos = {c: i for i, c in enumerate(ORDER)}
    tg, th, tc = g[tail_mask], hl[tail_mask], cls[tail_mask]
    placed = []                          # (x, log10y) of labels already placed
    for c in ORDER:                      # spread same-class points horizontally
        idx = [k for k in range(len(tg)) if tc[k] == c]
        m = len(idx)
        for rank, k in enumerate(sorted(idx, key=lambda m: -th[m])):
            off = (rank - (m - 1) / 2.0) * 0.30
            x = pos[c] + off
            axR.scatter([x], [th[k]], s=46, color=COL[c], edgecolor="black", lw=0.5, zorder=3)
            # choose label side/offset that avoids already-placed labels
            ly = np.log10(th[k])
            dy = 0.0
            for (px, ply) in placed:
                if abs(px - x) < 0.5 and abs((ly + dy) - ply) < 0.05:
                    dy += 0.06
            placed.append((x, ly + dy))
            axR.annotate(tg[k], (x, th[k]),
                         xytext=(0, 9 + dy * 90), textcoords="offset points",
                         fontsize=7.5, ha="center", va="bottom",
                         arrowprops=dict(arrowstyle="-", lw=0.4, color="0.5") if dy else None)
    axR.set_yscale("log")
    axR.set_xlim(-0.7, len(ORDER) - 0.3)
    axR.set_ylim(top=th.max() * 2.0)
    axR.set_xticks(range(len(ORDER)))
    axR.set_xticklabels([DISP[c] for c in ORDER], rotation=30, ha="right", fontsize=9)
    axR.set_ylabel("Half-life (days, log scale)")
    axR.set_title(f"Repeat proteins in the long-lived tail (n = {int(tail_mask.sum())})")

    fig.tight_layout()
    out = os.path.join(OUTDIR, "Fig_repeat_class_halflife_rev.png")
    fig.savefig(out, dpi=200)
    print("wrote", out, "| Other median %.2f days, tail thr %.1f days, tail n=%d"
          % (other_med, thr, int(tail_mask.sum())))


if __name__ == "__main__":
    main()
