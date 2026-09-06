"""
Tier 1-2  Repeat-class half-life.

Tests whether repeat-domain proteins (TPR, ANK, WD40, LRR, HEAT/ARM, KELCH/FBX,
SPECTRIN) are longer-lived than non-repeat proteins. The per-protein repeat
class is derived from live UniProt feature tables in the original run; that
classification is FROZEN here as Repeat_class_halflife_rev.tsv so no live fetch
is needed.

Input   : data/snapshots/Repeat_class_halflife_rev.tsv  (..., repeat_class)
Output  : results/Repeat_class_summary_rev.tsv
Expected: data/expected/Repeat_class_summary_rev.tsv
"""
import os
import numpy as np
from stats_numpy import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
SNAP = os.path.join(PACK, "data", "snapshots", "Repeat_class_halflife_rev.tsv")
OUT = os.path.join(PACK, "results", "Repeat_class_summary_rev.tsv")

FAMILIES = ["TPR", "ANK", "WD40", "LRR", "HEAT_ARM", "KELCH_FBX", "SPECTRIN"]


def load():
    cls, hl = [], []
    with open(SNAP, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            cls.append(f[idx["repeat_class"]])
            hl.append(float(f[idx["half_life_h"]]))
    return np.array(cls), np.array(hl)


def pct(a, q):
    return float(np.percentile(a, q))


def main():
    cls, hl = load()
    other = hl[cls == "OTHER"]
    rows = []
    print("=== Tier 1-2 repeat-class half-life (recomputed) ===")
    print(f"{'class':<11}{'n':>4}{'median':>9}{'med_other':>11}{'U':>10}{'p':>12}")
    for fam in FAMILIES:
        grp = hl[cls == fam]
        if len(grp) == 0:
            continue
        u, p = mannwhitneyu(grp, other)
        u = min(u, len(grp) * len(other) - u)   # report min(U1,U2) as published
        rows.append((fam, len(grp), np.median(grp), pct(grp, 25), pct(grp, 75),
                     np.median(other), u, p))
        print(f"{fam:<11}{len(grp):>4}{np.median(grp):>9.3f}"
              f"{np.median(other):>11.3f}{u:>10.0f}{p:>12.3e}")

    allrep = hl[cls != "OTHER"]
    u, p = mannwhitneyu(allrep, other)
    u = min(u, len(allrep) * len(other) - u)
    rows.append(("ALL_REPEAT", len(allrep), np.median(allrep), pct(allrep, 25),
                 pct(allrep, 75), np.median(other), u, p))
    print(f"{'ALL_REPEAT':<11}{len(allrep):>4}{np.median(allrep):>9.3f}"
          f"{np.median(other):>11.3f}{u:>10.0f}{p:>12.3e}")
    print("expected ALL_REPEAT: n=148 median=6.807 med_other=8.414 U=220946 p=1.984e-2")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as w:
        w.write("class\tn\tmedian_h\tP25_h\tP75_h\tmedian_other_h\tU\tp_two_sided\n")
        for r in rows:
            fam, n, med, p25, p75, mo, u, p = r
            w.write(f"{fam}\t{n}\t{med:.3f}\t{p25:.3f}\t{p75:.3f}\t{mo:.3f}"
                    f"\t{u:.0f}\t{p:.3e}\n")
    print(f"\nwrote {os.path.relpath(OUT, PACK)}")


if __name__ == "__main__":
    main()
