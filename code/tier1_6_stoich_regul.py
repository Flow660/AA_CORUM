"""
Tier 1-6  Stoichiometric vs Regulatory decomposition.

Tests whether textbook structural complexes (ribosome, proteasome, ...) show
co-stability. Consumes the FROZEN stratum assignment snapshot (NONE / STOICH /
REGUL / BOTH), so the ~50 complex-name regex patterns do not need to be
re-derived to reproduce the class membership.

Reproduces the published null finding: STOICH vs REGUL p = 0.92.

Input   : data/snapshots/corum_strat_assignment_rev.tsv
Output  : results/corum_strat_summary_rev.tsv
Expected: data/expected/corum_strat_summary_rev.tsv
"""
import os
import numpy as np
from stats_numpy import mannwhitneyu, cliffs_delta

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
SNAP = os.path.join(PACK, "data", "snapshots", "corum_strat_assignment_rev.tsv")
OUT = os.path.join(PACK, "results", "corum_strat_summary_rev.tsv")


def load():
    classes = {}
    with open(SNAP, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            cls = f[idx["corum_class"]]
            hl = float(f[idx["half_life_h"]])
            classes.setdefault(cls, []).append(hl)
    return {k: np.array(v) for k, v in classes.items()}


def main():
    cls = load()
    none = cls["NONE"]
    stoich = cls["STOICH"]
    regul = cls["REGUL"]
    both = cls["BOTH"]

    rows = []
    rows.append(("class_median", "NONE", len(none), "median_h",
                 f"{np.median(none):.4f}", "NA"))
    rows.append(("class_median", "STOICH", len(stoich), "median_h",
                 f"{np.median(stoich):.4f}", "NA"))
    rows.append(("class_median", "REGUL", len(regul), "median_h",
                 f"{np.median(regul):.4f}", "NA"))
    rows.append(("class_median", "BOTH", len(both), "median_h",
                 f"{np.median(both):.4f}", "NA"))

    for name, grp in [("STOICH_vs_NONE", stoich),
                      ("REGUL_vs_NONE", regul),
                      ("BOTH_vs_NONE", both)]:
        u, p = mannwhitneyu(grp, none)
        d = cliffs_delta(grp, none)
        rows.append(("mann_whitney", name, len(grp), "U", f"{u:.1f}", f"{p:.4e}"))
        rows.append(("cliffs_delta", name, len(grp), "delta", f"{d:.4f}", "NA"))

    # STOICH vs REGUL (headline comparison)
    u, p = mannwhitneyu(stoich, regul)
    d = cliffs_delta(stoich, regul)
    rows.append(("mann_whitney", "STOICH_vs_REGUL_vs_NONE", len(stoich) + len(regul),
                 "U", f"{u:.1f}", f"{p:.4e}"))
    rows.append(("cliffs_delta", "STOICH_vs_REGUL_vs_NONE", len(stoich) + len(regul),
                 "delta", f"{d:.4f}", "NA"))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as w:
        w.write("analysis\tlabel\tn\tstat\tvalue\tp\n")
        for r in rows:
            w.write("\t".join(str(x) for x in r) + "\n")

    print("=== Tier 1-6 STOICH vs REGUL (recomputed) ===")
    print(f"class sizes: NONE={len(none)}  STOICH={len(stoich)}  "
          f"REGUL={len(regul)}  BOTH={len(both)}")
    print(f"medians (h): NONE={np.median(none):.4f}  STOICH={np.median(stoich):.4f}"
          f"  REGUL={np.median(regul):.4f}  BOTH={np.median(both):.4f}")
    u, p = mannwhitneyu(stoich, regul)
    print(f"STOICH vs REGUL: U={u:.1f}  p={p:.4e}   (expected U=71473.0  p=9.1616e-1)")
    u2, p2 = mannwhitneyu(both, none)
    print(f"BOTH vs NONE   : U={u2:.1f}  p={p2:.4e}   (expected U=139420.0 p=2.8169e-9)")
    print(f"\nwrote {os.path.relpath(OUT, PACK)}")


if __name__ == "__main__":
    main()
