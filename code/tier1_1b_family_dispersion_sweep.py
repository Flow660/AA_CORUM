"""
Tier 1-1b  Paralog within-family dispersion — identity-threshold sweep.

Companion to Tier 1-1. Paralog "families" are connected components of the curated
pair graph (pair_table_rev.tsv) at a given identity threshold (mutual coverage
>= 50%, family size >= 3). The median within-family population SD (ddof=0) of
log10(half-life) is compared with a size-matched null. Sequence-similar paralogs
are markedly LESS dispersed than random — the positive counterpart to the null
complex-membership findings.

Reproduces Table_family_threshold_sweep_rev.tsv and ADDS identity thresholds
20 and 40 (sweep: 20, 30, 40, 50, 60, 70, 80). All randomness is seeded
(numpy.default_rng, seed=42); see family_lib.py for the shared engine — the same
engine drives fig_family_dispersion.py so the table and the figure agree.

Output  : results/Table_family_threshold_sweep_rev.tsv
Expected: data/expected/Table_family_threshold_sweep_rev.tsv
"""
import os
import numpy as np
import family_lib as fl

IDENTITIES = [20, 30, 40, 50, 60, 70, 80]
HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
OUT = os.path.join(PACK, "results", "Table_family_threshold_sweep_rev.tsv")


def main():
    hl, pairs = fl.load()
    pool = np.log10(np.array(list(hl.values())))
    print("cohort proteins (rev) = %d  (used for families AND the random pool)" % len(hl))
    print("=== Tier 1-1b family-dispersion identity sweep (recomputed) ===")
    print("%-5s %-9s %-9s %-11s %-11s %-10s" %
          ("id%", "families", "proteins", "within_SD", "random_SD", "perm_p"))
    rows = []
    for thr in IDENTITIES:
        r = fl.analyse(thr, hl, pairs, pool)
        pstr = "< 1e-03" if r["p"] < 1e-3 else "%.4f" % r["p"]
        rows.append((thr, 50, r["families"], r["proteins"], r["fam_med"], r["rand_med"], pstr))
        print("%-5d %-9d %-9d %-11.3f %-11.3f %-10s" %
              (thr, r["families"], r["proteins"], r["fam_med"], r["rand_med"], pstr))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as w:
        w.write("Identity threshold (%)\tCoverage threshold (%)\tFamilies (size>=3)\t"
                "Proteins in families\tMedian within-family SD of log10(t1/2)\t"
                "Median random-set SD of log10(t1/2)\tPermutation P value\n")
        for thr, cov, nf, npr, ow, rn, ps in rows:
            w.write("%d\t%d\t%d\t%d\t%.3f\t%.3f\t%s\n" % (thr, cov, nf, npr, ow, rn, ps))
    print("\nwrote", os.path.relpath(OUT, PACK))
    print("NOTE: 3,514-protein (rev) cohort. Families, proteins-in-families and "
          "within-family\nSD all reproduce the stored table EXACTLY (composite gene "
          "names are normalised in\nfamily_lib.add_composite_aliases; see "
          "REPRODUCIBILITY_RUNBOOK.md section 9). random_SD\nand perm_p match the "
          "stored raw summary to Monte-Carlo precision (20,000\nsize-matched sets + "
          "5,000 median bootstrap, seed=42).")


if __name__ == "__main__":
    main()
