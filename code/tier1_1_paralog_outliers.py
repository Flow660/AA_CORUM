"""
Tier 1-1  Outlier paralog pairs.

Finds paralog pairs (high sequence identity) whose half-lives differ strongly.

REPRODUCIBILITY NOTE
--------------------
This tier is reproduced from the FROZEN evaluated-pairs table
(pair_table_rev.tsv: gene1, gene2, pident, cov_min, abs_dlog2_t), NOT by
re-running DIAMOND. Two details are required to land on the published 31-pair set
and are the reason a naive re-run diverges:

  1. pair_table_rev.tsv was built on the PRE-dedup cohort (3,515 rows) and still
     contains the ambiguous `Mrc1;Slc39A12` composite as `Mrc1_Slc39A12`. The
     published rev set drops it, so we re-apply the rev-cohort gene filter
     (cyno_target_rev.tsv, 3,514 genes). That removes exactly the stale
     Mrc1 / Mrc1_Slc39A12 self-pair (pident 100, cov 1.0).
  2. The published "id80 / 2-fold" set applies pident >= 80, |dlog2(t1/2)| >= 1
     (i.e. >= 2-fold), AND an alignment-coverage QC of cov_min >= 0.5 — the
     coverage cut is not stated in the v1 prose but is recovered here.

With the saved table + these two rules, the result is EXACTLY the 31 published
pairs (0 extra, 0 missing).

Inputs : data/snapshots/pair_table_rev.tsv   (evaluated pairs)
         data/inputs/cyno_target_rev.tsv      (rev cohort: gene -> half-life)
Output : results/Outlier_pairs_id80_2fold_rev.tsv
Expected: data/expected/Outlier_pairs_id80_2fold_rev.tsv  (31 pairs)
"""
import os
import numpy as np

PIDENT_MIN = 80.0
DLOG2_MIN = 1.0      # 2-fold
COV_MIN = 0.5

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
PAIRTAB = os.path.join(PACK, "data", "snapshots", "pair_table_rev.tsv")
COHORT = os.path.join(PACK, "data", "inputs", "cyno_target_rev.tsv")
OUT = os.path.join(PACK, "results", "Outlier_pairs_id80_2fold_rev.tsv")


def load_cohort():
    by_gene = {}
    with open(COHORT, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            by_gene[f[idx["row_gene"]]] = (f[idx["accession"]],
                                           float(f[idx["half_life_h"]]))
    # Accept both spellings of composite gene names ('Rack1;Gnb2L1' vs
    # 'Rack1_Gnb2L1'); see family_lib.add_composite_aliases. The dedup-dropped
    # 'Mrc1;Slc39A12' is still excluded, so the stale self-pair remains filtered.
    from family_lib import add_composite_aliases
    add_composite_aliases(by_gene)
    return by_gene


def main():
    cohort = load_cohort()
    rows = []
    with open(PAIRTAB, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            g1, g2 = f[idx["gene1"]], f[idx["gene2"]]
            pident = float(f[idx["pident"]])
            cov_min = float(f[idx["cov_min"]])
            dlog2 = float(f[idx["abs_dlog2_t"]])
            # rev-cohort filter drops the stale Mrc1;Slc39A12 composite pair
            if g1 not in cohort or g2 not in cohort:
                continue
            if pident < PIDENT_MIN or dlog2 < DLOG2_MIN or cov_min < COV_MIN:
                continue
            a1, h1 = cohort[g1]
            a2, h2 = cohort[g2]
            if h1 >= h2:
                gl, gs, al, ash, hl_l, hl_s = g1, g2, a1, a2, h1, h2
            else:
                gl, gs, al, ash, hl_l, hl_s = g2, g1, a2, a1, h2, h1
            fold = hl_l / hl_s
            rows.append((gl, gs, al, ash, pident, cov_min, hl_l, hl_s, fold,
                         abs(np.log2(fold))))

    rows.sort(key=lambda r: r[9], reverse=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as w:
        w.write("gene_longer_t\tgene_shorter_t\tacc_longer\tacc_shorter\tpident"
                "\tcov_min\thalf_life_longer_h\thalf_life_shorter_h"
                "\tfold_difference\tabs_dlog2\n")
        for r in rows:
            w.write("\t".join(str(x) for x in r) + "\n")

    print("=== Tier 1-1 outlier paralog pairs (recomputed) ===")
    print(f"filter: rev-cohort genes, pident >= {PIDENT_MIN}, "
          f"|dlog2| >= {DLOG2_MIN} (2-fold), cov_min >= {COV_MIN}")
    print(f"outlier pairs = {len(rows)}   (expected 31)")
    print("top 5 by |dlog2|:")
    for r in rows[:5]:
        print(f"  {r[0]:>9} / {r[1]:<9} pident={r[4]:.1f} cov={r[5]:.3f} "
              f"fold={r[8]:.2f} dlog2={r[9]:.3f}")
    print(f"\nwrote {os.path.relpath(OUT, PACK)}")


if __name__ == "__main__":
    main()
