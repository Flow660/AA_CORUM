"""
Tier 1-4  Sequence-feature regression  (mechanistic-predictability ceiling).

Re-derives the OLS multiple regression of half-life on linear sequence
features, consuming the FROZEN feature snapshot so no live UniProt fetch is
needed. Reproduces the published null finding: whole-model R^2 = 0.0190, with
pI the only strongly significant feature.

Inputs  : data/snapshots/seq_features_rev.tsv   (gene, ..., pI, disorder_frac,
                                                  nend_class, pest_max, ...)
Outputs : results/seq_features_summary_rev.tsv
Expected: data/expected/Supp_Table_seq_features_rev.tsv
"""
import os
import numpy as np
from stats_numpy import spearman, ols

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
SNAP = os.path.join(PACK, "data", "snapshots", "seq_features_rev.tsv")
OUT = os.path.join(PACK, "results", "seq_features_summary_rev.tsv")
# Supplementary Table S2 -- the per-feature publication table. Regenerated here
# (added 2026-09-06); previously it existed only as a pinned acceptance target
# under data/expected/ and no script produced it.
OUT_SUPP = os.path.join(PACK, "results", "Supp_Table_seq_features_rev.tsv")

NEND_TYPE1 = {"R", "K", "H"}
NEND_TYPE2 = {"L", "F", "W", "Y", "I"}
NEND_NDEG = {"D", "E", "N", "Q"}


def load():
    rows = []
    with open(SNAP, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rows.append(f)
    col = lambda f, name: f[idx[name]]
    gene = [col(f, "gene") for f in rows]
    hl = np.array([float(col(f, "half_life_h")) for f in rows])
    length = np.array([float(col(f, "length")) for f in rows])
    pI = np.array([float(col(f, "pI")) for f in rows])
    disorder = np.array([float(col(f, "disorder_frac")) for f in rows])
    pest = np.array([float(col(f, "pest_max")) for f in rows])
    nend_cls = [col(f, "nend_class") for f in rows]
    return gene, hl, length, pI, disorder, pest, nend_cls


def main():
    gene, hl, length, pI, disorder, pest, nend_cls = load()
    n = len(hl)
    y = np.log10(hl)
    log_len = np.log10(length)

    # N-end one-hot dummies; STAB is the reference level.
    t1 = np.array([1.0 if c == "TYPE1" else 0.0 for c in nend_cls])
    t2 = np.array([1.0 if c == "TYPE2" else 0.0 for c in nend_cls])
    nd = np.array([1.0 if c == "NDEG" else 0.0 for c in nend_cls])

    # Continuous predictors are z-scored (standardized betas, as published);
    # binary N-end dummies are left on the raw 0/1 scale. Standardization
    # leaves R^2 and all p-values unchanged.
    def z(v):
        return (v - v.mean()) / v.std(ddof=0)

    X = np.column_stack([z(log_len), z(pI), z(disorder), z(pest), t1, t2, nd])
    fit = ols(X, y)

    names = ["intercept", "log_len", "pI", "disorder", "pest_max",
             "nend_TYPE1", "nend_TYPE2", "nend_NDEG"]

    # Spearman for the four continuous features
    sp = {}
    for nm, v in [("length (log10)", log_len), ("pI", pI),
                  ("disorder_frac", disorder), ("PEST_max", pest)]:
        sp[nm] = spearman(v, y)

    # Spearman for the three binary N-end dummies (Table S2 reports these too).
    sp_bin = {"nend_TYPE1": spearman(t1, y),
              "nend_TYPE2": spearman(t2, y),
              "nend_NDEG": spearman(nd, y)}

    # Median half-life per N-end class, in HOURS (published as `nend_median`).
    nend_med = []
    for cls_name in ["STAB", "TYPE2", "TYPE1", "NDEG"]:
        mask = np.array([c == cls_name for c in nend_cls])
        if mask.any():
            nend_med.append((cls_name, int(mask.sum()), float(np.median(hl[mask]))))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as w:
        w.write("analysis\tlabel\tn\tstat\tp\n")
        for nm in ["length (log10)", "pI", "disorder_frac", "PEST_max"]:
            rho, p = sp[nm]
            w.write(f"spearman\t{nm}\t{n}\t{rho:.4f}\t{p:.4e}\n")
        for cls_name, k, med in nend_med:
            w.write(f"nend_median\t{cls_name}\t{k}\t{med:.4f}\tNA\n")
        for nm, b, p in zip(names, fit["beta"], fit["p"]):
            tag = "intercept" if nm == "intercept" else nm
            w.write(f"OLS_beta\t{tag}\t{n}\t{b:.4f}\t{p:.4e}\n")
        w.write(f"OLS_R2\twhole_model\t{n}\t{fit['r2']:.4f}\tNA\n")

    # ---- Supplementary Table S2 (one row per feature)
    SUPP_ROWS = [
        ("log10(length)", "continuous", "length (log10)", "log_len"),
        ("pI", "continuous", "pI", "pI"),
        ("Disorder fraction", "continuous", "disorder_frac", "disorder"),
        ("PEST max", "continuous", "PEST_max", "pest_max"),
        ("N-end TYPE1 (R/K/H)", "binary", "nend_TYPE1", "nend_TYPE1"),
        ("N-end TYPE2 (L/F/W/Y/I)", "binary", "nend_TYPE2", "nend_TYPE2"),
        ("N-end NDEG (D/E/N/Q)", "binary", "nend_NDEG", "nend_NDEG"),
    ]
    by_name = {nm: (b, se, pv) for nm, b, se, pv
               in zip(names, fit["beta"], fit["se"], fit["p"])}
    with open(OUT_SUPP, "w", encoding="utf-8") as w:
        w.write("feature\ttype\tn\tspearman_rho\tspearman_p\t"
                "OLS_beta\tOLS_SE\tOLS_p\n")
        for label, kind, sp_key, ols_key in SUPP_ROWS:
            rho, sp_p = sp[sp_key] if kind == "continuous" else sp_bin[sp_key]
            b, se, pv = by_name[ols_key]
            w.write(f"{label}\t{kind}\t{n}\t{rho:+.4f}\t{sp_p:.3e}\t"
                    f"{b:+.4f}\t{se:.4f}\t{pv:.3e}\n")
        w.write(f"WHOLE_MODEL\tOLS\t{n}\tNA\tNA\tR2={fit['r2']:.4f}\tNA\tNA\n")

    print("=== Tier 1-4 sequence-feature regression (recomputed) ===")
    print(f"n = {n}")
    print(f"whole-model R^2 = {fit['r2']:.4f}   (expected 0.0190)")
    print(f"adjusted R^2    = {fit['adj_r2']:.4f}")
    print(f"{'feature':<14}{'beta':>10}{'SE':>9}{'p':>13}")
    for nm, b, s, p in zip(names, fit["beta"], fit["se"], fit["p"]):
        print(f"{nm:<14}{b:>10.4f}{s:>9.4f}{p:>13.3e}")
    print("\nSpearman (continuous):")
    for nm in ["length (log10)", "pI", "disorder_frac", "PEST_max"]:
        rho, p = sp[nm]
        print(f"  {nm:<16} rho={rho:+.4f}  p={p:.3e}")
    print("\nMedian half-life by N-end class (hours):")
    for cls_name, k, med in nend_med:
        print(f"  {cls_name:<6} n={k:<5} median={med:.4f}")
    print(f"\nwrote {os.path.relpath(OUT, PACK)}")
    print(f"wrote {os.path.relpath(OUT_SUPP, PACK)}")


if __name__ == "__main__":
    main()
