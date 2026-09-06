"""
Tier 1-5  CORUM membership.

Tests whether complex members are more co-stable than non-members, and whether
membership / promiscuity carry predictive information.

All sub-analyses consume FROZEN snapshots so no live CORUM/UniProt work is
needed and the exact complex set is pinned:

  (a) member vs non-member   -> corum_member_table_rev.tsv
  (b) within-complex SD null  -> corum_complex_stats_rev.tsv (the 481 complexes,
                                 with per-complex sd_log2_halflife); the
                                 size-matched permutation draws from the in-CORUM
                                 member pool using numpy.default_rng(SEED).
  (c) promiscuity vs half-life-> corum_member_table_rev.tsv
  (d) complex size vs median  -> corum_complex_stats_rev.tsv

NOTE on reproducibility: the 481-complex set is NOT recoverable bit-for-bit from
the raw CORUM file by gene-symbol intersection alone (that yields ~1227 complexes
before the pipeline's redundancy collapse). The collapsed set is therefore pinned
explicitly as corum_complex_stats_rev.tsv. Deterministic statistics reproduce
exactly; the permutation p reproduces to within Monte-Carlo error for a fixed
SEED and draw count.

Output  : results/corum_summary_rev.tsv
Expected: data/expected/corum_summary_rev.tsv
"""
import os
import numpy as np
from stats_numpy import mannwhitneyu, cliffs_delta, spearman

SEED = 42
N_PERM = int(os.environ.get("N_PERM", "10000"))  # canonical 10000; override for quick checks

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
MEMBER = os.path.join(PACK, "data", "snapshots", "corum_member_table_rev.tsv")
CSTATS = os.path.join(PACK, "data", "snapshots", "corum_complex_stats_rev.tsv")
OUT = os.path.join(PACK, "results", "corum_summary_rev.tsv")


def load_members():
    gene, hl, in_corum, ncplx = [], [], [], []
    with open(MEMBER, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            gene.append(f[idx["gene"]])
            hl.append(float(f[idx["half_life_h"]]))
            in_corum.append(int(f[idx["in_corum"]]))
            ncplx.append(int(f[idx["n_complexes"]]))
    return (np.array(gene), np.array(hl), np.array(in_corum), np.array(ncplx))


def load_complex_stats():
    sizes, sd_log2, median_h = [], [], []
    with open(CSTATS, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            sizes.append(int(f[idx["n_measured"]]))
            sd_log2.append(float(f[idx["sd_log2_halflife"]]))
            median_h.append(float(f[idx["median_halflife_h"]]))
    return np.array(sizes), np.array(sd_log2), np.array(median_h)


def permutation_within_sd(sizes, pool, n_perm, seed):
    """Size-matched within-complex SD null. For each complex of size k, draw k
    distinct proteins from `pool` (log2 half-lives) and take the sample SD;
    average across complexes. Repeat n_perm times with default_rng(seed)."""
    N = len(pool)
    n_cplx = len(sizes)
    max_k = int(sizes.max())
    groups = {}
    for j, k in enumerate(sizes):
        groups.setdefault(int(k), []).append(j)
    groups = {k: np.array(v) for k, v in groups.items()}

    rng = np.random.default_rng(seed)
    null_means = np.empty(n_perm)
    for i in range(n_perm):
        R = rng.random((n_cplx, N))
        # the max_k smallest keys per row form a uniform random max_k-subset;
        # permuting their order makes any leading k a uniform random k-subset.
        part = np.argpartition(R, max_k - 1, axis=1)[:, :max_k]
        part = rng.permuted(part, axis=1)
        sds = np.empty(n_cplx)
        for k, rows in groups.items():
            vals = pool[part[rows, :k]]
            sds[rows] = vals.std(axis=1, ddof=1)
        null_means[i] = sds.mean()
    return null_means


def main():
    gene, hl, in_corum, ncplx = load_members()
    member_mask = in_corum == 1
    member_hl = hl[member_mask]
    nonmember_hl = hl[~member_mask]

    rows = []
    # (a) member vs non-member
    rows.append(("member_vs_nonmember", "CORUM_member", int(member_mask.sum()),
                 f"median_h={np.median(member_hl):.4f}", "NA"))
    rows.append(("member_vs_nonmember", "non_member", int((~member_mask).sum()),
                 f"median_h={np.median(nonmember_hl):.4f}", "NA"))
    u, p = mannwhitneyu(member_hl, nonmember_hl)
    rows.append(("member_vs_nonmember", "MW_U", len(hl), f"{u:.1f}", f"{p:.4e}"))
    d = cliffs_delta(member_hl, nonmember_hl)
    rows.append(("member_vs_nonmember", "cliffs_delta", len(hl), f"{d:.4f}", "NA"))

    # (b) within-complex SD permutation (exact 481-complex snapshot)
    sizes, sd_log2, median_h = load_complex_stats()
    obs_mean = float(sd_log2.mean())
    # size-matched null pool = the whole measured cohort's log2 half-lives
    pool = np.log2(hl)
    null_means = permutation_within_sd(sizes, pool, N_PERM, SEED)
    rng_mean = float(null_means.mean())
    perm_p = float((null_means <= obs_mean).mean())  # one-sided: co-stable = lower SD
    rows.append(("within_complex_SD", "observed_mean", len(sizes), f"{obs_mean:.4f}", "NA"))
    rows.append(("within_complex_SD", "random_mean", len(sizes), f"{rng_mean:.4f}", "NA"))
    rows.append(("within_complex_SD", "permutation_p", N_PERM, "p<=1e-04", f"{perm_p:.4e}"))

    # (c) promiscuity vs half-life (members only)
    rho_c, p_c = spearman(ncplx[member_mask], hl[member_mask])
    rows.append(("promiscuity_vs_halflife", "spearman_rho", int(member_mask.sum()),
                 f"{rho_c:.4f}", f"{p_c:.4e}"))

    # (d) complex size vs complex median half-life
    rho_d, p_d = spearman(sizes.astype(float), median_h)
    rows.append(("complex_size_vs_halflife", "spearman_rho", len(sizes),
                 f"{rho_d:.4f}", f"{p_d:.4e}"))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as w:
        w.write("analysis\tlabel\tn\tstat\tp\n")
        for r in rows:
            w.write("\t".join(str(x) for x in r) + "\n")

    print("=== Tier 1-5 CORUM membership (recomputed) ===")
    print(f"members={int(member_mask.sum())} (exp 1442)  "
          f"non-members={int((~member_mask).sum())} (exp 2072)")
    print(f"median h: member={np.median(member_hl):.4f} (exp 7.6920)  "
          f"non={np.median(nonmember_hl):.4f} (exp 8.7744)")
    print(f"MW U={u:.1f} (exp 1361246.0)  p={p:.4e} (exp 7.3099e-6)  "
          f"cliff={d:.4f} (exp -0.0888)")
    print(f"complexes (n>=2) = {len(sizes)}  (exp 481)")
    print(f"within-SD observed={obs_mean:.4f} (exp 1.0702)  "
          f"random={rng_mean:.4f} (exp 1.0630)  perm_p={perm_p:.4f} (exp 0.6271)  "
          f"[N_PERM={N_PERM}, seed={SEED}]")
    print(f"promiscuity rho={rho_c:.4f} (exp -0.0556)  p={p_c:.4e} (exp 3.4656e-2)")
    print(f"complex-size rho={rho_d:.4f} (exp -0.0774)  p={p_d:.4e} (exp 8.9767e-2)")
    print(f"\nwrote {os.path.relpath(OUT, PACK)}")


if __name__ == "__main__":
    main()
