"""
Shared helpers for the paralog within-family dispersion analysis (Tier 1-1b and
its two figures), so the table and the figures use ONE seeded engine.

Method (recovered from the original outputs):
  * families = connected components of the curated pair graph (pair_table_rev),
    at identity >= thr, mutual coverage >= 50%, family size >= 3;
  * within-family dispersion = POPULATION SD (ddof=0) of log10(half-life);
  * null = SD of 20,000 size-matched random sets drawn (without replacement)
    from the full proteome;
  * one-sided p = bootstrap of the median over `n_families` resamples of the
    null pool (5,000 replicates), with floor 1/(n_boot+1).
All randomness uses numpy.default_rng(seed) with seed=42.
"""
import os
import numpy as np

SEED = 42
COV_MIN = 0.5
MIN_SIZE = 3
N_SETS = 20000
N_BOOT = 5000

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
PAIRTAB = os.path.join(PACK, "data", "snapshots", "pair_table_rev.tsv")
COHORT = os.path.join(PACK, "data", "inputs", "cyno_target_rev.tsv")


def canonical(g, hl):
    """Map a pair-table gene name onto its cohort spelling, or None if absent.

    The pair table joins composite gene names with '_' ('Rack1_Gnb2L1') while the
    cohort tsv uses ';' ('Rack1;Gnb2L1'). Resolving to the COHORT spelling (rather
    than adding alias keys to the cohort dict) keeps the half-life pool at exactly
    one entry per protein -- see the warning in add_composite_aliases().

    'Mrc1_Slc39A12' resolves to None: it is the row removed by the build_cohort
    dedup rule and is absent from the cohort under either spelling, so the stale
    self-pair stays filtered.
    """
    if g in hl:
        return g
    alt = g.replace("_", ";")
    return alt if alt in hl else None


def add_composite_aliases(d):
    """Accept BOTH spellings of a composite gene name (adopted 2026-08-31).

    ⚠️ MEMBERSHIP LOOKUP ONLY. This ADDS alias keys, so `len(d)` and `d.values()`
    will double-count the composite genes. Never build a sampling pool or a count
    from a dict that has been through this function -- use the canonical cohort
    (or `canonical()` above) for anything statistical. Tier 1-1 uses it safely
    because it only ever does per-pair lookups.

    The pair table joins composite gene names with '_' ('Rack1_Gnb2L1') while the
    cohort tsv uses ';' ('Rack1;Gnb2L1'). Without this, a plain `gene in cohort`
    membership test silently drops every pair touching a composite cohort gene.

    Four cohort genes are affected: H2Ac12;..., H2Bc17;..., H2Bc5;..., Rack1;Gnb2L1.

    Genes ABSENT from the cohort are never reintroduced -- in particular
    'Mrc1_Slc39A12', the row removed by the build_cohort dedup rule, stays out
    because 'Mrc1;Slc39A12' is not a cohort key under either spelling.

    This is the scientifically correct behaviour and is now applied consistently
    across Tier 1-1, Tier 1-1b and both paralog figures. It supersedes the
    original naive lookup; see REPRODUCIBILITY_RUNBOOK.md section 9.
    """
    for g in [k for k in d if ";" in k]:
        d[g.replace(";", "_")] = d[g]
    return d


def load():
    hl = {}
    with open(COHORT, encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t"); j = {c: i for i, c in enumerate(h)}
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t"); hl[f[j["row_gene"]]] = float(f[j["half_life_h"]])
    # NOTE: `hl` is returned CANONICAL (one entry per protein, n = 3,514) so that
    # `hl.values()` is a clean sampling pool. Composite-name normalisation happens
    # inside build_families() via canonical().
    pairs = []
    with open(PAIRTAB, encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t"); j = {c: i for i, c in enumerate(h)}
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            pairs.append((f[j["gene1"]], f[j["gene2"]], float(f[j["pident"]]), float(f[j["cov_min"]])))
    return hl, pairs


def build_families(hl, pairs, thr):
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for g1, g2, pid, cv in pairs:
        if pid < thr or cv < COV_MIN:
            continue
        c1, c2 = canonical(g1, hl), canonical(g2, hl)
        if c1 is not None and c2 is not None:
            a, b = find(c1), find(c2)
            if a != b:
                parent[a] = b
    comp = {}
    for g in list(parent):
        comp.setdefault(find(g), []).append(g)
    return [m for m in comp.values() if len(m) >= MIN_SIZE]


def within_family_sds(families, hl):
    return np.array([np.log10([hl[g] for g in f]).std(ddof=0) for f in families])


def random_set_sds(sizes, pool, n_sets=N_SETS, seed=SEED):
    """SD (ddof=0) of n_sets size-matched random draws (without replacement)
    from `pool` (log10 half-lives). Sizes are cycled to match the family-size
    distribution exactly."""
    sizes_full = np.resize(np.asarray(sizes), n_sets)
    N = len(pool)
    out = np.empty(n_sets)
    rng = np.random.default_rng(seed)
    order = np.argsort(sizes_full, kind="mergesort")
    for k in np.unique(sizes_full):
        idx = order[sizes_full[order] == k]
        m = len(idx)
        R = rng.random((m, N))
        part = np.argpartition(R, k - 1, axis=1)[:, :k]
        out[idx] = pool[part].std(axis=1, ddof=0)
    return out


def bootstrap_p(within_med, pool_sds, n_fam, n_boot=N_BOOT, seed=SEED + 1):
    rng = np.random.default_rng(seed)
    meds = np.median(rng.choice(pool_sds, size=(n_boot, n_fam), replace=True), axis=1)
    return (int((meds <= within_med).sum()) + 1) / (n_boot + 1)


def analyse(thr, hl, pairs, pool):
    fams = build_families(hl, pairs, thr)
    sizes = np.array([len(f) for f in fams])
    within = within_family_sds(fams, hl)
    pool_sds = random_set_sds(sizes, pool, seed=SEED + int(thr))
    fam_med = float(np.median(within))
    rand_med = float(np.median(pool_sds))
    p = bootstrap_p(fam_med, pool_sds, len(fams), seed=SEED + 1000 + int(thr))
    return dict(families=len(fams), proteins=int(sizes.sum()), within=within,
                pool_sds=pool_sds, fam_med=fam_med, rand_med=rand_med, p=p)
