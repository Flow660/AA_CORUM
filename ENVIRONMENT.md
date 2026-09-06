# Environment manifest — Cyno proteome half-life pipeline

This pack was validated in the Copilot Cowork sandbox with the following runtime.
Reproduction is robust to minor version drift, but record this for provenance.

| Component | Version used | Notes |
|---|---|---|
| Python | 3.12.9 | CPython |
| NumPy | 2.5.0 | the **only** required third-party package |
| matplotlib | 3.11.0 | figures only; not needed for the statistics |
| pandas | 3.0.3 | present but **not used** by the pipeline |
| SciPy | **absent** | intentionally — all stats are NumPy-only re-implementations |

## Hard requirements
- **NumPy only** for every statistical result. `code/stats_numpy.py` re-implements
  ranks, Spearman, Mann-Whitney (tie-corrected), Cliff's delta, OLS, and the
  Student-t p-value (via a regularized incomplete-beta evaluation) so that **no
  SciPy is required**. Do not "upgrade" to SciPy — results are pinned to these
  implementations.
- **Determinism**: every permutation uses `numpy.default_rng(42)`. Seed = 42.
- **No network access is needed.** All external data (UniProt features, domain
  calls, GO terms, CORUM membership) is captured in `data/snapshots/`. The
  pipeline never re-fetches a live database, so results do not drift with the
  state of UniProt / CORUM / GOA on the day you run it.

## Tools that are NOT available in the sandbox (and are not needed here)
- **DIAMOND** (the BLAST-like aligner used to build `cyno_pairs.tsv`). It is a
  compiled binary and cannot be installed. The DIAMOND output is therefore pinned
  as `data/inputs/cyno_pairs.tsv`; for a brand-new species you must produce the
  pairs file outside the sandbox and upload it.
- `pip install` / `apt` — package installation is blocked. Work with NumPy only.

## Permutation runtime
Tier 1-5's within-complex SD permutation is the only slow step. With the
canonical `N_PERM = 10000` and the whole-cohort pool it takes a few minutes in
the sandbox. For a quick sanity check use `N_PERM=2000` (≈45 s); the permutation
p has Monte-Carlo noise of order ±0.01 at that size, which does not change any
conclusion.
