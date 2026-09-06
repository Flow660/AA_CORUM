# RECOVERY — restore the Cyno pipeline after a Cowork refresh

A refresh wipes `input/` and `working/`. This pack lives in `output/` and syncs to
your OneDrive `Documents/Cowork`, so it survives. Four steps to reproduce:

1. **Re-upload the pack.** In the refreshed Cowork session, bring the whole
   `cyno_repro_pack/` folder back from OneDrive `Documents/Cowork` into the session.

2. **Run everything (no network needed):**
   ```bash
   cd cyno_repro_pack/code
   python3 run_all.py
   ```
   For a fast check (≈1 min) use `N_PERM=2000 python3 run_all.py`; the canonical
   run uses `N_PERM=10000` for the Tier 1-5 permutation (a few minutes).

3. **Verify.** Each tier prints `recomputed (expected …)`. Confirm against the
   acceptance table in `MANIFEST.md`. The headline checks:
   - Tier 1-4 whole-model **R² = 0.0190**
   - Tier 1-5 within-complex permutation **p ≈ 0.63**
   - Tier 1-6 STOICH vs REGUL **p = 0.92**
   - Tier 1-2 ALL_REPEAT **p = 1.984e-2**

4. **(Optional) check file integrity:**
   ```bash
   cd cyno_repro_pack
   sha256sum -c MANIFEST.sha256
   ```

That's it — no UniProt, CORUM, or GOA lookups are required; every external input
is frozen under `data/snapshots/`.

   - Tier 1-1 paralog outliers **31/31** (0 extra, 0 missing)

> Only Tier 1-3 (GO comparison) is not included — its outputs were never saved and
> it needs a live UniProt-GOA snapshot. See REPRODUCIBILITY_RUNBOOK.md §5.
