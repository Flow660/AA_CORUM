"""
Run the full Cyno half-life pipeline from the pinned snapshots and verify each
tier against the expected values.

Usage:
    python3 run_all.py            # canonical (N_PERM=10000 for Tier 1-5)
    N_PERM=2000 python3 run_all.py  # faster permutation for a quick check

Every tier reads only files under ../data/ (inputs + frozen snapshots), so this
runs with no network access and no live UniProt/CORUM/GOA fetch.
"""
import runpy
import os
import sys

TIERS = [
    ("Tier 1-1  paralog outliers   ", "tier1_1_paralog_outliers.py"),
    ("Tier 1-1b family dispersion  ", "tier1_1b_family_dispersion_sweep.py"),
    ("Tier 1-2  repeat-class        ", "tier1_2_repeat_class.py"),
    ("Tier 1-4  sequence regression ", "tier1_4_seq_regression.py"),
    ("Tier 1-5  CORUM membership    ", "tier1_5_corum_membership.py"),
    ("Tier 1-6  STOICH vs REGUL     ", "tier1_6_stoich_regul.py"),
]

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    for label, script in TIERS:
        print("\n" + "=" * 72)
        print(label)
        print("=" * 72)
        runpy.run_path(os.path.join(HERE, script), run_name="__main__")
    print("\nAll tiers complete. Compare results/ against data/expected/ and the")
    print("acceptance table in MANIFEST.md.")
    print("\nTier 1-3 (GO category comparison) is NOT included: its outputs were "
          "never\nsaved and it depends on live UniProt-GOA. See REPRODUCIBILITY_RUNBOOK.md.")


if __name__ == "__main__":
    main()
