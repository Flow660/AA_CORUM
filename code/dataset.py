"""
dataset.py — describes ONE species' worth of input files.

A `Dataset` is just a set of file paths plus the knowledge of which tiers each
path unlocks. The app uses it to run the pinned cyno reference or a user-supplied
species through exactly the same tier scripts, by overriding the path constants
the tier modules already expose.

WHAT A NEW SPECIES NEEDS
------------------------
Two files are produced OUTSIDE this repository and cannot be derived here:

  1. half-life workbook (.xlsx)  columns: gene name | (any) | half-life in HOURS
  2. protein FASTA               headers: >GENE|ACCESSION

From those, `build_cohort.py` produces the cohort tsv. A third external file,

  3. DIAMOND all-vs-all (.tsv)   outfmt 6: qseqid sseqid pident length qlen slen
                                 evalue bitscore

lets `build_pair_table.py` produce the evaluated pair table, which unlocks the
paralog tiers (1-1, 1-1b).

The remaining tiers each need one annotation table that must be derived from a
live UniProt / CORUM query for that species and then frozen:

  4. sequence features  -> Tier 1-4
  5. repeat classes     -> Tier 1-2
  6. CORUM member table + complex stats -> Tier 1-5
  7. CORUM stratification -> Tier 1-6

Nothing is fetched live. A tier whose file is missing is simply reported as
unavailable rather than silently producing a wrong answer.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
DATA = os.path.join(PACK, "data")

# key -> (label, required columns, which tiers it unlocks, how to obtain it)
SPEC = {
    "cohort": (
        "Cohort table",
        ["row_gene", "chosen_alias", "accession", "length", "half_life_h"],
        ["1-1", "1-1b"],
        "Built by build_cohort.py from the half-life workbook + FASTA.",
    ),
    "pairtab": (
        "Evaluated paralog pair table",
        ["gene1", "gene2", "pident", "cov_min", "abs_dlog2_t"],
        ["1-1", "1-1b"],
        "Built by build_pair_table.py from the raw DIAMOND all-vs-all output.",
    ),
    "seq_features": (
        "Sequence features",
        ["gene", "accession", "length", "half_life_h", "pI", "disorder_frac",
         "nend_class", "pest_max"],
        ["1-4"],
        "Derived from UniProt for your species, then frozen as a TSV.",
    ),
    "repeat_class": (
        "Repeat-domain classes",
        ["row_gene", "accession", "length", "half_life_h", "repeat_class"],
        ["1-2"],
        "Derived from UniProt domain annotations, then frozen as a TSV.",
    ),
    "corum_member": (
        "CORUM member table",
        ["gene", "accession", "half_life_h", "in_corum", "n_complexes"],
        ["1-5"],
        "Map your genes onto CORUM by UPPERCASE gene symbol, then freeze.",
    ),
    "corum_cstats": (
        "CORUM per-complex stats",
        ["complex_id", "complex_name", "n_measured", "median_halflife_h",
         "sd_log2_halflife", "measured_genes"],
        ["1-5"],
        "Collapse CORUM complexes to those with >=2 measured subunits.",
    ),
    "corum_strat": (
        "CORUM stratification",
        ["gene", "half_life_h", "n_complexes", "corum_class"],
        ["1-6"],
        "Assign NONE / STOICH / REGUL / BOTH from complex-name patterns.",
    ),
}

# tier -> the dataset keys it needs
TIER_NEEDS = {
    "1-1": ["cohort", "pairtab"],
    "1-1b": ["cohort", "pairtab"],
    "1-2": ["repeat_class"],
    "1-4": ["seq_features"],
    "1-5": ["corum_member", "corum_cstats"],
    "1-6": ["corum_strat"],
}

# where each key is injected: (module name, attribute)
INJECT = {
    "1-1": {"cohort": [("tier1_1_paralog_outliers", "COHORT")],
            "pairtab": [("tier1_1_paralog_outliers", "PAIRTAB")]},
    "1-1b": {"cohort": [("family_lib", "COHORT")],
             "pairtab": [("family_lib", "PAIRTAB")]},
    "1-2": {"repeat_class": [("tier1_2_repeat_class", "SNAP")]},
    "1-4": {"seq_features": [("tier1_4_seq_regression", "SNAP")]},
    "1-5": {"corum_member": [("tier1_5_corum_membership", "MEMBER")],
            "corum_cstats": [("tier1_5_corum_membership", "CSTATS")]},
    "1-6": {"corum_strat": [("tier1_6_stoich_regul", "SNAP")]},
}


class Dataset:
    def __init__(self, name, paths, reference=False):
        self.name = name
        self.paths = {k: v for k, v in paths.items() if v}
        self.reference = reference

    def has(self, key):
        p = self.paths.get(key)
        return bool(p) and os.path.exists(p)

    def missing_for(self, tier):
        return [k for k in TIER_NEEDS[tier] if not self.has(k)]

    def available_tiers(self):
        return [t for t in TIER_NEEDS if not self.missing_for(t)]

    def injections(self, tier):
        """[(module, attribute, path), ...] for this tier."""
        out = []
        for key, targets in INJECT[tier].items():
            path = self.paths.get(key)
            if not path:
                continue
            for mod, attr in targets:
                out.append((mod, attr, path))
        return out


def reference_dataset():
    """The pinned cynomolgus reference shipped with the repository."""
    return Dataset("Cynomolgus (pinned reference)", {
        "cohort": os.path.join(DATA, "inputs", "cyno_target_rev.tsv"),
        "pairtab": os.path.join(DATA, "snapshots", "pair_table_rev.tsv"),
        "seq_features": os.path.join(DATA, "snapshots", "seq_features_rev.tsv"),
        "repeat_class": os.path.join(DATA, "snapshots", "Repeat_class_halflife_rev.tsv"),
        "corum_member": os.path.join(DATA, "snapshots", "corum_member_table_rev.tsv"),
        "corum_cstats": os.path.join(DATA, "snapshots", "corum_complex_stats_rev.tsv"),
        "corum_strat": os.path.join(DATA, "snapshots", "corum_strat_assignment_rev.tsv"),
    }, reference=True)


def check_columns(path, key):
    """Return (ok, message) after checking the header of an uploaded file."""
    label, required, _, _ = SPEC[key]
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
    except Exception as exc:                                   # noqa: BLE001
        return False, f"{label}: could not read the file ({exc})"
    missing = [c for c in required if c not in header]
    if missing:
        return False, (f"{label}: missing required column(s) "
                       f"{', '.join(missing)}. Found: {', '.join(header[:10])}")
    return True, f"{label}: OK ({len(header)} columns)"


# --------------------------------------------------------------------- figures
# Each figure, the dataset keys it needs, and where those paths are injected.
# This is what makes the Figures page follow the ACTIVE dataset instead of only
# ever showing the pinned cynomolgus renderings.
FIGURES = {
    "seq_features": {
        "module": "fig_seq_features",
        "filename": "Fig_seq_features_rev.png",
        "title": "Tier 1-4 - sequence features",
        "needs": ["seq_features"],
        "inject": {"seq_features": [("fig_seq_features", "SNAP")]},
    },
    "repeat_class": {
        "module": "fig_repeat_class",
        "filename": "Fig_repeat_class_halflife_rev.png",
        "title": "Tier 1-2 - repeat-class half-life",
        "needs": ["repeat_class"],
        "inject": {"repeat_class": [("fig_repeat_class", "SNAP")]},
    },
    "family_dispersion": {
        "module": "fig_family_dispersion",
        "filename": "Fig_family_dispersion_threshold_sweep_rev.png",
        "title": "Tier 1-1b - family dispersion sweep",
        "needs": ["cohort", "pairtab"],
        "inject": {"cohort": [("family_lib", "COHORT")],
                   "pairtab": [("family_lib", "PAIRTAB")]},
    },
    "identity_boxplot": {
        "module": "fig_identity_boxplot",
        "filename": "Fig_identity_vs_halflife_boxplot_rev.png",
        "title": "Identity vs half-life difference",
        "needs": ["cohort", "pairtab"],
        "inject": {"cohort": [("family_lib", "COHORT")],
                   "pairtab": [("family_lib", "PAIRTAB")]},
    },
}


def figure_missing(dataset, fig_key):
    """Dataset keys this figure needs but does not have."""
    return [k for k in FIGURES[fig_key]["needs"] if not dataset.has(k)]


def figure_injections(dataset, fig_key):
    """[(module, attribute, path), ...] for this figure."""
    out = []
    for key, targets in FIGURES[fig_key]["inject"].items():
        path = dataset.paths.get(key)
        if not path:
            continue
        for mod, attr in targets:
            out.append((mod, attr, path))
    return out
