"""
app_api.py — thin adapter between the Streamlit app and the validated tier scripts.

DESIGN RULE (important for the paper)
-------------------------------------
This module contains **no statistics of its own**. Every number the app shows is
produced by the same tier scripts that `run_all.py` executes and that
`MANIFEST.md` pins with an acceptance table. The adapter only:

  * overrides module-level knobs (thresholds, N_PERM, identity list) so the app
    can be interactive;
  * redirects each script's output file to a scratch directory so the app never
    writes into the repository;
  * captures stdout and parses the produced TSV.

Because the tier modules expose their parameters as module-level constants, the
app can drive them without a fork. If a statistic ever needs changing, it changes
in ONE place — the tier script — and the app follows automatically. That is the
whole point of the reproducibility pack, extended to the UI.
"""

import contextlib
import importlib
import io
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
DATA = os.path.join(PACK, "data")

if HERE not in sys.path:
    sys.path.insert(0, HERE)

TIERS = {
    "1-1": ("tier1_1_paralog_outliers", "Paralog outliers"),
    "1-1b": ("tier1_1b_family_dispersion_sweep", "Family-dispersion identity sweep"),
    "1-2": ("tier1_2_repeat_class", "Repeat-class half-life"),
    "1-4": ("tier1_4_seq_regression", "Sequence-feature regression"),
    "1-5": ("tier1_5_corum_membership", "CORUM membership"),
    "1-6": ("tier1_6_stoich_regul", "STOICH vs REGUL"),
}


def read_tsv(path):
    """Return (header, rows) of a TSV as lists of strings."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return [], []
    return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]


def run_tier(tier_key, overrides=None, lib_overrides=None, dataset=None):
    """Execute a tier script with optional knob and data-path overrides.

    Parameters
    ----------
    tier_key : one of TIERS
    overrides : module-level constants to set on the tier module
                (e.g. {"PIDENT_MIN": 90.0} or {"N_PERM": 2000})
    lib_overrides : module-level constants to set on family_lib
                (e.g. {"N_SETS": 5000})
    dataset : a dataset.Dataset whose file paths are injected into the tier
              modules. When None the pinned cyno reference is used.

    Path injection is what lets a NEW SPECIES run through the identical,
    acceptance-tested tier code with no fork and no copy-paste.

    Returns
    -------
    dict with 'stdout', 'header', 'rows', 'out_path'
    """
    import dataset as ds_mod

    module_name, _ = TIERS[tier_key]
    mod = importlib.import_module(module_name)
    importlib.reload(mod)

    lib = importlib.import_module("family_lib")
    saved = []

    for k, v in (lib_overrides or {}).items():
        saved.append((lib, k, getattr(lib, k)))
        setattr(lib, k, v)

    if dataset is not None:
        missing = dataset.missing_for(tier_key)
        if missing:
            labels = ", ".join(ds_mod.SPEC[m][0] for m in missing)
            raise ValueError(
                f"Tier {tier_key} needs file(s) this dataset does not have: {labels}"
            )
        for target_mod, attr, path in dataset.injections(tier_key):
            tm = importlib.import_module(target_mod)
            saved.append((tm, attr, getattr(tm, attr)))
            setattr(tm, attr, path)

    for k, v in (overrides or {}).items():
        setattr(mod, k, v)

    # Redirect EVERY output constant (OUT, OUT_SUPP, ...) into a scratch dir so
    # the app never writes into the repository, and restore them afterwards.
    tmpdir = tempfile.mkdtemp(prefix="cyno_app_")
    out_attrs = [a for a in dir(mod)
                 if (a == "OUT" or a.startswith("OUT_"))
                 and isinstance(getattr(mod, a, None), str)]
    for a in out_attrs:
        saved.append((mod, a, getattr(mod, a)))
        setattr(mod, a, os.path.join(tmpdir, os.path.basename(getattr(mod, a))))
    out_path = mod.OUT

    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            mod.main()
    finally:
        for obj, attr, value in reversed(saved):
            setattr(obj, attr, value)

    header, rows = ([], [])
    if os.path.exists(out_path):
        header, rows = read_tsv(out_path)

    extras = {}
    for name in sorted(os.listdir(tmpdir)):
        path = os.path.join(tmpdir, name)
        if path != out_path and os.path.isfile(path):
            eh, er = read_tsv(path)
            extras[name] = {"header": eh, "rows": er, "path": path}

    return {"stdout": buf.getvalue(), "header": header, "rows": rows,
            "out_path": out_path, "extras": extras}


def build_cohort_from_uploads(xlsx_path, fasta_path, out_path):
    """Run build_cohort.py on user files. Returns (stdout, out_path)."""
    mod = importlib.import_module("build_cohort")
    importlib.reload(mod)
    rows, dropped, unresolved, n_raw = mod.build(xlsx_path, fasta_path)
    mod.write_tsv(rows, out_path)
    lines = [
        f"workbook rows                 = {n_raw}",
        f"cohort rows after dedup       = {len(rows)}",
        f"dropped (duplicate accession) = {len(dropped)}",
    ]
    for gene, acc, keeper in dropped:
        lines.append(f"    drop {gene!r} (accession {acc}) - redundant with {keeper!r}")
    if unresolved:
        lines.append(f"unresolved genes (no accession in FASTA) = {len(unresolved)}")
        for g in unresolved[:15]:
            lines.append(f"    {g}")
        if len(unresolved) > 15:
            lines.append(f"    ... and {len(unresolved) - 15} more")
    return "\n".join(lines), out_path


def build_pair_table_from_uploads(pairs_path, cohort_path, out_path):
    """Run build_pair_table.py on user files. Returns (stdout, out_path)."""
    mod = importlib.import_module("build_pair_table")
    importlib.reload(mod)
    rows, stats = mod.build(pairs_path, cohort_path)
    mod.write(rows, out_path)
    lines = [
        f"unordered pairs in DIAMOND output = {stats['raw_pairs']}",
        f"dropped (gene not in cohort)      = {stats['skipped_not_in_cohort']}",
        f"pairs written                     = {len(rows)}",
    ]
    if stats["malformed_lines"]:
        lines.append(f"malformed lines skipped           = {stats['malformed_lines']}")
    return "\n".join(lines), out_path


def load_cohort():
    """(genes, half_lives_hours) from the pinned cohort — for the overview page."""
    path = os.path.join(DATA, "inputs", "cyno_target_rev.tsv")
    header, rows = read_tsv(path)
    j = {c: i for i, c in enumerate(header)}
    genes = [r[j["row_gene"]] for r in rows]
    hl = [float(r[j["half_life_h"]]) for r in rows]
    lengths = [int(r[j["length"]]) for r in rows]
    return genes, hl, lengths


def expected_table(name):
    """Read an acceptance target from data/expected/."""
    path = os.path.join(DATA, "expected", name)
    if not os.path.exists(path):
        return [], []
    return read_tsv(path)


def verify_integrity():
    """Check every file listed in MANIFEST.sha256. Returns (n_ok, [failures])."""
    import hashlib
    manifest = os.path.join(PACK, "MANIFEST.sha256")
    if not os.path.exists(manifest):
        return 0, ["MANIFEST.sha256 not found"]
    ok, bad = 0, []
    with open(manifest, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            digest, _, rel = line.partition("  ")
            target = os.path.join(PACK, rel)
            if not os.path.exists(target):
                bad.append(f"{rel}: MISSING")
                continue
            h = hashlib.sha256(open(target, "rb").read()).hexdigest()
            if h == digest:
                ok += 1
            else:
                bad.append(f"{rel}: CHECKSUM MISMATCH")
    return ok, bad


def environment_report():
    """Runtime provenance for the sidebar."""
    import numpy
    rep = {"Python": sys.version.split()[0], "NumPy": numpy.__version__}
    try:
        import matplotlib
        rep["matplotlib"] = matplotlib.__version__
    except Exception:
        rep["matplotlib"] = "not installed"
    try:
        import scipy  # noqa: F401
        rep["SciPy"] = "PRESENT (must not be used for statistics)"
    except Exception:
        rep["SciPy"] = "absent (correct)"
    return rep


# ---------------------------------------------------- automatic annotation

def auto_seq_features(fasta_path, cohort_path, out_path, repeat_class_path=None):
    """Generate the Tier 1-4 feature table from sequence alone. No network."""
    import argparse as _a
    mod = importlib.import_module("annotate")
    importlib.reload(mod)
    args = _a.Namespace(fasta=fasta_path, cohort=cohort_path, out=out_path,
                        repeat_class=repeat_class_path)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mod.cmd_seq_features(args)
    return buf.getvalue(), out_path


def auto_repeat_class(cohort_path, out_path, uniprot_tsv=None):
    """Generate the Tier 1-2 repeat-class table. Fetches UniProt unless a
    pre-downloaded TSV is supplied."""
    import argparse as _a
    mod = importlib.import_module("annotate")
    importlib.reload(mod)
    args = _a.Namespace(cohort=cohort_path, out=out_path, uniprot_tsv=uniprot_tsv)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mod.cmd_repeat_class(args)
    return buf.getvalue(), out_path


def auto_corum(corum_path, cohort_path, out_member, out_stats, out_strat,
               organism="Human"):
    """Generate the Tier 1-5 and 1-6 tables from a CORUM download. No network."""
    import argparse as _a
    mod = importlib.import_module("annotate")
    importlib.reload(mod)
    args = _a.Namespace(corum=corum_path, cohort=cohort_path,
                        out_member=out_member, out_stats=out_stats,
                        out_strat=out_strat, organism=organism)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mod.cmd_corum(args)
    return buf.getvalue(), (out_member, out_stats, out_strat)


def render_figure(fig_key, dataset=None, lib_overrides=None, outdir=None):
    """Render one figure FROM THE ACTIVE DATASET.

    The figure scripts hold their input path and their OUTDIR as module-level
    constants, so the same injection trick used for the tiers works here: point
    them at this dataset's files, point OUTDIR at a scratch directory, run
    main(), then restore everything. Nothing is written into the repository.

    Returns {"path": <png>, "stdout": ...}.
    """
    import dataset as ds_mod

    spec = ds_mod.FIGURES[fig_key]
    if dataset is not None:
        missing = ds_mod.figure_missing(dataset, fig_key)
        if missing:
            labels = ", ".join(ds_mod.SPEC[m][0] for m in missing)
            raise ValueError(
                f"{spec['title']} needs file(s) this dataset does not have: {labels}")

    mod = importlib.import_module(spec["module"])
    importlib.reload(mod)

    saved = []
    lib = importlib.import_module("family_lib")
    for k, v in (lib_overrides or {}).items():
        saved.append((lib, k, getattr(lib, k)))
        setattr(lib, k, v)

    if dataset is not None:
        for target_mod, attr, path in ds_mod.figure_injections(dataset, fig_key):
            tm = importlib.import_module(target_mod)
            saved.append((tm, attr, getattr(tm, attr)))
            setattr(tm, attr, path)

    tmpdir = outdir or tempfile.mkdtemp(prefix="cyno_fig_")
    os.makedirs(tmpdir, exist_ok=True)
    saved.append((mod, "OUTDIR", getattr(mod, "OUTDIR")))
    mod.OUTDIR = tmpdir

    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            mod.main()
    finally:
        for obj, attr, value in reversed(saved):
            setattr(obj, attr, value)

    path = os.path.join(tmpdir, spec["filename"])
    return {"path": path if os.path.exists(path) else None,
            "stdout": buf.getvalue()}
