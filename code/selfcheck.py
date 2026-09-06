"""
selfcheck.py — verify the installation before (or instead of) launching Streamlit.

Run this first. It exercises everything the app depends on WITHOUT importing
Streamlit, so it separates "my data/code setup is broken" from "Streamlit itself
won't start" — two very different problems.

    python3 selfcheck.py            # fast checks (~1 min)
    python3 selfcheck.py --full     # also runs Tier 1-1b and 1-5 (a few minutes)

Exit code 0 = everything passed.
"""

import argparse
import os
import pathlib
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

PASS, FAIL, WARN = "  PASS", "  FAIL", "  WARN"
results = []


def check(name, fn, warn_only=False):
    try:
        detail = fn()
        results.append((PASS, name, detail or ""))
    except Exception as exc:                                     # noqa: BLE001
        results.append((WARN if warn_only else FAIL, name, f"{type(exc).__name__}: {exc}"))
        if os.environ.get("SELFCHECK_TRACE"):
            traceback.print_exc()


# ------------------------------------------------------------------- the checks

def c_python():
    v = sys.version_info
    if v < (3, 10):
        raise RuntimeError(f"Python {v.major}.{v.minor} is too old; need 3.10+")
    return f"Python {v.major}.{v.minor}.{v.micro}"


def c_numpy():
    import numpy
    return f"NumPy {numpy.__version__}"


def c_matplotlib():
    import matplotlib
    return f"matplotlib {matplotlib.__version__} (figures only)"


def c_streamlit():
    import streamlit
    ver = streamlit.__version__
    parts = [int(x) for x in ver.split(".")[:2]]
    if parts < [1, 39]:
        raise RuntimeError(
            f"Streamlit {ver} is too old — axis labels need 1.39+. "
            "Run: pip install --upgrade streamlit")
    return f"Streamlit {ver}"


def c_no_scipy():
    try:
        import scipy                                             # noqa: F401
    except ImportError:
        return "SciPy absent (correct — statistics are NumPy-only)"
    return "SciPy is installed; it is NOT used and must not be substituted"


def c_data_files():
    import dataset as ds
    d = ds.reference_dataset()
    missing = [k for k in d.paths if not d.has(k)]
    if missing:
        raise RuntimeError("missing reference data: " + ", ".join(missing))
    return f"all {len(d.paths)} reference data files present"


def c_integrity():
    import app_api as api
    ok, bad = api.verify_integrity()
    if bad:
        raise RuntimeError(f"{len(bad)} problem(s), first: {bad[0]}")
    return f"all {ok} files match MANIFEST.sha256"


def c_cohort():
    import app_api as api
    genes, hl, _ = api.load_cohort()
    if len(genes) != 3514:
        raise RuntimeError(f"expected 3514 proteins, got {len(genes)}")
    return f"{len(genes)} proteins loaded"


def c_tier(key, needle, expect, overrides=None):
    def inner():
        import app_api as api
        import dataset as ds
        res = api.run_tier(key, overrides=overrides or {},
                           dataset=ds.reference_dataset())
        line = next((l for l in res["stdout"].splitlines() if needle in l), None)
        if line is None:
            raise RuntimeError(f"no output line containing {needle!r}")
        if expect not in line:
            raise RuntimeError(f"expected {expect!r} in: {line.strip()}")
        return line.strip()[:78]
    return inner


def c_builders():
    import app_api as api
    import tempfile
    tmp = tempfile.mkdtemp()
    xlsx = os.path.join(PACK, "data", "inputs", "Monkey_half-life.xlsx")
    fasta = os.path.join(PACK, "data", "inputs", "cyno_target_rev.fasta")
    if not os.path.exists(xlsx):
        raise RuntimeError("Monkey_half-life.xlsx not present")
    coh = os.path.join(tmp, "c.tsv")
    api.build_cohort_from_uploads(xlsx, fasta, coh)
    pairs = os.path.join(PACK, "data", "inputs", "cyno_pairs.tsv")
    pt = os.path.join(tmp, "p.tsv")
    log, _ = api.build_pair_table_from_uploads(pairs, coh, pt)
    n = sum(1 for _ in open(pt, encoding="utf-8")) - 1
    if n != 11522:
        raise RuntimeError(f"pair table has {n} rows, expected 11522")
    return f"cohort + pair table rebuilt from raw inputs ({n} pairs)"


def c_annotate():
    import app_api as api
    import tempfile
    tmp = tempfile.mkdtemp()
    out = os.path.join(tmp, "sf.tsv")
    api.auto_seq_features(
        os.path.join(PACK, "data", "inputs", "cyno_target_rev.fasta"),
        os.path.join(PACK, "data", "inputs", "cyno_target_rev.tsv"), out)
    n = sum(1 for _ in open(out, encoding="utf-8")) - 1
    return f"sequence features generated for {n} proteins"


def c_encoding_hygiene():
    """Guard against the cp932 class of bug.

    On Japanese (or any non-UTF-8) Windows, a bare open() uses the OS locale
    codec and dies on UTF-8 data with e.g.
        'cp932' codec can't decode byte 0x93 ... illegal multibyte sequence
    Every text open() in this project must therefore pass encoding="utf-8",
    and no binary open() may pass an encoding at all.
    """
    import re
    bad = []
    for path in sorted(pathlib.Path(HERE).glob("*.py")) + [pathlib.Path(PACK) / "app.py"]:
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            # (?<![\w.]) skips urlopen(); an empty argument list is prose such as
            # "a bare open() uses the OS locale codec", not a real call.
            for m in re.finditer(r'(?<![\w.])open\(([^()]+)\)', line):
                inner = m.group(1)
                binary = '"rb"' in inner or '"wb"' in inner
                has_enc = "encoding=" in inner
                if binary and has_enc:
                    bad.append(f"{path.name}:{i} binary open with encoding=")
                elif not binary and not has_enc:
                    bad.append(f"{path.name}:{i} text open without encoding=")
    if bad:
        raise RuntimeError(f"{len(bad)} problem(s), first: {bad[0]}")
    return "every open() declares an explicit encoding"


def c_figures():
    """The figure scripts must still run after any I/O change."""
    import subprocess
    import tempfile
    out = tempfile.mkdtemp()
    env = dict(os.environ, FIG_OUTDIR=out)
    r = subprocess.run([sys.executable, os.path.join(HERE, "fig_repeat_class.py")],
                       capture_output=True, text=True, env=env, cwd=HERE, timeout=600)
    if r.returncode != 0:
        raise RuntimeError((r.stderr or r.stdout).strip().splitlines()[-1][:120])
    made = os.listdir(out)
    if not made:
        raise RuntimeError("figure script produced no file")
    return f"figure script OK ({made[0]})"


def c_app_syntax():
    import py_compile
    py_compile.compile(os.path.join(PACK, "app.py"), doraise=True)
    return "app.py compiles"


# -------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true",
                    help="also run the slow tiers (1-1b, 1-5)")
    args = ap.parse_args()

    print("=" * 74)
    print("Self-check — proteome half-life pipeline")
    print("=" * 74)

    print("\n[1] Runtime")
    check("Python version", c_python)
    check("NumPy", c_numpy)
    check("matplotlib", c_matplotlib, warn_only=True)
    check("Streamlit", c_streamlit, warn_only=True)
    check("SciPy absent", c_no_scipy, warn_only=True)

    print("\n[2] Data")
    check("Reference data files", c_data_files)
    check("File integrity", c_integrity, warn_only=True)
    check("Cohort loads", c_cohort)

    print("\n[3] Analysis (reference dataset)")
    check("Tier 1-1 paralog outliers", c_tier("1-1", "outlier pairs", "31"))
    check("Tier 1-2 repeat class", c_tier("1-2", "ALL_REPEAT", "1.984e-02"))
    check("Tier 1-4 sequence regression", c_tier("1-4", "whole-model", "0.0190"))
    check("Tier 1-6 STOICH vs REGUL", c_tier("1-6", "class sizes", "STOICH=127"))
    if args.full:
        check("Tier 1-5 CORUM (N_PERM=500)",
              c_tier("1-5", "members=", "1442", {"N_PERM": 500}))
        check("Tier 1-1b family dispersion",
              c_tier("1-1b", "50  ", "321", None))

    print("\n[4] Builders and generators")
    check("build_cohort + build_pair_table", c_builders)
    check("annotate seq-features", c_annotate)

    print("\n[5] App and portability")
    check("app.py syntax", c_app_syntax)
    check("Encoding hygiene (Windows-safe)", c_encoding_hygiene)
    if args.full:
        check("Figure scripts run", c_figures)

    print("\n" + "=" * 74)
    for status, name, detail in results:
        print(f"{status}  {name:38s} {detail}")
    print("=" * 74)

    failed = [r for r in results if r[0] == FAIL]
    warned = [r for r in results if r[0] == WARN]
    print(f"\n{len(results) - len(failed) - len(warned)} passed, "
          f"{len(warned)} warning(s), {len(failed)} failed")

    if failed:
        print("\nSomething is genuinely broken — fix the FAIL lines above.")
        print("Re-run with SELFCHECK_TRACE=1 for full tracebacks.")
        return 1
    if any(w[1] == "Streamlit" for w in warned):
        print("\nEverything works, but Streamlit is missing or too old.")
        print("Install it with:  pip install -r ../requirements.txt")
        print("Then launch with: cd .. && streamlit run app.py")
        return 0
    print("\nAll good. Launch the app with:")
    print("    cd .. && streamlit run app.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
