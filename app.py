"""
Proteome half-life — interactive analysis app.

Run:  streamlit run app.py

Ships with the pinned cynomolgus reference dataset, and accepts your own files
so the identical pipeline can be run on another species (mouse, rat, ...).

Every statistic shown here is computed by the tier scripts in code/ — the same
ones run_all.py executes and MANIFEST.md pins with an acceptance table. The app
adds interactivity and file handling, not analysis. See code/app_api.py.
"""

import os
import sys
import tempfile

import numpy as np
import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "code")
if CODE not in sys.path:
    sys.path.insert(0, CODE)

import app_api as api          # noqa: E402
import dataset as ds           # noqa: E402

st.set_page_config(page_title="Proteome half-life pipeline", layout="wide")

FIGDIR = os.path.join(HERE, "figures")
DOC = os.path.join(HERE, "REPRODUCIBILITY_RUNBOOK.md")

TIER_TITLES = {
    "1-1": "Tier 1-1 · Paralog outliers",
    "1-1b": "Tier 1-1b · Family dispersion",
    "1-2": "Tier 1-2 · Repeat class",
    "1-4": "Tier 1-4 · Sequence regression",
    "1-5": "Tier 1-5 · CORUM membership",
    "1-6": "Tier 1-6 · STOICH vs REGUL",
}


# ------------------------------------------------------------------ session state

def workdir():
    if "workdir" not in st.session_state:
        st.session_state.workdir = tempfile.mkdtemp(prefix="halflife_app_")
    return st.session_state.workdir


def user_paths():
    if "user_paths" not in st.session_state:
        st.session_state.user_paths = {}
    return st.session_state.user_paths


def save_upload(upload, filename):
    """Persist an uploaded file into the session working directory."""
    path = os.path.join(workdir(), filename)
    with open(path, "wb") as fh:
        fh.write(upload.getbuffer())
    return path


def active_dataset():
    if st.session_state.get("source") == "custom":
        return ds.Dataset("My species", dict(user_paths()))
    return ds.reference_dataset()


# ----------------------------------------------------------------------- helpers

@st.cache_data(show_spinner=False)
def cached_tier(tier_key, overrides_items, lib_items, path_items, ds_name):
    dataset = ds.Dataset(ds_name, dict(path_items))
    return api.run_tier(tier_key,
                        overrides=dict(overrides_items),
                        lib_overrides=dict(lib_items),
                        dataset=dataset)


def run(tier_key, overrides=(), lib_overrides=()):
    d = active_dataset()
    return cached_tier(tier_key, tuple(overrides), tuple(lib_overrides),
                       tuple(sorted(d.paths.items())), d.name)


def safe_chart(kind, data, x_label=None, y_label=None):
    """st.bar_chart / st.line_chart with axis labels.

    `x_label` / `y_label` were only added in Streamlit 1.39. requirements.txt
    pins >= 1.39, but fall back gracefully so the app still renders on an older
    install rather than crashing with a TypeError.
    """
    fn = st.bar_chart if kind == "bar" else st.line_chart
    try:
        fn(data, x_label=x_label, y_label=y_label)
    except TypeError:
        fn(data)


def table(header, rows, max_rows=None):
    if not header:
        st.info("No rows.")
        return
    shown = rows if max_rows is None else rows[:max_rows]
    st.dataframe({c: [r[i] if i < len(r) else "" for r in shown]
                  for i, c in enumerate(header)}, use_container_width=True)


def download_tsv(header, rows, filename):
    if not header:
        return
    body = "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows)
    st.download_button("Download TSV", body, file_name=filename,
                       mime="text/tab-separated-values")


def guard(tier_key):
    """Return the dataset if the tier can run, else render an explanation."""
    d = active_dataset()
    missing = d.missing_for(tier_key)
    if not missing:
        return d
    st.error(f"**{TIER_TITLES[tier_key]} cannot run with the current dataset.**")
    st.markdown("This tier needs the following file(s), which are not loaded:")
    for key in missing:
        label, cols, _, how = ds.SPEC[key]
        st.markdown(f"- **{label}** — {how}")
        st.caption("Required columns: `" + "`, `".join(cols) + "`")
    st.info("Load them on the **Data & Setup** page, or switch back to the "
            "pinned cynomolgus reference in the sidebar.")
    return None


# ------------------------------------------------------------------------ sidebar

st.sidebar.title("Proteome half-life")

st.session_state.setdefault("source", "reference")
st.sidebar.radio(
    "Dataset",
    ["reference", "custom"],
    format_func=lambda k: {"reference": "Cynomolgus (pinned reference)",
                           "custom": "My species (upload files)"}[k],
    key="source",
)

d = active_dataset()
avail = sorted(d.available_tiers())
if avail:
    st.sidebar.success(f"Tiers available: {', '.join(avail)}")
else:
    st.sidebar.warning("No tiers available yet — load files on **Data & Setup**.")

page = st.sidebar.radio("Section", [
    "Overview",
    "Data & Setup",
    "Input requirements",
    TIER_TITLES["1-1"],
    TIER_TITLES["1-1b"],
    TIER_TITLES["1-2"],
    TIER_TITLES["1-4"],
    TIER_TITLES["1-5"],
    TIER_TITLES["1-6"],
    "Figures",
    "Data export",
    "Reproducibility",
])

with st.sidebar.expander("Runtime"):
    for k, v in api.environment_report().items():
        st.write(f"**{k}** — {v}")
st.sidebar.caption("Statistics are NumPy-only by design; SciPy must not be used.")


# ----------------------------------------------------------------------- overview

if page == "Overview":
    st.title("Proteome half-life pipeline")
    st.markdown(
        "Does protein half-life track **sequence features**, **paralogy**, or "
        "**complex membership**? The pipeline ships with a pinned cynomolgus "
        "reference and runs unchanged on any species you supply files for."
    )

    genes, hl = [], []
    if d.reference:
        genes, hl, _ = api.load_cohort()
    elif d.has("cohort"):
        header, rows = api.read_tsv(d.paths["cohort"])
        j = {c: i for i, c in enumerate(header)}
        genes = [r[j["row_gene"]] for r in rows]
        hl = [float(r[j["half_life_h"]]) for r in rows]

    if genes:
        hl = np.array(hl)
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Proteins", f"{len(genes):,}")
        c2.metric("Median half-life", f"{np.median(hl):.2f} h")
        c3.metric("= in days", f"{np.median(hl) / 24:.2f} d")
        c4.metric("Range", f"{hl.min():.1f}–{hl.max():.0f} h")

        st.subheader("Half-life distribution")
        counts, _ = np.histogram(np.log10(hl), bins=48)
        safe_chart("bar", {"proteins": counts.tolist()},
                   x_label="log10(half-life, hours) — binned", y_label="proteins")
        st.caption("Half-life is handled internally in **hours**; figures convert "
                   "to days (day = h/24). `log2` for effect sizes, `log10` for "
                   "correlations.")
    else:
        st.info("No cohort loaded yet. Go to **Data & Setup**.")

    st.subheader("Reference findings (cynomolgus, n = 3,514)")
    st.markdown(
        "- **Sequence features barely explain half-life** — R² = 0.0190 (Tier 1-4).\n"
        "- **Complex membership is not a strong predictor** — within-complex "
        "permutation p ≈ 0.63 (Tier 1-5).\n"
        "- **Stoichiometric vs regulatory subunits do not differ** — p = 0.92 (Tier 1-6).\n"
        "- **Repeat-domain proteins are shorter-lived** — p = 1.984e-2 (Tier 1-2).\n"
        "- **Sequence-similar paralogs are co-stable** — permutation P < 1e-3 "
        "(Tier 1-1b), yet individual pairs differ up to ~9-fold (Tier 1-1)."
    )
    st.info("These null results motivate the cross-species comparison — which is "
            "what the custom-species mode is for.")


# ------------------------------------------------------------------ data & setup

elif page == "Data & Setup":
    st.title("Data & Setup")

    if d.reference:
        st.success("**Using the pinned cynomolgus reference.** All six tiers are "
                   "available and reproduce the published values.")
        st.markdown("Switch the sidebar **Dataset** to *My species* to load your "
                    "own files.")
        st.subheader("Reference files in use")
        for key, path in sorted(d.paths.items()):
            label = ds.SPEC[key][0]
            mark = "✅" if d.has(key) else "❌"
            st.markdown(f"{mark} **{label}** — `{os.path.relpath(path, HERE)}`")
    else:
        st.markdown(
            "Load your species' files below. **Three of them must be produced "
            "outside this app** (see *Input requirements*); the rest are built "
            "here or derived from UniProt/CORUM for your species."
        )

        up = user_paths()

        # ---- Step 1: cohort
        st.header("Step 1 — Cohort")
        st.caption("Needed by Tier 1-1 and 1-1b, and by the pair-table builder.")
        mode = st.radio("How do you want to provide the cohort?",
                        ["Build it from a half-life workbook + FASTA",
                         "Upload a ready-made cohort TSV"],
                        key="cohort_mode")

        if mode.startswith("Build"):
            c1, c2 = st.columns(2)
            xlsx = c1.file_uploader("Half-life workbook (.xlsx)", type=["xlsx"],
                                    key="u_xlsx")
            fasta = c2.file_uploader("Protein FASTA (>GENE|ACCESSION)",
                                     type=["fasta", "fa", "faa", "txt"],
                                     key="u_fasta")
            st.caption("Workbook: column 1 = gene name, last column = half-life in "
                       "**hours**. FASTA headers must be `>GENE|ACCESSION`.")
            if xlsx and fasta:
                if st.button("Build cohort", type="primary"):
                    xp = save_upload(xlsx, "halflife.xlsx")
                    fp = save_upload(fasta, "proteins.fasta")
                    out = os.path.join(workdir(), "cohort.tsv")
                    try:
                        log, path = api.build_cohort_from_uploads(xp, fp, out)
                        up["cohort"] = path
                        st.success("Cohort built.")
                        st.code(log)
                    except Exception as exc:                       # noqa: BLE001
                        st.error(f"Could not build the cohort: {exc}")
        else:
            f = st.file_uploader("Cohort TSV", type=["tsv", "txt"], key="u_cohort")
            if f:
                p = save_upload(f, "cohort.tsv")
                ok, msg = ds.check_columns(p, "cohort")
                (st.success if ok else st.error)(msg)
                if ok:
                    up["cohort"] = p

        # ---- Step 2: pair table
        st.header("Step 2 — Paralog pair table")
        st.caption("Needed by Tier 1-1 and 1-1b.")
        if "cohort" not in up:
            st.info("Load the cohort first.")
        else:
            pmode = st.radio("How do you want to provide the pair table?",
                             ["Build it from raw DIAMOND output",
                              "Upload a ready-made pair table TSV"],
                             key="pair_mode")
            if pmode.startswith("Build"):
                f = st.file_uploader("DIAMOND all-vs-all (.tsv)",
                                     type=["tsv", "txt"], key="u_diamond")
                st.caption("outfmt 6 with: `qseqid sseqid pident length qlen slen "
                           "evalue bitscore`. IDs must be `GENE|ACCESSION`.")
                if f and st.button("Build pair table", type="primary"):
                    dp = save_upload(f, "diamond.tsv")
                    out = os.path.join(workdir(), "pair_table.tsv")
                    try:
                        log, path = api.build_pair_table_from_uploads(
                            dp, up["cohort"], out)
                        up["pairtab"] = path
                        st.success("Pair table built.")
                        st.code(log)
                    except Exception as exc:                       # noqa: BLE001
                        st.error(f"Could not build the pair table: {exc}")
            else:
                f = st.file_uploader("Pair table TSV", type=["tsv", "txt"],
                                     key="u_pt")
                if f:
                    p = save_upload(f, "pair_table.tsv")
                    ok, msg = ds.check_columns(p, "pairtab")
                    (st.success if ok else st.error)(msg)
                    if ok:
                        up["pairtab"] = p

        # ---- Step 3: annotation tables
        st.header("Step 3 — Annotation tables")
        st.caption("Each unlocks one tier. **Most can now be generated for you** "
                   "— you no longer have to build them by hand.")

        st.subheader("Generate automatically")

        # --- Tier 1-4: fully offline from the FASTA
        with st.expander("Sequence features  →  Tier 1-4   ·   no network needed",
                         expanded=True):
            st.markdown(
                "Computed **from your FASTA alone**. Validated against the cyno "
                "reference: pI 3513/3514 exact, N-end class 3514/3514 exact, "
                "disorder fraction exact."
            )
            st.caption("`pest_max` uses a fresh PEST-FIND implementation and will "
                       "differ from the cyno snapshot. PEST is the weakest term in "
                       "the model (β = −0.0028, p = 0.65), so this does not change "
                       "the conclusion — auto-generated features give R² = 0.0189 "
                       "against a published 0.0190.")
            fa = st.file_uploader("Protein FASTA", type=["fasta", "fa", "faa", "txt"],
                                  key="a_fasta")
            if "cohort" not in up:
                st.info("Load the cohort first (Step 1).")
            elif fa and st.button("Generate sequence features", type="primary"):
                fp = save_upload(fa, "annot.fasta")
                out = os.path.join(workdir(), "seq_features.tsv")
                try:
                    log, path = api.auto_seq_features(
                        fp, up["cohort"], out, up.get("repeat_class"))
                    up["seq_features"] = path
                    st.success("Sequence features generated.")
                    st.code(log)
                except Exception as exc:                          # noqa: BLE001
                    st.error(f"Could not generate: {exc}")

        # --- Tier 1-2: UniProt domains
        with st.expander("Repeat-domain classes  →  Tier 1-2   ·   fetches UniProt"):
            st.markdown(
                "Maps UniProt domain/repeat annotations onto TPR, ANK, WD40, LRR, "
                "HEAT/ARM, KELCH/FBX and SPECTRIN classes."
            )
            st.caption("Fetching queries UniProt once over the network, then the "
                       "result is frozen as a file.")
            st.markdown(
                "**No internet on this machine?** Export the annotation from the "
                "UniProt website on any machine and upload the TSV below:\n\n"
                "1. Go to <https://www.uniprot.org/id-mapping> and paste your "
                "accession list (the `accession` column of your cohort).\n"
                "2. Map from **UniProtKB AC/ID** to **UniProtKB**.\n"
                "3. On the results page choose **Customize columns** and tick "
                "*Repeat*, *Domain [FT]*, *Domain [CC]* and *Protein families*.\n"
                "4. **Download → Format: TSV** (uncompressed), then upload it here.\n\n"
                "The first column must be the accession; every other column is "
                "scanned for repeat keywords, so extra columns are harmless."
            )
            uni = st.file_uploader("Pre-downloaded UniProt TSV (optional)",
                                   type=["tsv", "txt"], key="a_uniprot")
            if "cohort" not in up:
                st.info("Load the cohort first (Step 1).")
            elif st.button("Generate repeat classes"):
                out = os.path.join(workdir(), "repeat_class.tsv")
                up_tsv = save_upload(uni, "uniprot.tsv") if uni else None
                try:
                    with st.spinner("Querying UniProt…" if not up_tsv else "Working…"):
                        log, path = api.auto_repeat_class(up["cohort"], out, up_tsv)
                    up["repeat_class"] = path
                    st.success("Repeat classes generated.")
                    st.code(log)
                except Exception as exc:                          # noqa: BLE001
                    st.error(f"Could not generate: {exc}. If this is a network "
                             "error, download the UniProt TSV separately and "
                             "upload it above.")

        # --- Tier 1-5 / 1-6: CORUM
        with st.expander("CORUM tables  →  Tier 1-5 and 1-6   ·   no network needed"):
            st.markdown(
                "Generates all three CORUM tables at once. CORUM is a **human** "
                "resource and is reused for every species, keyed on the "
                "**uppercase gene symbol**."
            )
            st.caption("Download `corum_humanComplexes.txt` from "
                       "https://mips.helmholtz-muenchen.de/corum/download and "
                       "upload it here. It is deliberately not redistributed with "
                       "this repository.")
            st.warning(
                "Reproduces the pinned cyno stratification for **~97 %** of genes "
                "(3,409/3,514). Use it for a new species — do **not** use it to "
                "regenerate the cyno reference, which stays authoritative."
            )
            cf = st.file_uploader("corum_humanComplexes.txt", type=["txt", "tsv"],
                                  key="a_corum")
            if "cohort" not in up:
                st.info("Load the cohort first (Step 1).")
            elif cf and st.button("Generate CORUM tables"):
                cp = save_upload(cf, "corum.txt")
                m = os.path.join(workdir(), "corum_member.tsv")
                sfp = os.path.join(workdir(), "corum_cstats.tsv")
                t = os.path.join(workdir(), "corum_strat.tsv")
                try:
                    log, _ = api.auto_corum(cp, up["cohort"], m, sfp, t)
                    up["corum_member"], up["corum_cstats"], up["corum_strat"] = m, sfp, t
                    st.success("CORUM tables generated.")
                    st.code(log)
                    st.caption("Download them here, or any time from the "
                               "**Data export** page.")
                    for lbl, pth in [("corum_member_table.tsv", m),
                                     ("corum_complex_stats.tsv", sfp),
                                     ("corum_strat_assignment.tsv (Table S5)", t)]:
                        with open(pth, encoding="utf-8") as fh:
                            st.download_button(f"Download {lbl}", fh.read(),
                                               file_name=os.path.basename(pth),
                                               mime="text/tab-separated-values",
                                               key=f"dlc_{os.path.basename(pth)}")
                except Exception as exc:                          # noqa: BLE001
                    st.error(f"Could not generate: {exc}")

        st.subheader("Or upload your own")
        st.caption("Supplying a file overrides anything generated above.")
        for key in ["repeat_class", "seq_features", "corum_member",
                    "corum_cstats", "corum_strat"]:
            label, cols, tiers, how = ds.SPEC[key]
            done = "  ✅" if key in up else ""
            with st.expander(f"{label}  →  Tier {', '.join(tiers)}{done}"):
                st.caption(how)
                st.caption("Required columns: `" + "`, `".join(cols) + "`")
                f = st.file_uploader(label, type=["tsv", "txt"], key=f"u_{key}")
                if f:
                    p = save_upload(f, f"{key}.tsv")
                    ok, msg = ds.check_columns(p, key)
                    (st.success if ok else st.error)(msg)
                    if ok:
                        up[key] = p

        # ---- status
        st.header("Status")
        cur = ds.Dataset("My species", dict(up))
        for tier in ["1-1", "1-1b", "1-2", "1-4", "1-5", "1-6"]:
            miss = cur.missing_for(tier)
            if miss:
                names = ", ".join(ds.SPEC[m][0] for m in miss)
                st.markdown(f"❌ **Tier {tier}** — needs: {names}")
            else:
                st.markdown(f"✅ **Tier {tier}** — ready")
        if st.button("Clear all loaded files"):
            st.session_state.user_paths = {}
            st.rerun()


# ---------------------------------------------------------- input requirements

elif page == "Input requirements":
    st.title("What you need to prepare")
    st.markdown(
        "To run a **new species** you must supply the files below. Nothing is "
        "fetched from the internet — every external dependency is frozen, which "
        "is what keeps results stable over time."
    )

    st.header("① Produced outside this app — you must create these")
    st.markdown(
        "| # | File | Format | How to produce it |\n"
        "|---|---|---|---|\n"
        "| 1 | **Half-life workbook** | `.xlsx` | Column 1 = gene name "
        "(`;`-joined if peptide-ambiguous), last column = **half-life in hours**. "
        "Your proteomics turnover measurement. |\n"
        "| 2 | **Protein FASTA** | `.fasta` | Headers **must** be "
        "`>GENE|ACCESSION`, e.g. `>Aak1|A0A2K5V4D5`. One entry per gene in the "
        "workbook. |\n"
        "| 3 | **DIAMOND all-vs-all** | `.tsv` | DIAMOND is a compiled binary and "
        "cannot run inside this app. Produce it separately (command below). |\n"
    )
    st.code(
        "diamond makedb --in <species>.fasta -d <species>\n"
        "diamond blastp -q <species>.fasta -d <species> -o <species>_pairs.tsv \\\n"
        "        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \\\n"
        "        --max-target-seqs 50 --evalue 1e-3",
        language="bash")
    st.caption("The 8 columns and their order matter — the builder reads them "
               "positionally.")

    st.header("② Built for you by the app")
    st.markdown(
        "| Built file | From | Unlocks |\n|---|---|---|\n"
        "| Cohort table | workbook + FASTA | Tier 1-1, 1-1b |\n"
        "| Paralog pair table | DIAMOND output + cohort | Tier 1-1, 1-1b |\n"
    )
    st.success("**With files 1–3 alone you get Tier 1-1 and Tier 1-1b** — the "
               "paralog analyses. They need no external annotation.")

    st.header("③ Annotation tables — mostly generated for you")
    st.success(
        "**You no longer have to build these by hand.** The app generates all "
        "four from inputs you already have. Only CORUM needs one extra download."
    )
    st.markdown(
        "| Tier | Table | What it needs | Network? |\n|---|---|---|---|\n"
        "| 1-4 | Sequence features | your FASTA (already uploaded) | **No** |\n"
        "| 1-2 | Repeat-domain classes | UniProt domain annotation | Once, then frozen |\n"
        "| 1-5 | CORUM member table + complex stats | `corum_humanComplexes.txt` | **No** |\n"
        "| 1-6 | CORUM stratification | same CORUM file | **No** |\n"
    )
    st.markdown(
        "Download CORUM once from "
        "<https://mips.helmholtz-muenchen.de/corum/download>. It is a **human** "
        "resource, reused for every species via the uppercase gene symbol, so "
        "the same file works for mouse, rat and cyno alike. It is deliberately "
        "not redistributed with this repository."
    )
    st.caption(
        "Validation against the cyno reference — sequence features: pI 3513/3514 "
        "and N-end class 3514/3514 exact, disorder fraction exact, and the "
        "auto-generated table gives Tier 1-4 R² = 0.0189 against a published "
        "0.0190. CORUM stratification: ~97 % of genes match. You can still upload "
        "your own tables to override any of it."
    )
    for key in ["repeat_class", "seq_features", "corum_member", "corum_cstats",
                "corum_strat"]:
        label, cols, tiers, how = ds.SPEC[key]
        st.markdown(f"**{label}** → Tier {', '.join(tiers)}")
        st.caption(how + "  \nColumns: `" + "`, `".join(cols) + "`")

    st.header("Conventions you must follow")
    st.markdown(
        "- **Half-life is in HOURS** in every input. Figures convert to days.\n"
        "- **Gene symbols are the join key.** CORUM mapping uses the "
        "**uppercase** symbol; UniProt accession is the secondary key.\n"
        "- **Composite gene names** use `;` in the cohort (`Rack1;Gnb2L1`) and "
        "`_` in FASTA/DIAMOND (`Rack1_Gnb2L1`). Both are matched automatically.\n"
        "- **Duplicate accessions are deduplicated**: if two rows resolve to the "
        "same protein, the unambiguous (non-composite) row is kept."
    )
    st.warning("Tier 1-3 (GO comparison) is **descoped** and not part of this "
               "pipeline.")


# ---------------------------------------------------------------------- tier 1-1

elif page == TIER_TITLES["1-1"]:
    st.title(TIER_TITLES["1-1"])
    st.markdown(
        "Paralog pairs that are **sequence-similar but half-life-divergent**. "
        "The published cyno set uses `pident ≥ 80`, `≥ 2-fold`, `cov_min ≥ 0.5` "
        "and yields **31 pairs**."
    )
    if guard("1-1"):
        c1, c2, c3 = st.columns(3)
        pident = c1.slider("Min % identity", 30.0, 100.0, 80.0, 1.0)
        fold = c2.slider("Min fold difference", 1.5, 10.0, 2.0, 0.5)
        cov = c3.slider("Min mutual coverage", 0.0, 1.0, 0.5, 0.05)

        res = run("1-1", (("PIDENT_MIN", pident),
                          ("DLOG2_MIN", float(np.log2(fold))),
                          ("COV_MIN", cov)))
        n = len(res["rows"])
        a, b = st.columns([1, 3])
        a.metric("Outlier pairs", n)
        canonical = (abs(pident - 80) < 1e-9 and abs(fold - 2) < 1e-9
                     and abs(cov - 0.5) < 1e-9)
        if canonical and d.reference:
            (b.success if n == 31 else b.error)(
                "Canonical thresholds — expected 31 pairs." if n == 31
                else f"Canonical thresholds but {n} pairs — investigate.")
        table(res["header"], res["rows"])
        download_tsv(res["header"], res["rows"], "outlier_pairs.tsv")
        with st.expander("Script output"):
            st.code(res["stdout"])


# --------------------------------------------------------------------- tier 1-1b

elif page == TIER_TITLES["1-1b"]:
    st.title(TIER_TITLES["1-1b"])
    st.markdown(
        "Families are connected components of the pair graph (coverage ≥ 50 %, "
        "size ≥ 3). Within-family dispersion is the **population SD (ddof = 0)** "
        "of log10(half-life), compared with size-matched random sets."
    )
    if guard("1-1b"):
        ids = st.multiselect("Identity thresholds (%)",
                             [20, 30, 40, 50, 60, 70, 80],
                             default=[20, 30, 40, 50, 60, 70, 80])
        speed = st.select_slider("Null-model resolution",
                                 options=["fast", "standard", "canonical"],
                                 value="fast")
        n_sets = {"fast": 2000, "standard": 8000, "canonical": 20000}[speed]
        n_boot = {"fast": 1000, "standard": 3000, "canonical": 5000}[speed]
        if speed != "canonical":
            st.warning(
                f"Reduced null model ({n_sets:,} sets / {n_boot:,} bootstrap). "
                "Families, proteins and within-family SD are deterministic and "
                "unaffected; only `random_SD` and `perm_p` gain Monte-Carlo "
                "noise. Use **canonical** for quotable values."
            )
        if not ids:
            st.info("Select at least one identity threshold.")
        else:
            res = run("1-1b", (("IDENTITIES", tuple(sorted(ids))),),
                      (("N_SETS", n_sets), ("N_BOOT", n_boot)))
            table(res["header"], res["rows"])
            download_tsv(res["header"], res["rows"], "family_threshold_sweep.tsv")
            if res["rows"]:
                st.subheader("Within-family vs random-set dispersion")
                safe_chart("line",
                           {"within-family SD": [float(r[4]) for r in res["rows"]],
                            "random-set SD": [float(r[5]) for r in res["rows"]]},
                           x_label="identity threshold (row order)",
                           y_label="SD of log10(t½)")
            with st.expander("Script output"):
                st.code(res["stdout"])


# ---------------------------------------------------------------------- tier 1-2

elif page == TIER_TITLES["1-2"]:
    st.title(TIER_TITLES["1-2"])
    st.markdown("Do repeat-domain proteins (TPR, ANK, WD40, LRR, HEAT/ARM, "
                "KELCH/FBX, SPECTRIN) differ in half-life? Tie-corrected "
                "Mann–Whitney.")
    if guard("1-2"):
        res = run("1-2")
        table(res["header"], res["rows"])
        download_tsv(res["header"], res["rows"], "repeat_class_summary.tsv")
        st.code(res["stdout"])
        if d.reference:
            st.caption("Published: ALL_REPEAT n = 148, median 6.807 h, "
                       "p = 1.984e-2.")


# ---------------------------------------------------------------------- tier 1-4

elif page == TIER_TITLES["1-4"]:
    st.title(TIER_TITLES["1-4"])
    st.markdown(
        "OLS of log10(half-life) on log10(length), pI, disorder fraction, PEST "
        "max and three N-end-rule dummies. The feature set is **pinned to the "
        "snapshot** — a longer prose feature list does not reproduce R² = 0.0190."
    )
    if guard("1-4"):
        res = run("1-4")
        st.subheader("Summary")
        table(res["header"], res["rows"])
        download_tsv(res["header"], res["rows"], "seq_features_summary_rev.tsv")

        supp = res.get("extras", {}).get("Supp_Table_seq_features_rev.tsv")
        if supp:
            st.subheader("Supplementary Table S2")
            st.caption("One row per feature: Spearman rho/p alongside the OLS "
                       "beta, SE and p. This is the publication supplementary "
                       "table.")
            table(supp["header"], supp["rows"])
            download_tsv(supp["header"], supp["rows"],
                         "Supp_Table_seq_features_rev.tsv")
        st.code(res["stdout"])
        if d.reference:
            st.warning("Whole-model R² = 0.0190 — sequence features explain under "
                       "2 % of the variance. This is a null result, and it is the "
                       "point.")


# ---------------------------------------------------------------------- tier 1-5

elif page == TIER_TITLES["1-5"]:
    st.title(TIER_TITLES["1-5"])
    st.markdown(
        "Are complex members longer- or shorter-lived, and are subunits of one "
        "complex more similar than chance? The size-matched permutation null "
        "draws from the **whole cohort**, not the member subset — that choice "
        "moves p from ≈0.44 to ≈0.63."
    )
    if guard("1-5"):
        n_perm = st.select_slider("Permutations", options=[500, 2000, 10000],
                                  value=2000)
        if n_perm != 10000:
            st.warning(f"N_PERM = {n_perm:,}. Monte-Carlo error ≈ ±0.01–0.02; the "
                       "canonical value is **10,000** (cyno perm p = 0.6290). "
                       "Deterministic rows are unaffected.")
        else:
            st.info("Canonical setting — a few minutes on first run.")
        res = run("1-5", (("N_PERM", n_perm),))
        table(res["header"], res["rows"])
        download_tsv(res["header"], res["rows"], "corum_summary.tsv")
        st.code(res["stdout"])


# ---------------------------------------------------------------------- tier 1-6

elif page == TIER_TITLES["1-6"]:
    st.title(TIER_TITLES["1-6"])
    st.markdown("Proteins are stratified NONE / STOICH / REGUL / BOTH from "
                "complex-name patterns. If structural role drove turnover, "
                "STOICH and REGUL should differ.")
    if guard("1-6"):
        res = run("1-6")
        table(res["header"], res["rows"])
        download_tsv(res["header"], res["rows"], "corum_strat_summary.tsv")
        st.code(res["stdout"])
        if d.reference:
            st.warning("STOICH vs REGUL p = 0.92 — no difference. Another null "
                       "result.")


# ----------------------------------------------------------------------- figures

elif page == "Data export":
    st.title("Data export")
    st.markdown(
        "Every table the current dataset is built from, in one place. These are "
        "the **inputs** the tiers consume — including the ones that are also "
        "publication supplementary tables. Tier *results* are downloaded from "
        "each tier's own page."
    )
    if d.reference:
        st.caption("Showing the pinned cynomolgus reference tables.")
    else:
        st.caption("Showing the files you loaded or generated in this session.")

    SUPP = {
        "corum_strat": "Supplementary Table S5",
        "corum_cstats": "per-complex statistics (the 481-complex set)",
        "seq_features": "per-protein sequence features (source of Table S2)",
        "pairtab": "evaluated paralog pairs",
    }
    ORDER = ["cohort", "pairtab", "seq_features", "repeat_class",
             "corum_member", "corum_cstats", "corum_strat"]

    shown = 0
    for key in ORDER:
        label, cols, tiers, _ = ds.SPEC[key]
        if not d.has(key):
            st.markdown(f"❌ **{label}** — not loaded (Tier {', '.join(tiers)})")
            continue
        shown += 1
        path = d.paths[key]
        note = f" · {SUPP[key]}" if key in SUPP else ""
        st.markdown(f"✅ **{label}** — Tier {', '.join(tiers)}{note}")
        try:
            with open(path, encoding="utf-8") as fh:
                body = fh.read()
            n_rows = max(body.count("\n") - 1, 0)
            st.caption(f"`{os.path.basename(path)}` · {n_rows:,} rows · "
                       + ", ".join(cols[:4]) + (", …" if len(cols) > 4 else ""))
            st.download_button(f"Download {os.path.basename(path)}", body,
                               file_name=os.path.basename(path),
                               mime="text/tab-separated-values",
                               key=f"exp_{key}")
        except Exception as exc:                                    # noqa: BLE001
            st.error(f"Could not read the file: {exc}")

    if shown == 0:
        st.info("No tables loaded yet. Go to **Data & Setup**.")

    st.divider()
    st.subheader("Where these live in the repository")
    st.markdown(
        "| Table | Path |\n|---|---|\n"
        "| Cohort | `data/inputs/cyno_target_rev.tsv` |\n"
        "| Paralog pairs | `data/snapshots/pair_table_rev.tsv` |\n"
        "| Sequence features | `data/snapshots/seq_features_rev.tsv` |\n"
        "| Repeat classes | `data/snapshots/Repeat_class_halflife_rev.tsv` |\n"
        "| CORUM member table | `data/snapshots/corum_member_table_rev.tsv` |\n"
        "| CORUM complex stats | `data/snapshots/corum_complex_stats_rev.tsv` |\n"
        "| CORUM stratification (Table S5) | `data/snapshots/corum_strat_assignment_rev.tsv` |\n"
    )


elif page == "Figures":
    st.title("Figures")

    if d.reference:
        st.caption("Rendered from the pinned cynomolgus reference dataset.")
    else:
        st.info("Figures are rendered from **your** loaded files, so they follow "
                "your data. A figure whose input is missing is listed as "
                "unavailable rather than shown with reference data.")

    speed = st.select_slider(
        "Null-model resolution (family-dispersion figure only)",
        options=["fast", "standard", "canonical"], value="fast")
    n_sets = {"fast": 2000, "standard": 8000, "canonical": 20000}[speed]
    n_boot = {"fast": 1000, "standard": 3000, "canonical": 5000}[speed]

    order = ["repeat_class", "seq_features", "family_dispersion", "identity_boxplot"]
    todo, blocked = [], []
    for key in order:
        miss = ds.figure_missing(d, key)
        (blocked if miss else todo).append((key, miss))

    if blocked:
        st.subheader("Unavailable")
        for key, miss in blocked:
            names = ", ".join(ds.SPEC[m][0] for m in miss)
            st.markdown(f"- **{ds.FIGURES[key]['title']}** — needs: {names}")
        st.caption("Load the missing file(s) on the Data & Setup page.")

    if todo and st.button("Render figures", type="primary"):
        st.session_state["render_figs"] = True

    if st.session_state.get("render_figs") and todo:
        for key, _ in todo:
            spec = ds.FIGURES[key]
            st.subheader(spec["title"])
            try:
                with st.spinner(f"Rendering {spec['title']}…"):
                    res = api.render_figure(
                        key, dataset=d,
                        lib_overrides={"N_SETS": n_sets, "N_BOOT": n_boot})
                if not res["path"]:
                    st.error("The figure script produced no image.")
                    continue
                st.image(res["path"], use_container_width=True)
                with open(res["path"], "rb") as fh:
                    st.download_button(f"Download {spec['filename']}", fh.read(),
                                       file_name=spec["filename"],
                                       mime="image/png", key=f"dl_{key}")
                if res["stdout"].strip():
                    with st.expander("Script output"):
                        st.code(res["stdout"])
            except Exception as exc:                                # noqa: BLE001
                st.error(f"Could not render: {exc}")
    elif todo:
        st.caption("Press **Render figures** to generate them from the current "
                   "dataset.")

    with st.expander("Pinned reference renderings (shipped with the repository)"):
        st.caption("These are the published cynomolgus figures, for comparison.")
        for key in order:
            path = os.path.join(FIGDIR, ds.FIGURES[key]["filename"])
            if os.path.exists(path):
                st.image(path, caption=ds.FIGURES[key]["title"],
                         use_container_width=True)
            else:
                st.warning(f"Not found: `figures/{ds.FIGURES[key]['filename']}`")


# ---------------------------------------------------------------- reproducibility

elif page == "Reproducibility":
    st.title("Reproducibility")

    st.subheader("File integrity")
    if st.button("Verify checksums (MANIFEST.sha256)"):
        ok, bad = api.verify_integrity()
        if bad:
            st.error(f"{ok} file(s) OK, {len(bad)} problem(s):")
            for b in bad:
                st.write("- " + b)
        else:
            st.success(f"All {ok} files match their recorded SHA-256.")

    st.subheader("Acceptance targets (cynomolgus reference)")
    st.markdown(
        "| Tier | Headline | Status |\n|---|---|---|\n"
        "| 1-1 | 31 outlier pairs; Srsf6/Srsf4 fold 9.22 | exact |\n"
        "| 1-1b | families, proteins and within-family SD | exact |\n"
        "| 1-2 | ALL_REPEAT n=148, p = 1.984e-2 | exact |\n"
        "| 1-4 | R² = 0.0190; pI β = +0.0380 | exact |\n"
        "| 1-5 | 1442/2072; medians 7.6920/8.7744; perm p ≈ 0.63 | exact + MC |\n"
        "| 1-6 | 2072/127/1132/183; STOICH vs REGUL p = 0.92 | exact |\n"
    )
    st.caption("Full published-vs-recomputed table: `MANIFEST.md`.")

    st.subheader("Rules that keep results stable")
    st.markdown(
        "- Every permutation uses `numpy.default_rng(42)`.\n"
        "- Statistics are **NumPy-only**; `code/stats_numpy.py` is authoritative.\n"
        "- External databases are **frozen** — nothing is fetched live.\n"
        "- Half-life is internally in hours; figures convert to days.\n"
        "- Cross-species keys are **uppercase gene symbol**, accession secondary."
    )

    if os.path.exists(DOC):
        with st.expander("REPRODUCIBILITY_RUNBOOK.md"):
            with open(DOC, encoding="utf-8") as fh:
                st.markdown(fh.read())
