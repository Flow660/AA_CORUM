# Cynomolgus proteome half-life — analysis pipeline

Reproducible analysis of protein half-life in the cynomolgus monkey (*Macaca
fascicularis*) proteome (**n = 3,514 proteins**), testing whether turnover is
predicted by **sequence features**, **paralogy**, or **protein-complex
membership** — with an interactive Streamlit app and a fully pinned,
offline-reproducible pipeline.

**Owner:** Kosugi, Youhei — Preclinical Biomarkers, Takeda

---

## Summary of findings

| Tier | Question | Headline result |
|---|---|---|
| **1-1** | Do sequence-similar paralogs share a half-life? | **31** paralog pairs are ≥ 2-fold divergent despite ≥ 80 % identity; top pair Srsf6/Srsf4 at **9.22-fold** |
| **1-1b** | Are paralog families co-stable? | Yes — within-family dispersion is far below size-matched random sets, **permutation P < 1e-3** at every identity threshold 20–70 |
| **1-2** | Are repeat-domain proteins different? | Repeat proteins are **shorter-lived**: ALL_REPEAT n = 148, median 6.807 h vs 8.414 h, **p = 1.984e-2** |
| **1-4** | Do sequence features predict half-life? | **Barely.** Whole-model **R² = 0.0190** — under 2 % of variance. pI β = +0.0380 (p = 1.07e-9), disorder β = −0.0294 |
| **1-5** | Do complex subunits share a half-life? | Members are modestly shorter-lived (7.69 h vs 8.77 h, MW p = 7.3e-6), but **within-complex similarity is not better than chance: permutation p ≈ 0.63** |
| **1-6** | Does structural role matter? | **No.** Stoichiometric vs regulatory subunits: **p = 0.92** |

**The negative results are the point.** Sequence features, complex membership and
structural role all fail to explain protein half-life. This motivates the
cross-species work (Tier 1-7).

> Tier 1-3 (GO category comparison) was **descoped**. It is not part of this
> repository and is not a reproducibility gap.

---

## Quick start

```bash
git clone <this-repo>
cd <this-repo>
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

**Check the installation first**

```bash
cd code
python3 selfcheck.py          # ~1 min
python3 selfcheck.py --full   # also runs the slow tiers
```

`selfcheck.py` imports **no Streamlit**, so it cleanly separates "my data or code
setup is broken" from "Streamlit itself won't start". It verifies the runtime,
the data files, the checksums, four tiers against their expected values, both
builders and the annotation generator. Exit code 0 means everything passed.

**Interactive app**

```bash
streamlit run app.py          # from the repository root
```

Then open <http://localhost:8501> (Streamlit usually opens it for you). Stop it
with `Ctrl-C`.

Requires **Streamlit ≥ 1.39** — `st.bar_chart`/`st.line_chart` axis labels were
added in that release. `selfcheck.py` checks the version explicitly and tells you
if it is too old.

### If the app does not start

| Symptom | Cause | Fix |
|---|---|---|
| `streamlit: command not found` | not installed, or a different Python env is active | `pip install -r requirements.txt` |
| `ModuleNotFoundError: app_api` | launched from the wrong directory | run `streamlit run app.py` from the **repository root**, not from `code/` |
| `TypeError: ... x_label` | Streamlit older than 1.39 | `pip install --upgrade streamlit` |
| Tiers show "cannot run" | you are in *My species* mode with files missing | switch the sidebar to the pinned reference, or load the files |
| Port already in use | another Streamlit is running | `streamlit run app.py --server.port 8502` |
| `'cp932' codec can't decode byte ...` | fixed in this version — every file is read as UTF-8 regardless of the OS locale | update to the current code; `selfcheck.py` enforces it |

### A note on Windows and character encoding

All data files here are **UTF-8**. On a Japanese (or any non-UTF-8) Windows,
Python's bare `open()` defaults to the OS codec — `cp932` — and dies on UTF-8
bytes with `illegal multibyte sequence`. Every file operation in this project
therefore declares `encoding="utf-8"` explicitly, and `selfcheck.py` has an
**Encoding hygiene** check that fails the build if a bare `open()` ever creeps
back in.

**Command line — run every tier and print recomputed vs expected**

```bash
cd code
python3 run_all.py                 # canonical (N_PERM=10000, a few minutes)
N_PERM=2000 python3 run_all.py     # fast check (~1 min)
```

**Rebuild the cohort from the raw workbook**

```bash
cd code
python3 build_cohort.py            # verifies byte-identical output
```

**Regenerate figures**

```bash
cd code
FIG_OUTDIR=../figures python3 fig_repeat_class.py
FIG_OUTDIR=../figures python3 fig_seq_features.py
FIG_OUTDIR=../figures python3 fig_family_dispersion.py
FIG_OUTDIR=../figures python3 fig_identity_boxplot.py
```

**Verify file integrity**

```bash
sha256sum -c MANIFEST.sha256
```

No network access is required at any point.

---

## What to store in the repository

This is the section to get right before pushing. Total tracked size is
**≈ 5.7 MB**, largest single file 2.4 MB — **plain Git is fine, Git LFS is not
needed.**

### ✅ Commit these — 47 files

| Path | Files | Why it must be committed |
|---|---|---|
| `README.md` | 1 | This file |
| `LICENSE` | 1 | BSD 3-Clause — reuse terms for the code |
| `CITATION.cff` | 1 | Machine-readable citation for GitHub/Zenodo |
| `requirements.txt` | 1 | Pins the runtime |
| `.gitignore` | 1 | Keeps regenerated output out of history |
| `app.py` | 1 | Streamlit entry point |
| `MANIFEST.sha256` | 1 | SHA-256 of every stable artifact — the integrity check |
| `MANIFEST.md` | 1 | **Acceptance table**: published vs recomputed, per tier |
| `REPRODUCIBILITY_RUNBOOK.md` | 1 | Methods, recovered rules, and every pinned decision |
| `ENVIRONMENT.md` | 1 | Versions, seed, no-network / no-SciPy constraints |
| `RECOVERY.md` | 1 | 4-step restore procedure |
| `code/` | 19 | All analysis code (see below) |
| `data/inputs/` | 4 | Pinned raw inputs |
| `data/snapshots/` | 6 | **Frozen external-database captures** |
| `data/expected/` | 8 | Acceptance targets (the published values) |

**`code/` (19 files)**

```
stats_numpy.py                      NumPy-only ranks, Spearman, Mann-Whitney,
                                    Cliff's delta, OLS, Student-t p  ← authoritative
family_lib.py                       shared seeded engine for family dispersion
build_cohort.py                     raw .xlsx -> cohort tsv (upstream data prep)
build_pair_table.py                 raw DIAMOND -> evaluated pair table
annotate.py                         generates the Tier 1-2/1-4/1-5/1-6 tables
selfcheck.py                        verifies the install without Streamlit
dataset.py                          per-species file set + which tiers it unlocks
app_api.py                          adapter: Streamlit -> tier scripts
run_all.py                          runs every tier, prints recomputed vs expected
tier1_1_paralog_outliers.py
tier1_1b_family_dispersion_sweep.py
tier1_2_repeat_class.py
tier1_4_seq_regression.py
tier1_5_corum_membership.py
tier1_6_stoich_regul.py
fig_repeat_class.py
fig_seq_features.py
fig_family_dispersion.py
fig_identity_boxplot.py
```

**`data/inputs/` (4 files) — pinned raw inputs**

| File | Size | Notes |
|---|---|---|
| `Monkey_half-life.xlsx` | 140 KB | Raw half-life workbook, 3,515 rows |
| `cyno_target_rev.tsv` | 136 KB | Derived cohort, n = 3,514 |
| `cyno_target_rev.fasta` | 2.4 MB | Also the frozen gene → cyno accession map |
| `cyno_pairs.tsv` | 1.7 MB | DIAMOND all-vs-all output |

**`data/snapshots/` (6 files) — the reproducibility keystone**

These freeze the state of live databases at the time of the original run. Without
them the pipeline is **not** reproducible, because UniProt / CORUM change daily.

```
seq_features_rev.tsv            UniProt-derived features per protein
Repeat_class_halflife_rev.tsv   UniProt domain -> repeat class
corum_member_table_rev.tsv      member / non-member, promiscuity
corum_strat_assignment_rev.tsv  NONE / STOICH / REGUL / BOTH
corum_complex_stats_rev.tsv     the 481-complex set, per-complex SD
pair_table_rev.tsv              11,527 evaluated paralog pairs
```

**`data/expected/` (8 files)** — the published values each tier is checked
against. `run_all.py` prints `recomputed (expected …)` for every one.

### ❌ Do not commit

| Path | Why |
|---|---|
| `results/` | Regenerated by `run_all.py` — reproducible, not source of truth. Already in `.gitignore` |
| `__pycache__/`, `*.pyc` | Build artifacts |
| `.venv/`, `venv/` | Environment |
| `working/`, `baseline/`, `*.log` | Scratch |

### 🤔 Optional — your call

| Path | Recommendation |
|---|---|
| `figures/` (4 PNGs) | **Commit if** the journal needs the exact submitted renderings; otherwise leave them regenerated. Currently commented out in `.gitignore` — uncomment that line to exclude them |
| `Cyno_pipeline_runbook_v2_reproducible.pdf` | A rendering of `REPRODUCIBILITY_RUNBOOK.md`. Redundant with the Markdown; drop it unless you want the typeset version archived |

---

## Repository layout

```
.
├── README.md                       ← you are here
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── .gitignore
├── app.py                          ← streamlit run app.py
├── MANIFEST.sha256                 ← sha256sum -c MANIFEST.sha256
├── MANIFEST.md                     ← acceptance table
├── REPRODUCIBILITY_RUNBOOK.md      ← methods + every recovered rule
├── ENVIRONMENT.md                  ← versions, seed, constraints
├── RECOVERY.md                     ← restore procedure
├── code/                           ← 19 files, all analysis code
├── data/
│   ├── inputs/                     ← 4 pinned raw inputs
│   ├── snapshots/                  ← 6 frozen external-DB captures
│   └── expected/                   ← 8 acceptance targets
├── figures/                        ← 4 PNGs (regenerable)
└── results/                        ← gitignored; created by run_all.py
```

`code/` and `data/` are siblings and the scripts resolve paths relative to
themselves. **Do not rename or reorganise these two directories** — it will break
every path and the checksum manifest.

---

## The app

`streamlit run app.py` gives twelve sections: an overview, a **Data & Setup**
page for loading your own species' files, an **Input requirements** page, one
page per tier, a figure gallery that **re-renders from the active dataset**, a
**Data export** page that offers every input table as a download (including
Supplementary Table S5), and a reproducibility page with a live checksum
verifier.

The sidebar **Dataset** switch chooses between the pinned cynomolgus reference
and your own uploads, and continuously reports which tiers are currently
available.

Tier 1-1 exposes interactive identity / fold-change / coverage thresholds, Tier
1-1b lets you pick identity thresholds and null-model resolution, and Tier 1-5
lets you trade permutation count against runtime.

**The app contains no statistics of its own.** `code/app_api.py` drives the same
tier scripts that `run_all.py` runs, overriding module-level constants — including
the **data-path constants**, which is how a new species runs through unchanged
code with no fork. Output goes to a scratch directory, so the app never writes
into the repository, and every overridden value is restored afterwards.

> When you lower a Monte-Carlo setting the app tells you so and marks the values
> as non-quotable. Deterministic results are unaffected by those controls.

---

## Reproducibility guarantees

- **Seed** — every permutation uses `numpy.default_rng(42)`.
- **NumPy only** — `code/stats_numpy.py` re-implements ranks, Spearman,
  tie-corrected Mann–Whitney, Cliff's delta, OLS and the Student-t p-value.
  **Do not substitute SciPy**; published values are pinned to these
  implementations.
- **No live database access** — everything external is frozen in
  `data/snapshots/`, so results do not drift with the state of UniProt or CORUM
  on the day you run them.
- **Units** — half-life is handled internally in **hours**; figures convert to
  days (day = h/24). `log2` for effect sizes, `log10` for correlations.
- **Keys** — cross-species comparisons key on **uppercase gene symbol**, with
  UniProt accession as the secondary key.

Validated on Python 3.12 / NumPy 2.5 / matplotlib 3.11, with SciPy absent.
Minor version drift is tolerated; see `ENVIRONMENT.md`.

### Acceptance

`run_all.py` prints recomputed vs expected for every tier. Deterministic rows
must match to the printed precision; the single permutation row must fall within
Monte-Carlo error (± ~0.01 at N_PERM = 10,000). The full table is in
`MANIFEST.md`.

---

## Methodological notes worth knowing

Four details were **not** recoverable from prose and are now pinned in code.
Each one silently changes the answer if you get it wrong.

1. **The Tier 1-4 feature set is smaller than it looks.** The model that yields
   R² = 0.0190 uses only log10(length), pI, disorder fraction, PEST max and
   three N-end dummies — not the longer list of MW / GRAVY / instability index /
   TM count. It must come from the frozen `seq_features_rev.tsv`.

2. **The Tier 1-5 permutation pool is the whole cohort.** The size-matched null
   draws from all 3,514 proteins, not the 1,442 members. Member pool → null
   1.073, p ≈ 0.44; whole-cohort pool → 1.063, p ≈ 0.63. One unstated choice,
   opposite conclusions about how close to significance the result sits.

3. **Tier 1-1 needs the saved intermediate plus a coverage cut.** Rebuilding the
   pair table from raw DIAMOND diverges (24–35 pairs) because the coverage / HSP /
   best-hit rules were never written down. Using the frozen `pair_table_rev.tsv`
   with the rev-cohort gene filter and `pident ≥ 80, 2-fold, cov_min ≥ 0.5`
   yields exactly the 31 published pairs.

4. **Composite gene names must be normalised.** The pair table joins composite
   genes with `_` (`Rack1_Gnb2L1`) while the cohort uses `;` (`Rack1;Gnb2L1`).
   A naive membership test silently drops every pair touching a composite gene.
   Normalising is what makes Tier 1-1b reproduce the published families,
   proteins-in-families and within-family SD exactly. See
   `REPRODUCIBILITY_RUNBOOK.md` §9.

### The cohort dedup rule

`build_cohort.py` reduces 3,515 workbook rows to the 3,514-protein cohort by
dropping exactly one row. The workbook contains **seven** semicolon composites
and six are kept, so the rule is *not* "drop composite gene names". It is
**collision by accession**:

```
Mrc1           -> A0A2K5TN91   (half-life  5.996546499 h)
Mrc1;Slc39A12  -> A0A2K5TN91   (half-life 15.40101596  h)   <- same protein
```

`Mrc1` exists as its own confident row, so the composite is a redundant second
measurement of the same protein. Where several rows resolve to one accession,
the confident (non-composite) assignment is kept. The pre-dedup pair table still
carries the pair `Mrc1 / Mrc1_Slc39A12` at `pident = 100.00, cov_min = 1.000` —
a self-pair, which is the corroborating evidence.

---

## Data provenance and licensing

The **code** is released under the **BSD 3-Clause** licence (`LICENSE`).
Third-party data redistributed here keeps its own terms — **verify these before
making the repository public**:

| Data | Source | Note |
|---|---|---|
| CORUM complexes | CORUM (Helmholtz Munich), release 5.x, Human | **NOT redistributed** — the raw file is deliberately excluded (see below). Only derived, aggregated tables are shipped |
| `data/snapshots/*` (UniProt-derived) | UniProt | UniProt data is CC BY 4.0 — attribute it |
| `cyno_target_rev.fasta`, `cyno_pairs.tsv` | UniProt sequences; DIAMOND all-vs-all | DIAMOND is a compiled binary and is **not** run in this repository; its output is pinned |
| `Monkey_half-life.xlsx` | In-house measurements | **Confirm internal clearance to publish before pushing** |

### Why the raw CORUM file is not in this repository

`corum_humanComplexes.txt` (CORUM 5.x, Human) was **removed on purpose**. CORUM's
public pages give citation instructions but do not state explicit redistribution
terms, so shipping a verbatim copy would need clearance we do not have.

Removing it costs nothing: **no script in this repository ever read that file.**
The analyses consume the derived, aggregated tables in `data/snapshots/`
(`corum_member_table_rev.tsv`, `corum_complex_stats_rev.tsv`,
`corum_strat_assignment_rev.tsv`), not the raw database. Reproducibility is
unaffected and the repository is ~6.3 MB smaller.

For provenance, the exact file used was:

```
CORUM 5.x, Human complexes (corum_humanComplexes.txt)
7,813 lines
SHA-256: 5e40556eb59bb767f3396ba8c080f64c777903ffca9102345fdc35662a2e0282
```

Download it from <https://mips.helmholtz-muenchen.de/corum/download> and cite:
Tsitsiridis G, et al. *CORUM: the comprehensive resource of mammalian protein
complexes–2022.* Nucleic Acids Res. 2022. doi:10.1093/nar/gkac1015

> ⚠️ **Two release dates are still unrecorded** — the exact CORUM release and the
> UniProt release used for the snapshots. Add them to `ENVIRONMENT.md` before
> submission; reviewers will reasonably ask.

---

## Known limitations

- **Tier 1-3 (GO comparison) is descoped** and absent from this repository.
- **DIAMOND is not run here.** For a new species you must produce
  `<species>_pairs.tsv` outside this repository and add it to `data/inputs/`.
- **Snapshots are cyno-specific.** A new species needs its own UniProt features,
  domains and complex mappings, fetched live and then frozen the same way.
- **`random_SD` and `perm_p` in Tier 1-1b are Monte-Carlo** and are additionally
  sensitive to family ordering, because `null_sds()` cycles the family-size array
  with `np.resize(sizes, n_sets)`. Deterministic columns are unaffected.

## Running another species

The app runs any species through the **identical, acceptance-tested tier code**
— open it and switch **Dataset → My species**, or use the command line.

### What you must prepare

**Produced outside this repository (3 files):**

| # | File | Format |
|---|---|---|
| 1 | Half-life workbook | `.xlsx` — column 1 = gene name (`;`-joined if peptide-ambiguous), last column = **half-life in HOURS** |
| 2 | Protein FASTA | headers **must** be `>GENE\|ACCESSION` |
| 3 | DIAMOND all-vs-all | `.tsv`, outfmt 6 (command below) |

```bash
diamond makedb --in <species>.fasta -d <species>
diamond blastp -q <species>.fasta -d <species> -o <species>_pairs.tsv \
        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
        --max-target-seqs 50 --evalue 1e-3
```

**Built for you from those three:**

```bash
cd code
python3 build_cohort.py     --xlsx <species>.xlsx --fasta <species>.fasta \
                            --out ../data/inputs/<species>_target.tsv --no-verify
python3 build_pair_table.py --pairs <species>_pairs.tsv \
                            --cohort ../data/inputs/<species>_target.tsv \
                            --out ../data/snapshots/<species>_pair_table.tsv
```

> **Files 1–3 alone unlock Tier 1-1 and Tier 1-1b.** The paralog analyses need no
> external annotation.

**The remaining tiers are generated for you** by `code/annotate.py` — you do not
have to build these tables by hand.

| Tier | Table | Needs | Network? |
|---|---|---|---|
| 1-4 | Sequence features | your FASTA | **No** |
| 1-2 | Repeat-domain classes | UniProt domain annotation | Once, then frozen |
| 1-5 | CORUM member table + complex stats | `corum_humanComplexes.txt` | **No** |
| 1-6 | CORUM stratification | same CORUM file | **No** |

```bash
python3 annotate.py seq-features --fasta <species>.fasta \
        --cohort ../data/inputs/<species>_target.tsv --out <species>_seq_features.tsv
python3 annotate.py repeat-class --cohort ../data/inputs/<species>_target.tsv \
        --out <species>_repeat_class.tsv          # add --uniprot-tsv to work offline
python3 annotate.py corum --corum corum_humanComplexes.txt \
        --cohort ../data/inputs/<species>_target.tsv \
        --out-member <s>_member.tsv --out-stats <s>_cstats.tsv --out-strat <s>_strat.tsv
```

CORUM is a **human** resource reused for every species via the uppercase gene
symbol, so one download serves mouse, rat and cyno alike. Download it from
<https://mips.helmholtz-muenchen.de/corum/download>.

#### Getting the UniProt annotation for Tier 1-2

`repeat-class` classifies each protein into TPR / ANK / WD40 / LRR / HEAT_ARM /
KELCH_FBX / SPECTRIN / OTHER from UniProt domain and repeat annotations. Either
let the script fetch it, or export it yourself and pass `--uniprot-tsv`:

1. <https://www.uniprot.org/id-mapping> — paste the `accession` column of your cohort.
2. Map **UniProtKB AC/ID → UniProtKB**.
3. **Customize columns** → tick *Repeat*, *Domain [FT]*, *Domain [CC]*, *Protein families*.
4. **Download → TSV**, uncompressed.

The reader takes column 1 as the accession and keyword-scans every other column,
so extra columns are harmless and column order beyond the first does not matter.

### How much of this is validated?

Measured against the pinned cyno snapshots (`python3 annotate.py validate`):

| Feature | Result |
|---|---|
| `nend_residue`, `nend_class` | **3514/3514 exact** (Met-excision rule recovered) |
| `pI` | **3513/3514 exact** — Bjellqvist/Biopython pKa set with residue-specific terminal pKa |
| `disorder_frac` | **exact** — it is the fraction of {A,R,G,Q,S,P,E,K}, a composition metric rather than a disorder predictor |
| `length` | exact |
| `pest_max` | **not reproducible** — the original region-splitting rule was not recovered |
| CORUM stratification | **~97 %** of genes (3,409/3,514) |

Despite `pest_max` differing, a **fully auto-generated** feature table gives Tier
1-4 **R² = 0.0189** against a published 0.0190, with pI β = +0.0379 (vs +0.0380)
and disorder β = −0.0290 (vs −0.0294). PEST is the weakest term in the model
(β = −0.0028, p = 0.65), so it does not move the result.

> **Use the generator for new species, not to regenerate the cyno reference.**
> The snapshots in `data/snapshots/` stay authoritative — the published values are
> pinned to them. For a new species there is no table to match, so what matters is
> that the method is deterministic, documented and applied identically everywhere.

You can still supply any table yourself; an uploaded file overrides the generated
one. A tier whose table is missing is reported as unavailable, never approximated.

### Conventions

- **Half-life is in HOURS** in every input.
- **Gene symbol is the join key**; CORUM mapping uses the **uppercase** symbol,
  UniProt accession is secondary.
- **Composite gene names** use `;` in the cohort and `_` in FASTA/DIAMOND. Both
  spellings are matched automatically.
- Rows resolving to the **same accession** are deduplicated, keeping the
  unambiguous (non-composite) row.

---

## Citation

See `CITATION.cff`. Please cite both the paper and the pinned release of this
repository (tag a release and archive it via Zenodo so it carries a DOI).

## License

Code: **BSD 3-Clause** — see `LICENSE`. Note clause 3: the name of the copyright
holder and of contributors may not be used to endorse or promote derived products
without prior written permission.

Third-party data retains its original terms — see **Data provenance and
licensing**.
