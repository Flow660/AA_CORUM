"""
build_cohort.py — the upstream data-prep step (Monkey_half-life.xlsx -> cyno_target_rev.tsv)

This closes the last unscripted stage of the pipeline. Until now the runbook said
the cohort was "derived from Monkey_half-life.xlsx (3,515 -> 3,514)" but never
recorded HOW, so the upstream step was not reproducible.

WHAT THE RAW WORKBOOK CONTAINS
------------------------------
Sheet1, 3 columns, 3,515 data rows:
    Gene Name  -- mouse-style gene symbol; may be a semicolon-joined composite
                  when the peptide evidence could not separate several genes
    Entry      -- MOUSE UniProt accession (e.g. Aak1 -> Q3UHJ0). NOT used in the
                  cohort: the cohort is keyed on the CYNO ortholog. Empty for
                  every composite row.
    Monkey     -- half-life in HOURS (the pipeline's internal unit)

WHAT THE COHORT NEEDS
---------------------
    row_gene, chosen_alias, accession, length, half_life_h

`accession` and `length` are CYNO values and originally came from a live UniProt
ortholog lookup. That lookup is frozen here as `data/inputs/cyno_target_rev.fasta`
(headers `>gene|accession`, sequence length = `length`), exactly the same
"snapshot an external database" pattern the rest of the pack uses. No network.

THE DEDUP RULE (3,515 -> 3,514) -- recovered, and now enforced in code
----------------------------------------------------------------------
The workbook holds SEVEN semicolon composites, and only ONE is dropped, so the
rule is emphatically NOT "drop composite gene names":

    H2Ac12;H2Ac7;H2Ac14;H2Ac17;H2Aj   kept  -> H2AC12
    H2Bc17;H2Bc11;H2Bc21              kept  -> H2BC17
    H2Bc5;H2Bc18;H2Bc15;H2Bc12L       kept  -> H2BC5
    Magohb;Magoh                      kept  -> MAGOHB
    Rack1;Gnb2L1                      kept  -> RACK1
    Snrpb;Snrpn                       kept  -> SNRPB
    Mrc1;Slc39A12                     DROPPED

The runbook's phrase is "excluded by accession", and that is literally the rule:

    Mrc1           -> A0A2K5TN91   (half-life  5.996546499 h)
    Mrc1;Slc39A12  -> A0A2K5TN91   (half-life 15.40101596  h)   <-- same protein

The workbook carries `Mrc1` as its own confident row, so the composite is a
REDUNDANT second measurement of the very same cyno protein. The six kept
composites each resolve to an accession no other row claims, so nothing collides
and they survive.

Independent corroboration from a frozen snapshot: `data/snapshots/pair_table_rev.tsv`
was built on the pre-dedup 3,515-row cohort and contains

    Mrc1    Mrc1_Slc39A12    pident=100.00    cov_min=1.000

i.e. a 100 % identity / full-coverage SELF-PAIR -- two rows, one protein. This is
the stale self-pair the runbook blames for the Tier 1-1 divergence.

    RULE: resolve every row to its cyno accession; where several rows share one
          accession, keep the CONFIDENT (non-composite) gene assignment and drop
          the ambiguous composite(s).

Because the frozen FASTA is post-dedup, it cannot by itself demonstrate the
collision. PREDEDUP_ACCESSIONS below pins the composite's resolution (recovered
from the self-pair above) so the rule actually FIRES and is auditable, rather
than being silently bypassed by a missing key.

OUTPUT FORMAT (byte-exact against data/inputs/cyno_target_rev.tsv)
------------------------------------------------------------------
    - row order      : workbook order (already alphabetical), dropped rows removed
    - chosen_alias   : first gene before ';', uppercased
    - half_life_h    : Python repr() of the float -- full precision, no rounding
    - line endings   : LF, with a trailing newline

DEPENDENCIES: standard library only. The .xlsx is parsed directly with zipfile +
xml.etree so the step needs neither openpyxl nor pandas, consistent with
ENVIRONMENT.md ("NumPy only", SciPy absent, no pip install).

USAGE
    python3 build_cohort.py                    # rebuild + verify against pinned
    python3 build_cohort.py --xlsx <f> --fasta <f> --out <f>
    python3 build_cohort.py --no-verify        # new species: skip the comparison
"""

import argparse
import os
import sys
import zipfile
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

DEFAULT_XLSX = os.path.join(DATA, "inputs", "Monkey_half-life.xlsx")
DEFAULT_FASTA = os.path.join(DATA, "inputs", "cyno_target_rev.fasta")
DEFAULT_OUT = os.path.join(HERE, "..", "results", "cyno_target_rev.tsv")
PINNED = os.path.join(DATA, "inputs", "cyno_target_rev.tsv")

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

# Pre-dedup cyno accession for composite rows that the post-dedup FASTA no longer
# carries. Recovered from the Mrc1 / Mrc1_Slc39A12 self-pair (pident 100.00,
# cov_min 1.000) in data/snapshots/pair_table_rev.tsv. Keeps the collision visible.
PREDEDUP_ACCESSIONS = {
    "Mrc1;Slc39A12": "A0A2K5TN91",
}

HEADER = ["row_gene", "chosen_alias", "accession", "length", "half_life_h"]


# ---------------------------------------------------------------- xlsx reading

def _col_index(ref):
    """'B7' -> 1 (zero-based column)."""
    n = 0
    for ch in ref:
        if ch.isalpha():
            n = n * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return n - 1


def read_xlsx(path):
    """Return [(gene, entry, half_life_hours), ...] from the first worksheet.

    Uses zipfile + ElementTree so no third-party reader is required.
    """
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{NS}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))

        root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        rows = []
        for r in root.iter(f"{NS}row"):
            cells = {}
            for c in r.findall(f"{NS}c"):
                ref = c.get("r") or ""
                idx = _col_index(ref)
                ctype = c.get("t")
                if ctype == "inlineStr":
                    v = "".join(t.text or "" for t in c.iter(f"{NS}t"))
                else:
                    ve = c.find(f"{NS}v")
                    if ve is None or ve.text is None:
                        continue
                    v = ve.text
                    if ctype == "s":
                        v = shared[int(v)]
                cells[idx] = v
            rows.append(cells)

    if not rows:
        raise SystemExit(f"no rows parsed from {path}")

    out = []
    for cells in rows[1:]:                       # drop the header row
        gene = (cells.get(0) or "").strip()
        if not gene:
            continue
        entry = (cells.get(1) or "").strip() or None
        hl = cells.get(2)
        if hl is None or hl == "":
            continue
        out.append((gene, entry, float(hl)))
    return out


# --------------------------------------------------------------- fasta reading

def read_fasta_map(path):
    """'>gene|accession' -> {gene_with_semicolons: (accession, seq_length)}.

    FASTA headers join composite genes with '_'; the cohort uses ';'. The map is
    returned keyed BOTH ways so either spelling resolves.
    """
    acc, seqlen, gene = {}, {}, None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                head = line[1:]
                gene, a = head.split("|", 1)
                acc[gene] = a.strip()
                seqlen[gene] = 0
            elif gene is not None:
                seqlen[gene] += len(line)

    out = {}
    for g, a in acc.items():
        out[g] = (a, seqlen[g])
        out[g.replace("_", ";")] = (a, seqlen[g])
    return out


# ------------------------------------------------------------------- the build

def build(xlsx, fasta):
    raw = read_xlsx(xlsx)
    fmap = read_fasta_map(fasta)

    resolved, unresolved = [], []
    for gene, _entry, hl in raw:
        hit = fmap.get(gene) or fmap.get(gene.replace(";", "_"))
        if hit is None:
            a = PREDEDUP_ACCESSIONS.get(gene)
            hit = (a, None) if a else None
        if hit is None:
            unresolved.append(gene)
            continue
        accession, length = hit
        resolved.append({
            "gene": gene,
            "alias": gene.split(";")[0].upper(),
            "accession": accession,
            "length": length,
            "half_life_h": hl,
            "composite": ";" in gene,
        })

    # --- dedup by accession: keep the confident (non-composite) assignment
    by_acc = {}
    for r in resolved:
        by_acc.setdefault(r["accession"], []).append(r)

    dropped = []
    for accession, group in by_acc.items():
        if len(group) < 2:
            continue
        confident = [r for r in group if not r["composite"]]
        if not confident:
            # no unambiguous row: keep the first in workbook order, drop the rest
            confident = group[:1]
        keep = confident[0]
        for r in group:
            if r is not keep:
                r["_drop"] = True
                dropped.append((r["gene"], accession, keep["gene"]))

    kept = [r for r in resolved if not r.get("_drop")]

    missing_len = [r["gene"] for r in kept if r["length"] is None]
    if missing_len:
        raise SystemExit(
            "kept rows with no sequence length (not in the FASTA snapshot): "
            + ", ".join(missing_len)
        )

    return kept, dropped, unresolved, len(raw)


def write_tsv(rows, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for r in rows:
            fh.write("\t".join([
                r["gene"], r["alias"], r["accession"],
                str(r["length"]), repr(r["half_life_h"]),
            ]) + "\n")


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", default=DEFAULT_XLSX)
    ap.add_argument("--fasta", default=DEFAULT_FASTA)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--pinned", default=PINNED)
    ap.add_argument("--no-verify", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.xlsx):
        raise SystemExit(
            f"workbook not found: {args.xlsx}\n"
            "Place Monkey_half-life.xlsx in data/inputs/ (it is a pinned raw input)."
        )

    kept, dropped, unresolved, n_raw = build(args.xlsx, args.fasta)

    print("=== build_cohort (Monkey_half-life.xlsx -> cohort tsv) ===")
    print(f"workbook rows      = {n_raw}")
    print(f"cohort rows        = {len(kept)}   (expected 3514)")

    print(f"dropped (duplicate accession) = {len(dropped)}")
    for gene, accession, keeper in dropped:
        print(f"    drop {gene!r}  accession={accession}  "
              f"-- redundant with confident row {keeper!r}")

    if unresolved:
        print(f"unresolved genes (no accession in the FASTA snapshot) = {len(unresolved)}")
        for g in unresolved[:20]:
            print(f"    {g}")

    composites = [r["gene"] for r in kept if r["composite"]]
    print(f"composite gene rows KEPT = {len(composites)}")
    for g in composites:
        print(f"    keep {g}")

    write_tsv(kept, args.out)
    print(f"\nwrote {args.out}")

    if args.no_verify:
        return 0

    if not os.path.exists(args.pinned):
        print(f"NOTE: no pinned file at {args.pinned}; skipping verification.")
        return 0

    with open(args.out, "rb") as a, open(args.pinned, "rb") as b:
        got, want = a.read(), b.read()

    if got == want:
        print("VERIFY: byte-identical to the pinned data/inputs/cyno_target_rev.tsv  [PASS]")
        return 0

    print("VERIFY: MISMATCH against the pinned cohort  [FAIL]")
    gl, wl = got.decode().split("\n"), want.decode().split("\n")
    print(f"    rebuilt lines={len(gl)}  pinned lines={len(wl)}")
    shown = 0
    for i, (x, y) in enumerate(zip(gl, wl)):
        if x != y:
            print(f"    line {i + 1}:\n        rebuilt {x}\n        pinned  {y}")
            shown += 1
            if shown >= 5:
                break
    return 1


if __name__ == "__main__":
    sys.exit(main())
