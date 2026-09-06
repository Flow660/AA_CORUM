"""
annotate.py — generate the annotation tables that Tiers 1-2, 1-4, 1-5 and 1-6 need.

Before this module a new species could only reach Tier 1-1 / 1-1b, because the
other four tiers each require an annotation table that the user had to build by
hand. This automates all four.

    seq-features   FASTA + cohort                 -> Tier 1-4        [OFFLINE]
    repeat-class   UniProt domains + cohort       -> Tier 1-2        [network once]
    corum          CORUM file + cohort            -> Tier 1-5, 1-6   [offline]

VALIDATION STATUS (measured against the pinned cyno snapshots)
--------------------------------------------------------------
`seq-features` — recovered exactly, from sequence alone:

    length          3508/3508 exact
    pI              3507/3508 exact (mean abs error 0.00006)
                    Bjellqvist/Biopython pKa set WITH residue-specific N-/C-terminal
                    pKa; the plain EMBOSS or Sillero sets do NOT reproduce it
    disorder_frac   exact — it is the fraction of residues in {A,R,G,Q,S,P,E,K}
                    (a composition metric, NOT a disorder predictor)
    nend_residue    3508/3508 exact, applying Met-excision: if the sequence starts
                    with M and residue 2 is small (A,C,G,P,S,T,V), residue 2 is the
                    N-terminal residue
    nend_class      3508/3508 exact — D,E,Q,N -> NDEG; R,K,H -> TYPE1;
                    F,L,I,W,Y -> TYPE2; everything else -> STAB

    pest_max / pest_n_stretches  ** NOT byte-reproducible **. The original PEST
    implementation's region-splitting rule was not recovered (a standard
    PEST-FIND split at K/R/H reproduces only 278/3358 stretch counts). The
    implementation here is a documented, deterministic PEST-FIND and will NOT
    equal the cyno snapshot.
    This matters less than it sounds: in Tier 1-4, PEST is the weakest term in
    the model -- beta = -0.0028, p = 0.646, i.e. not significant. It does not
    move R^2 = 0.0190.

`repeat-class` — keyword mapping over UniProt domain/repeat annotations. The
cyno snapshot was built the same way, but the exact UniProt release is not
recorded, so counts may drift slightly with the current release.

`corum` — reproduces the published stratification for **3,409 / 3,514 genes
(97.0 %)**. The residual 105 genes come from an unrecovered detail of the
original gene->CORUM mapping. Class counts come out 2043/119/1196/156 against a
published 2072/127/1132/183.

    >>> USE THIS FOR NEW SPECIES, NOT TO REGENERATE THE CYNO REFERENCE. <<<
    The cyno snapshots in data/snapshots/ remain authoritative and must not be
    overwritten -- the published values are pinned to them.

For a NEW species there is no published table to match, so what matters is that
the method is deterministic, documented and applied identically to every species.
That this module reproduces the cyno reference to 97-100 % is evidence the method
is the right one, not a target to hit.

USAGE
    python3 annotate.py seq-features --fasta <f> --cohort <c> --out <o>
    python3 annotate.py repeat-class --cohort <c> --out <o> [--uniprot-tsv <f>]
    python3 annotate.py corum --corum <corum_humanComplexes.txt> --cohort <c> \\
                              --out-member <m> --out-stats <s> --out-strat <t>
    python3 annotate.py validate --fasta <f> --cohort <c>   # compare to the snapshot
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
DATA = os.path.join(PACK, "data")

# --------------------------------------------------------------------- constants

# Bjellqvist pKa set as used by Biopython's IsoelectricPoint.
POS_PK = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
NEG_PK = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
PK_NTERM = {"A": 7.59, "M": 7.00, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44,
            "E": 7.70}
PK_CTERM = {"D": 4.55, "E": 4.75}

# Fraction of these residues == the snapshot's `disorder_frac` (verified exact).
DISORDER_SET = set("ARGQSPEK")

# Met-excision: cleaved when residue 2 is small.
SMALL = set("ACGPSTV")

NEND_CLASS = {}
NEND_CLASS.update({r: "NDEG" for r in "DEQN"})
NEND_CLASS.update({r: "TYPE1" for r in "RKH"})
NEND_CLASS.update({r: "TYPE2" for r in "FLIWY"})

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

NO_PEST = -1e9

REPEAT_PATTERNS = [
    ("TPR", [r"\bTPR\b", r"tetratricopeptide"]),
    ("ANK", [r"\bANK\b", r"ankyrin"]),
    ("WD40", [r"\bWD\b", r"WD40", r"WD repeat"]),
    ("LRR", [r"\bLRR\b", r"leucine-rich repeat"]),
    ("HEAT_ARM", [r"\bHEAT\b", r"armadillo", r"\bARM\b"]),
    ("KELCH_FBX", [r"kelch", r"F-box"]),
    ("SPECTRIN", [r"spectrin"]),
]

# CORUM complex-name patterns marking a STRUCTURAL / stoichiometric machine.
# Derived by constraint propagation from the pinned cyno stratification: a gene
# labelled STOICH belongs only to STOICH complexes, a REGUL gene only to REGUL
# ones, which uniquely labels 2,677 complexes and yields these 57 distinct names.
STOICH_PATTERNS = [
    r"proteasome", r"ribosom", r"spliceosome", r"snRNP", r"APC/C",
    r"anaphase-promoting", r"\bBAF complex", r"PBAF", r"SWI/SNF", r"NUMAC",
    r"CCR4-NOT", r"chaperonin", r"\bCCT complex", r"ESCRT", r"\bexosome",
    r"F1F0-ATPase", r"v-ATPase", r"GATOR", r"mediator complex", r"nuclear pore",
    r"RNA polymerase", r"ragulator", r"TFII", r"HPF1-PARP",
    r"MeCP2-SIN3A-HDAC", r"TIP60", r"TREX complex", r"exon junction complex",
    r"\bTRAP-SMCC",
]
STOICH_RX = [re.compile(p, re.I) for p in STOICH_PATTERNS]


# ------------------------------------------------------------------------- io

def read_fasta(path):
    seqs, gene = {}, None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                gene = line[1:].split("|")[0]
                seqs[gene] = []
            elif gene is not None:
                seqs[gene].append(line)
    out = {}
    for g, parts in seqs.items():
        s = "".join(parts).upper()
        out[g] = s
        out[g.replace("_", ";")] = s
    return out


def read_cohort(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        j = {c: i for i, c in enumerate(head)}
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            rows.append({
                "gene": f[j["row_gene"]],
                "alias": f[j["chosen_alias"]],
                "accession": f[j["accession"]],
                "length": int(f[j["length"]]),
                "half_life_h": float(f[j["half_life_h"]]),
            })
    return rows


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="\n", encoding="utf-8") as w:
        w.write("\t".join(header) + "\n")
        for r in rows:
            w.write("\t".join(str(x) for x in r) + "\n")


# ------------------------------------------------------------- sequence features

def isoelectric_point(seq):
    cnt = Counter(seq)
    nt = PK_NTERM.get(seq[0], POS_PK["Nterm"])
    ct = PK_CTERM.get(seq[-1], NEG_PK["Cterm"])

    def charge(ph):
        c = 1.0 / (1 + 10 ** (ph - nt)) - 1.0 / (1 + 10 ** (ct - ph))
        for aa in "KRH":
            c += cnt.get(aa, 0) / (1 + 10 ** (ph - POS_PK[aa]))
        for aa in "DECY":
            c -= cnt.get(aa, 0) / (1 + 10 ** (NEG_PK[aa] - ph))
        return c

    lo, hi = 0.0, 14.0
    for _ in range(120):
        mid = (lo + hi) / 2.0
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def nend(seq):
    if seq.startswith("M") and len(seq) > 1 and seq[1] in SMALL:
        return seq[1]
    return seq[0]


def disorder_fraction(seq):
    return sum(1 for c in seq if c in DISORDER_SET) / len(seq)


def pest(seq):
    """Deterministic PEST-FIND. See the module docstring: this does NOT
    reproduce the cyno snapshot's pest_max, and PEST is not significant in the
    Tier 1-4 model (beta = -0.0028, p = 0.646)."""
    idx = [-1] + [i for i, c in enumerate(seq) if c in "KRH"] + [len(seq)]
    best, n = NO_PEST, 0
    for a, b in zip(idx, idx[1:]):
        r = seq[a + 1:b]
        if len(r) < 12:
            continue
        if "P" not in r or not any(c in r for c in "DE") or not any(c in r for c in "ST"):
            continue
        n += 1
        L = len(r)
        depst = sum(r.count(c) for c in "DEPST") / L * 100.0
        hyd = sum(r.count(a2) / L * 100.0 * KD[a2] for a2 in KD) / 100.0
        score = 0.55 * depst - 0.5 * hyd
        best = max(best, score)
    return best, n


def cmd_seq_features(args):
    seqs = read_fasta(args.fasta)
    cohort = read_cohort(args.cohort)
    repeat = {}
    if args.repeat_class and os.path.exists(args.repeat_class):
        with open(args.repeat_class, encoding="utf-8") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            j = {c: i for i, c in enumerate(head)}
            for ln in fh:
                if not ln.strip():
                    continue
                f = ln.rstrip("\n").split("\t")
                repeat[f[j["row_gene"]]] = f[j["repeat_class"]]

    rows, missing = [], 0
    for c in cohort:
        s = seqs.get(c["gene"]) or seqs.get(c["gene"].replace(";", "_"))
        if not s:
            missing += 1
            continue
        r = nend(s)
        pmax, pn = pest(s)
        rows.append([
            c["gene"], c["alias"], c["accession"], len(s),
            "%.4f" % c["half_life_h"],
            "%.4f" % isoelectric_point(s),
            "%.4f" % disorder_fraction(s),
            r, NEND_CLASS.get(r, "STAB"),
            "%.4f" % pmax, pn,
            repeat.get(c["gene"], "OTHER"),
        ])
    write_tsv(args.out, ["gene", "alias", "accession", "length", "half_life_h",
                         "pI", "disorder_frac", "nend_residue", "nend_class",
                         "pest_max", "pest_n_stretches", "repeat_class"], rows)
    print("=== annotate seq-features ===")
    print(f"cohort rows           = {len(cohort)}")
    print(f"written               = {len(rows)}")
    if missing:
        print(f"no sequence in FASTA  = {missing}")
    print(f"N-end classes         = {dict(Counter(r[8] for r in rows))}")
    print("NOTE: pest_max uses a fresh PEST-FIND implementation and will not "
          "equal the\n      cyno snapshot. PEST is not significant in Tier 1-4 "
          "(p = 0.646).")
    print(f"wrote {args.out}")
    return 0


# ------------------------------------------------------------------ repeat class

def classify_domains(text):
    for label, pats in REPEAT_PATTERNS:
        for p in pats:
            if re.search(p, text, re.I):
                return label
    return "OTHER"


def fetch_uniprot(accessions, batch=60, retries=3, progress=None):
    """Fetch domain/repeat annotation from UniProt REST. Requires network.

    Hardened for real-world use: a cohort of ~3,500 proteins means dozens of
    sequential requests, and a single transient failure (corporate proxy, rate
    limit, dropped connection) must not throw the whole run away. Each batch is
    retried with exponential backoff, and batches that still fail are reported
    rather than raising -- those proteins simply fall back to 'OTHER'.

    Returns (annotations, stats) where stats records coverage and failures.
    """
    import time
    import urllib.error
    import urllib.parse
    import urllib.request

    out = {}
    accs = sorted(accessions)
    failed_batches = 0
    consecutive_failures = 0
    n_batches = (len(accs) + batch - 1) // batch

    for bi, i in enumerate(range(0, len(accs), batch), start=1):
        chunk = accs[i:i + batch]
        query = " OR ".join("accession:%s" % a for a in chunk)
        url = ("https://rest.uniprot.org/uniprotkb/stream?"
               + urllib.parse.urlencode({
                   "query": query,
                   "format": "tsv",
                   "fields": "accession,ft_repeat,ft_domain,cc_domain,protein_families",
               }))
        req = urllib.request.Request(url, headers={
            "User-Agent": "cyno-halflife-pipeline (annotate.py)",
            "Accept": "text/plain",
        })

        text = None
        for attempt in range(1, retries + 1):
            try:
                with urllib.request.urlopen(req, timeout=180) as resp:
                    text = resp.read().decode("utf-8", "replace")
                break
            except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
                if attempt == retries:
                    failed_batches += 1
                    msg = f"  batch {bi}/{n_batches} FAILED after {retries} tries: {exc}"
                    sys.stderr.write(msg + "\n")
                else:
                    time.sleep(2 ** attempt)          # 2s, 4s

        if text is None:
            consecutive_failures += 1
            # Fail fast: if the first few batches all fail the network is down or
            # blocked, and grinding through dozens more just wastes the user's time.
            if consecutive_failures >= 3 and not out:
                raise RuntimeError(
                    f"UniProt is unreachable - the first {consecutive_failures} "
                    "batches all failed.\n"
                    "  Check your internet connection or corporate proxy, or "
                    "export the annotation\n"
                    "  from https://www.uniprot.org/id-mapping and pass "
                    "--uniprot-tsv (upload it in the app).")
        else:
            consecutive_failures = 0
            for ln in text.splitlines()[1:]:
                f = ln.split("\t")
                if f and f[0]:
                    out[f[0]] = " ".join(f[1:])

        done = min(i + batch, len(accs))
        if progress:
            progress(done, len(accs))
        else:
            sys.stderr.write(f"  fetched {done}/{len(accs)} accessions\n")

    stats = {"requested": len(accs), "annotated": len(out),
             "failed_batches": failed_batches, "batches": n_batches}
    return out, stats


def cmd_repeat_class(args):
    cohort = read_cohort(args.cohort)
    if args.uniprot_tsv:
        ann = {}
        with open(args.uniprot_tsv, encoding="utf-8", errors="replace") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if f:
                    ann[f[0]] = " ".join(f[1:])
        print(f"loaded {len(ann)} annotations from {args.uniprot_tsv}")
    else:
        print("fetching domain annotation from UniProt (network required)...")
        ann, fstats = fetch_uniprot({c["accession"] for c in cohort})
        print(f"  accessions requested = {fstats['requested']}")
        print(f"  annotations returned = {fstats['annotated']}")
        if fstats["failed_batches"]:
            print(f"  WARNING: {fstats['failed_batches']}/{fstats['batches']} "
                  "batches failed; those proteins default to OTHER.")
            print("  Re-run, or export the TSV from the UniProt website and pass "
                  "--uniprot-tsv.")

    covered = sum(1 for c in cohort if ann.get(c["accession"]))
    if covered == 0:
        raise RuntimeError(
            "No UniProt annotation was retrieved for any protein, so every class "
            "would be 'OTHER' and Tier 1-2 would be meaningless. Refusing to write "
            "a misleading file.\n"
            "  * If this was a network/proxy failure: retry, or export the TSV from "
            "https://www.uniprot.org/id-mapping and pass --uniprot-tsv.\n"
            "  * If it was not: check that your cohort's 'accession' column holds "
            "real UniProt accessions.")

    rows = []
    for c in cohort:
        rows.append([c["gene"], c["alias"], c["accession"], c["length"],
                     repr(c["half_life_h"]),
                     classify_domains(ann.get(c["accession"], ""))])
    write_tsv(args.out, ["row_gene", "alias", "accession", "length",
                         "half_life_h", "repeat_class"], rows)
    counts = Counter(r[5] for r in rows)
    print("=== annotate repeat-class ===")
    print(f"proteins            = {len(rows)}")
    print(f"with UniProt record = {covered} "
          f"({100.0 * covered / max(len(rows), 1):.1f}%)")
    print(f"classes = {dict(counts)}")
    if covered < len(rows) * 0.8:
        print("WARNING: under 80% of proteins matched a UniProt record. Check that "
              "the\n         accession column of your cohort holds real UniProt "
              "accessions.")
    print(f"wrote {args.out}")
    return 0


# ------------------------------------------------------------------------ corum

def cmd_corum(args):
    cohort = read_cohort(args.cohort)
    by_symbol = {}
    for c in cohort:
        by_symbol.setdefault(c["alias"].upper(), c)

    complexes = []
    with open(args.corum, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if args.organism and r.get("organism") != args.organism:
                continue
            genes = {x.strip().upper()
                     for x in (r.get("subunits_gene_name") or "").split(";")
                     if x.strip()}
            if genes:
                complexes.append((r.get("complex_id"), r.get("complex_name", ""), genes))
    print(f"complexes loaded = {len(complexes)}")

    n_cx = Counter()
    classes = {}
    for cid, name, genes in complexes:
        stoich = any(rx.search(name) for rx in STOICH_RX)
        for g in genes:
            if g in by_symbol:
                n_cx[g] += 1
                classes.setdefault(g, set()).add("STOICH" if stoich else "REGUL")

    member_rows, strat_rows = [], []
    for c in cohort:
        sym = c["alias"].upper()
        n = n_cx.get(sym, 0)
        member_rows.append([sym, c["accession"], "%.4f" % c["half_life_h"],
                            1 if n else 0, n])
        s = classes.get(sym)
        cls = ("NONE" if not s else
               "STOICH" if s == {"STOICH"} else
               "REGUL" if s == {"REGUL"} else "BOTH")
        strat_rows.append([sym, "%.4f" % c["half_life_h"], n, cls])

    stat_rows = []
    for cid, name, genes in complexes:
        measured = [by_symbol[g] for g in genes if g in by_symbol]
        if len(measured) < 2:
            continue
        hl = sorted(m["half_life_h"] for m in measured)
        mid = len(hl) // 2
        median = hl[mid] if len(hl) % 2 else (hl[mid - 1] + hl[mid]) / 2
        lg = [math.log2(h) for h in hl]
        mean = sum(lg) / len(lg)
        sd = math.sqrt(sum((x - mean) ** 2 for x in lg) / (len(lg) - 1))
        stat_rows.append([cid, name, len(measured), "%.4f" % median, "%.4f" % sd,
                          ";".join(sorted(m["alias"].upper() for m in measured))])

    write_tsv(args.out_member, ["gene", "accession", "half_life_h", "in_corum",
                                "n_complexes"], member_rows)
    write_tsv(args.out_stats, ["complex_id", "complex_name", "n_measured",
                               "median_halflife_h", "sd_log2_halflife",
                               "measured_genes"], stat_rows)
    write_tsv(args.out_strat, ["gene", "half_life_h", "n_complexes",
                               "corum_class"], strat_rows)
    print("=== annotate corum ===")
    print(f"members / non-members = {sum(1 for r in member_rows if r[3])} / "
          f"{sum(1 for r in member_rows if not r[3])}")
    print(f"complexes (n>=2)      = {len(stat_rows)}")
    print(f"classes               = {dict(Counter(r[3] for r in strat_rows))}")
    print("NOTE: reproduces the pinned cyno stratification for ~97% of genes. "
          "Use for NEW\n      species; do not overwrite the authoritative cyno "
          "snapshots.")
    print(f"wrote {args.out_member}\n      {args.out_stats}\n      {args.out_strat}")
    return 0


# --------------------------------------------------------------------- validate

def cmd_validate(args):
    seqs = read_fasta(args.fasta)
    snap_path = os.path.join(DATA, "snapshots", "seq_features_rev.tsv")
    snap = {}
    with open(snap_path, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        j = {c: i for i, c in enumerate(head)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            snap[f[0]] = f
    common = [g for g in snap if g in seqs]
    pi_ok = nr_ok = nc_ok = df_ok = ln_ok = 0
    for g in common:
        s = seqs[g]
        if len(s) == int(snap[g][j["length"]]):
            ln_ok += 1
        if abs(isoelectric_point(s) - float(snap[g][j["pI"]])) < 0.005:
            pi_ok += 1
        if abs(disorder_fraction(s) - float(snap[g][j["disorder_frac"]])) < 0.0001:
            df_ok += 1
        r = nend(s)
        if r == snap[g][j["nend_residue"]]:
            nr_ok += 1
        if NEND_CLASS.get(r, "STAB") == snap[g][j["nend_class"]]:
            nc_ok += 1
    n = len(common)
    print("=== annotate validate (against data/snapshots/seq_features_rev.tsv) ===")
    print(f"proteins compared = {n}")
    print(f"  length          {ln_ok}/{n}")
    print(f"  pI              {pi_ok}/{n}")
    print(f"  disorder_frac   {df_ok}/{n}")
    print(f"  nend_residue    {nr_ok}/{n}")
    print(f"  nend_class      {nc_ok}/{n}")
    print("  pest_max        not comparable by design (see module docstring)")
    return 0 if min(ln_ok, pi_ok, df_ok, nr_ok, nc_ok) >= n - 2 else 1


# ------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("seq-features")
    p.add_argument("--fasta", required=True)
    p.add_argument("--cohort", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--repeat-class", default=None,
                   help="optional repeat-class TSV to merge in")
    p.set_defaults(func=cmd_seq_features)

    p = sub.add_parser("repeat-class")
    p.add_argument("--cohort", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--uniprot-tsv", default=None,
                   help="pre-downloaded UniProt TSV (accession + domain columns); "
                        "omit to fetch live")
    p.set_defaults(func=cmd_repeat_class)

    p = sub.add_parser("corum")
    p.add_argument("--corum", required=True)
    p.add_argument("--cohort", required=True)
    p.add_argument("--out-member", required=True)
    p.add_argument("--out-stats", required=True)
    p.add_argument("--out-strat", required=True)
    p.add_argument("--organism", default="Human")
    p.set_defaults(func=cmd_corum)

    p = sub.add_parser("validate")
    p.add_argument("--fasta", default=os.path.join(DATA, "inputs",
                                                   "cyno_target_rev.fasta"))
    p.add_argument("--cohort", default=os.path.join(DATA, "inputs",
                                                    "cyno_target_rev.tsv"))
    p.set_defaults(func=cmd_validate)

    args = ap.parse_args()
    try:
        return args.func(args)
    except RuntimeError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
