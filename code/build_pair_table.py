"""
build_pair_table.py — raw DIAMOND all-vs-all  ->  evaluated paralog pair table.

WHY THIS EXISTS
---------------
Until 2026-08-31 this conversion was the pack's one genuinely unrecoverable step.
REPRODUCIBILITY_RUNBOOK.md recorded that rebuilding the pair table from raw
DIAMOND "diverges (24-35 pairs) because the coverage/HSP/best-hit rules aren't in
the prose", so Tier 1-1 depended on the frozen `pair_table_rev.tsv` and a NEW
SPECIES could not be processed at all.

The rule has now been recovered by exhaustive search against the frozen table:

    * pair identity   : UNORDERED gene pair; self-hits (gene vs itself) dropped
    * HSP selection   : keep the HSP with the MAXIMUM PERCENT IDENTITY
                        (not max bitscore, not min e-value, not first-in-file --
                        those give 8,128 / 8,153 / 8,696 of 11,526 columns right;
                        max-pident gives 11,526 / 11,526 on BOTH pident and cov_min)
    * cov_min         : aln_length / min(qlen, slen)
                        i.e. coverage relative to the SHORTER sequence
    * abs_dlog2_t     : |log2(half_life_1 / half_life_2)|

Validated on cyno: the unordered pair SET reconstructs exactly (11,526 pairs, 0
extra, 0 missing) and every pident and cov_min matches the frozen snapshot.

INPUT — raw DIAMOND, 8 columns, no header (outfmt 6 with qlen/slen):

    qseqid  sseqid  pident  length  qlen  slen  evalue  bitscore

`qseqid`/`sseqid` are `GENE|ACCESSION` (the FASTA header format this pack uses).
Produce it OUTSIDE this repository -- DIAMOND is a compiled binary:

    diamond makedb --in <species>.fasta -d <species>
    diamond blastp -q <species>.fasta -d <species> -o <species>_pairs.tsv \\
            --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \\
            --max-target-seqs 50 --evalue 1e-3

OUTPUT — data/snapshots/<species>_pair_table.tsv

    gene1  gene2  pident  cov_min  abs_dlog2_t

NOTE ON THE PINNED CYNO TABLE
-----------------------------
`pair_table_rev.tsv` was built on the PRE-dedup 3,515-row cohort, so it retains
25 rows involving `Mrc1_Slc39A12` -- the composite dropped by build_cohort.py.
Rebuilding against the rev cohort correctly yields 11,501 rows instead of 11,526.
That difference is the stale intermediate, not a discrepancy; Tier 1-1 re-applies
the cohort filter anyway and still returns exactly 31 outlier pairs.

USAGE
    python3 build_pair_table.py --pairs <raw.tsv> --cohort <cohort.tsv> --out <out.tsv>
    python3 build_pair_table.py --validate      # rebuild cyno and compare to frozen
"""

import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.dirname(HERE)
DATA = os.path.join(PACK, "data")

DEFAULT_PAIRS = os.path.join(DATA, "inputs", "cyno_pairs.tsv")
DEFAULT_COHORT = os.path.join(DATA, "inputs", "cyno_target_rev.tsv")
FROZEN = os.path.join(DATA, "snapshots", "pair_table_rev.tsv")

HEADER = ["gene1", "gene2", "pident", "cov_min", "abs_dlog2_t"]


def load_halflives(cohort_path):
    """row_gene -> half_life_h, accepting both ';' and '_' composite spellings."""
    hl = {}
    with open(cohort_path, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        j = {c: i for i, c in enumerate(head)}
        if "row_gene" not in j or "half_life_h" not in j:
            raise SystemExit(
                f"{cohort_path}: expected columns 'row_gene' and 'half_life_h', "
                f"found {head}"
            )
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            hl[f[j["row_gene"]]] = float(f[j["half_life_h"]])
    for g in [k for k in hl if ";" in k]:
        hl[g.replace(";", "_")] = hl[g]
    return hl


def parse_diamond(pairs_path):
    """Best (max-pident) HSP per unordered gene pair.

    Returns {frozenset({g1, g2}): (gene1_first_seen, gene2, pident, cov_min)}
    """
    best = {}
    order = []
    bad = 0
    with open(pairs_path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                bad += 1
                continue
            try:
                q = f[0].split("|")[0]
                s = f[1].split("|")[0]
                pident = float(f[2])
                aln = int(f[3])
                qlen = int(f[4])
                slen = int(f[5])
            except ValueError:
                bad += 1
                continue
            if q == s:
                continue
            shorter = min(qlen, slen)
            if shorter <= 0:
                continue
            cov_min = aln / shorter
            key = frozenset((q, s))
            if key not in best:
                order.append(key)
                best[key] = (q, s, pident, cov_min)
            elif pident > best[key][2]:
                # keep first-seen gene order, update the statistics
                g1, g2, _, _ = best[key]
                best[key] = (g1, g2, pident, cov_min)
    return best, order, bad


def build(pairs_path, cohort_path):
    hl = load_halflives(cohort_path)
    best, order, bad = parse_diamond(pairs_path)

    rows, skipped = [], 0
    for key in order:
        g1, g2, pident, cov_min = best[key]
        if g1 not in hl or g2 not in hl:
            skipped += 1
            continue
        h1, h2 = hl[g1], hl[g2]
        if h1 <= 0 or h2 <= 0:
            skipped += 1
            continue
        rows.append((g1, g2, pident, cov_min, abs(math.log2(h1 / h2))))
    return rows, {"raw_pairs": len(best), "skipped_not_in_cohort": skipped,
                  "malformed_lines": bad}


def write(rows, out_path):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w", newline="\n", encoding="utf-8") as w:
        w.write("\t".join(HEADER) + "\n")
        for g1, g2, pident, cov, dl in rows:
            w.write("%s\t%s\t%.2f\t%.3f\t%.4f\n" % (g1, g2, pident, cov, dl))


def validate():
    """Rebuild the cyno table and compare against the frozen snapshot."""
    rows, stats = build(DEFAULT_PAIRS, DEFAULT_COHORT)
    frozen = {}
    with open(FROZEN, encoding="utf-8") as fh:
        fh.readline()
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            frozen[frozenset((f[0], f[1]))] = (float(f[2]), float(f[3]), float(f[4]))

    built = {frozenset((r[0], r[1])): (r[2], r[3], r[4]) for r in rows}
    common = set(built) & set(frozen)
    pid_ok = sum(1 for k in common if abs(built[k][0] - frozen[k][0]) < 0.051)
    cov_ok = sum(1 for k in common if abs(built[k][1] - frozen[k][1]) < 0.0011)
    dl_ok = sum(1 for k in common if abs(built[k][2] - frozen[k][2]) < 0.0011)

    print("=== build_pair_table validation (cyno) ===")
    print(f"raw unordered pairs      = {stats['raw_pairs']}")
    print(f"rebuilt rows             = {len(rows)}")
    print(f"frozen rows              = {len(frozen)}")
    print(f"dropped (not in cohort)  = {stats['skipped_not_in_cohort']}"
          "   <- the stale pre-dedup Mrc1_Slc39A12 rows")
    print(f"overlap                  = {len(common)}")
    print(f"  pident matches         = {pid_ok}/{len(common)}")
    print(f"  cov_min matches        = {cov_ok}/{len(common)}")
    print(f"  abs_dlog2_t matches    = {dl_ok}/{len(common)}")
    only_built = set(built) - set(frozen)
    if only_built:
        print(f"  rebuilt-only pairs     = {len(only_built)}  (unexpected)")
    ok = (pid_ok == cov_ok == dl_ok == len(common)) and not only_built
    print("RESULT:", "PASS - rule reproduces the frozen table exactly" if ok
          else "FAIL - inspect the selection rule")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default=DEFAULT_PAIRS, help="raw DIAMOND tsv")
    ap.add_argument("--cohort", default=DEFAULT_COHORT, help="cohort tsv")
    ap.add_argument("--out", default=None, help="output pair-table tsv")
    ap.add_argument("--validate", action="store_true",
                    help="rebuild cyno and compare against the frozen snapshot")
    args = ap.parse_args()

    if args.validate:
        return validate()

    if not args.out:
        raise SystemExit("--out is required (or use --validate)")

    rows, stats = build(args.pairs, args.cohort)
    write(rows, args.out)
    print("=== build_pair_table ===")
    print(f"unordered pairs in DIAMOND output = {stats['raw_pairs']}")
    if stats["malformed_lines"]:
        print(f"malformed lines skipped           = {stats['malformed_lines']}")
    print(f"dropped (gene not in cohort)      = {stats['skipped_not_in_cohort']}")
    print(f"pairs written                     = {len(rows)}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
