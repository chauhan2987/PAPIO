#!/usr/bin/env python3
"""
fusion_truth_match.py - compare pipeline fusion calls with a published fusion list
(e.g. Talsania et al. 2022, Additional file 9, column "RNA-seq fusion_genes").

Handles
  * gene order: A-B and B-A both count (orientation reported separately)
  * renamed genes: old symbols translated to current ones with the HGNC table (optional --hgnc)
  * intergenic sides written as lists: "BICC1-EEF1A1P2(60251),RP11-226P1.2(4655)"
  * hyphenated clone names: "SMG1P2-RP11-231C14.4" is SMG1P2 + RP11-231C14.4
  * several fusions in one cell: "SHQ1-HNRNPA3P8(27392),RP11-481N16.1(197455), ATXN7-ATXN7"
  * the same for the pipeline's own Gene5/Gene3 columns (comma lists)

Usage
  python3 fusion_truth_match.py --truth talsania_fusions.csv --calls rna_1_R_fusions_annotated.tsv \
         [--hgnc hgnc_complete_set.txt] [--column "RNA-seq fusion_genes"] [--out talsania_result.tsv]
Python 3 standard library only.
"""
import argparse, csv, re, sys

CLONE_PREFIX = re.compile(r"^(RP\d+|CT[A-D]|XXBAC|XXYAC|LL\d+|KB|WI2|CITF\d+)$", re.I)

def clean_gene(g):
    """strip '(distance)', spaces and quotes; upper case"""
    return re.sub(r"\(\d+\)", "", g).strip().strip('"').upper()

def split_side(side):
    """'A(12),B(34)' -> {'A','B'}"""
    return {clean_gene(x) for x in side.split(",") if clean_gene(x)}

def split_fusion(text):
    """'A(1),B(2)-C' -> ({A,B},{C}); chooses the '-' that is not inside a clone name like RP11-231C14.4"""
    text = text.strip()
    for m in re.finditer("-", text):
        left, right = text[:m.start()], text[m.end():]
        last_left = re.split(r"[,(]", left)[-1].strip()
        if CLONE_PREFIX.match(last_left):          # the '-' belongs to a clone name
            continue
        l, r = split_side(left), split_side(right)
        if l and r:
            return l, r
    return None

def parse_truth_cell(cell):
    """one cell may hold several fusions separated by ', ' (comma + space)"""
    fusions = []
    for part in re.split(r",\s+", cell.strip().strip('"')):
        f = split_fusion(part)
        if f:
            fusions.append(f)
    return fusions

def load_hgnc(path):
    """previous and alias symbols -> current approved symbol (ambiguous aliases are skipped)"""
    prev, alias, current = {}, {}, set()
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            sym = r["symbol"].upper(); current.add(sym)
            for col, store in (("prev_symbol", prev), ("alias_symbol", alias)):
                for old in (r.get(col) or "").strip('"').split("|"):
                    old = old.strip().upper()
                    if old:
                        store.setdefault(old, set()).add(sym)
    m = {}
    for store in (alias, prev):                    # previous symbols take priority over aliases
        for old, syms in store.items():
            if old not in current and len(syms) == 1:
                m[old] = next(iter(syms))
    return m

def norm(g, hgnc):
    return hgnc.get(g, g)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True, help="CSV with a fusion-gene column (Talsania Additional file 9)")
    ap.add_argument("--column", default="RNA-seq fusion_genes", help="name of the fusion-gene column")
    ap.add_argument("--calls", required=True, nargs="+", help="pipeline *_fusions_annotated.tsv (one or more)")
    ap.add_argument("--hgnc", help="hgnc_complete_set.txt, to translate renamed genes")
    ap.add_argument("--out", default="talsania_result.tsv")
    a = ap.parse_args()

    hgnc = load_hgnc(a.hgnc) if a.hgnc else {}

    # ---- truth
    truth = []
    with open(a.truth, newline="", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        if a.column not in rd.fieldnames:
            sys.exit(f"ERROR: column '{a.column}' not in {a.truth}; columns are: {rd.fieldnames}")
        for r in rd:
            cell = (r[a.column] or "").strip()
            fus = [({norm(g, hgnc) for g in l}, {norm(g, hgnc) for g in rr}) for l, rr in parse_truth_cell(cell)]
            if not fus:
                print(f"note: could not parse '{cell}' ({r.get('SV_id', '')})", file=sys.stderr)
            truth.append({"id": r.get("SV_id", ""), "type": r.get("SV_type", ""), "text": cell, "fusions": fus})

    # ---- calls
    calls = []
    for path in a.calls:
        with open(path) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#") and l.strip()]
        hdr = {h: i for i, h in enumerate(rows[0])}
        for k in ("Gene5", "Gene3"):
            if k not in hdr:
                sys.exit(f"ERROR: no {k} column in {path}; columns are: {rows[0]}")
        for f in rows[1:]:
            g5 = {norm(clean_gene(x), hgnc) for x in f[hdr["Gene5"]].split(",") if clean_gene(x)}
            g3 = {norm(clean_gene(x), hgnc) for x in f[hdr["Gene3"]].split(",") if clean_gene(x)}
            calls.append({"g5": g5, "g3": g3,
                          "verdict": f[hdr["Verdict"]] if "Verdict" in hdr else "",
                          "tier": f[hdr["TIER"]] if "TIER" in hdr else "",
                          "name": f[hdr["FusionName"]] if "FusionName" in hdr else ""})

    # ---- match: a truth fusion is found when one call joins a gene from each side (either orientation)
    found = 0
    with open(a.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["SV_id", "SV_type", "truth_fusion", "status", "matched_call", "orientation", "Verdict", "TIER",
                    "partner_genes_seen_in_calls"])
        seen_genes = set().union(*[c["g5"] | c["g3"] for c in calls]) if calls else set()
        for t in truth:
            hit = None
            for l, r in t["fusions"]:
                for c in calls:
                    same = (c["g5"] & l) and (c["g3"] & r)
                    swap = (c["g5"] & r) and (c["g3"] & l)
                    if same or swap:
                        hit = (c, "same" if same else "swapped"); break
                if hit: break
            genes = set().union(*[l | r for l, r in t["fusions"]]) if t["fusions"] else set()
            seen = ",".join(sorted(genes & seen_genes)) or "-"
            if hit:
                found += 1
                c, o = hit
                w.writerow([t["id"], t["type"], t["text"], "FOUND", c["name"], o, c["verdict"], c["tier"], seen])
            else:
                w.writerow([t["id"], t["type"], t["text"], "missed", "", "", "", "", seen])
    n = len(truth)
    print(f"truth fusions: {n}   found: {found}   recall: {found / n:.3f}" if n else "no truth rows")
    print(f"details: {a.out}")

if __name__ == "__main__":
    main()