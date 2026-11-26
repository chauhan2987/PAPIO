#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Merge fusion callers into a consensus table.
#          Called by the PAPIO process FUSION_CONSENSUS; moved from an inline heredoc.
# Usage  : fusion_consensus.py <min_callers> <readthrough_kb> <known_arg> <sample_id> [genes_main.tsv] [genes_ctat.tsv]
#          genes_*.tsv: gene_id<TAB>gene_name from the RNA arm's GTF (main) and the CTAT library's GTF.
#          Callers are matched on Ensembl gene IDs (version removed), so Arriba (main GTF) and
#          STAR-Fusion / FusionInspector (CTAT GTF) agree even when a gene symbol differs between
#          GENCODE releases. Names are reported with the main GTF's symbols; the FusionInspector
#          candidate list uses the CTAT symbols, which FusionInspector looks up.
#------------------------------------------------------------------------------
import sys, csv, re

MIN_CALLERS = int(sys.argv[1])
READTHRU_BP = int(sys.argv[2]) * 1000
KNOWN_PATH  = sys.argv[3]
SAMPLE      = sys.argv[4]
MAIN_MAP    = sys.argv[5] if len(sys.argv) > 5 else "NONE"
CTAT_MAP    = sys.argv[6] if len(sys.argv) > 6 else "NONE"
CHR_STYLE   = sys.argv[7] if len(sys.argv) > 7 else "asis"   # "chr" | "plain" | "asis": reference naming

def ref_chrom(c):
    # FusionCatcher reports Ensembl names (12, MT); write them as the reference names them.
    c = str(c or "").strip()
    if CHR_STYLE == "chr" and c and not c.lower().startswith("chr"):
        return "chrM" if c in ("MT", "M") else "chr" + c
    if CHR_STYLE == "plain" and c.lower().startswith("chr"):
        c = c[3:]; return "MT" if c in ("M", "MT") else c
    return c

def norm_chrom(c):
    if c is None: return ""
    c = str(c).strip()
    c = re.sub(r'^chr', '', c, flags=re.I)
    if c in ('M', 'MT'): return 'MT'
    return c

def norm_gene(g):
    if g is None: return ""
    g = str(g).strip().upper()
    # drop STAR-Fusion ENSG suffix "SYMBOL^ENSG000..." and parenthetical notes
    g = g.split('^')[0]
    g = re.sub(r'\(.*?\)', '', g)
    return g.strip()

def clean_id(x):
    # "ENSG00000157764.14" / "ENSG00000182378.15_PAR_Y" -> "ENSG00000157764" / "ENSG00000182378"
    x = str(x or "").strip()
    return re.sub(r'(\.\d+)?(_PAR_Y)?$', '', x) if re.match(r'^ENS[A-Z]*G\d+', x) else ""

def load_map(path):
    id2sym = {}
    if path == "NONE": return id2sym
    try:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and clean_id(f[0]) and f[1]:
                    id2sym.setdefault(clean_id(f[0]), norm_gene(f[1]))
    except FileNotFoundError:
        pass
    return id2sym

ID2SYM_MAIN = load_map(MAIN_MAP)
ID2SYM_CTAT = load_map(CTAT_MAP)
# symbol -> gene ID, only where the symbol names exactly one gene (main GTF first, then CTAT)
SYM2ID = {}
for m in (ID2SYM_MAIN, ID2SYM_CTAT):
    seen = {}
    for gid, sym in m.items(): seen.setdefault(sym, set()).add(gid)
    for sym, ids in seen.items():
        if len(ids) == 1 and sym not in SYM2ID: SYM2ID[sym] = next(iter(ids))

def partner(sym, gid=""):
    # One fusion partner: (match token, reported symbol, CTAT symbol, all names seen)
    sym = norm_gene(sym)
    ids = [clean_id(x) for x in re.split(r'[,;]', str(gid or "")) if clean_id(x)]
    i = ids[0] if len(ids) == 1 else SYM2ID.get(sym, "")
    if not i:
        return ("SYM:" + sym, sym, sym, {sym})
    main = ID2SYM_MAIN.get(i) or ID2SYM_CTAT.get(i) or sym
    ctat = ID2SYM_CTAT.get(i) or main
    return (i, main, ctat, {sym, main, ctat})

def canon_key(p5, p3):
    return "--".join(sorted([p5[0], p3[0]]))   # order-insensitive for matching

def split_pair(line):
    # Separators: TAB, '--', comma, semicolon, space. A single hyphen is a separator
    # ONLY when the line has exactly one ('BCR-ABL1'); symbols such as NKX2-1 / HLA-A
    # contain hyphens, so 'A--B' or TAB is required for those.
    parts = [p for p in re.split(r'[\t,; ]+|--', line) if p]
    if len(parts) == 1 and parts[0].count('-') == 1:
        parts = parts[0].split('-')
    return parts

def parse_bp(s):
    # "chr7:55019278:+" or "7:55019278" -> (chrom, pos)
    if not s: return ("", "")
    parts = str(s).split(':')
    if len(parts) >= 2:
        return (str(parts[0]).strip(), parts[1])     # chromosome as the caller wrote it (reference naming)
    return ("", "")

def bp_strand(s):
    # gene strand of a "chr:pos:+" breakpoint (STAR-Fusion, FusionCatcher); "" when absent
    parts = str(s or "").split(':')
    return parts[2] if len(parts) >= 3 and parts[2] in ("+", "-") else ""

CONF_RANK = {"high": 3, "medium": 2, "low": 1}

def add(store, p5, p3, bp5, bp3, caller, support, extra, strands=("", "")):
    key, g5, g3 = canon_key(p5, p3), p5[1], p3[1]
    d = store.setdefault(key, {
        "genes": (g5, g3), "ctat": (p5[2], p3[2]), "bp5": bp5, "bp3": bp3,
        "callers": set(), "support": 0, "arriba_conf": "", "type": "", "known": False,
        "alias5": set(), "alias3": set(), "strands": ("", ""), "arriba_best": None
    })
    d["callers"].add(caller)
    d["alias5"] |= p5[3]; d["alias3"] |= p3[3]
    try: sup = int(float(support))
    except Exception: sup = 0
    d["support"] += sup
    # Breakpoints and directional gene order: Arriba's best row for the pair (highest confidence,
    # then most reads; Arriba lists a pair once per breakpoint, in no confidence order), else the
    # first caller that gave breakpoints. The pair's Arriba confidence is that of its best row.
    if caller == "arriba":
        rank = (CONF_RANK.get((extra.get("arriba_conf") or "").lower(), 0), sup)
        if d["arriba_best"] is None or rank > d["arriba_best"]:
            d["arriba_best"] = rank
            d["genes"] = (g5, g3); d["ctat"] = (p5[2], p3[2])
            d["bp5"] = bp5 or d["bp5"]; d["bp3"] = bp3 or d["bp3"]
            if extra.get("arriba_conf"): d["arriba_conf"] = extra["arriba_conf"]
            if extra.get("type"):        d["type"] = extra["type"]
            if all(strands):             d["strands"] = strands
    elif d["arriba_best"] is None and not d["bp5"] and bp5:
        d["genes"] = (g5, g3); d["ctat"] = (p5[2], p3[2])
        d["bp5"] = bp5; d["bp3"] = bp3 or d["bp3"]
    if all(strands) and not all(d["strands"]): d["strands"] = strands

store = {}

# ---- Arriba (tab; header starts with '#gene1') ----
try:
    with open("arriba.tsv") as fh:
        rows = [l.rstrip("\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].lstrip("#").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in rows[1:]:
            if line.startswith("#"): continue
            f = line.split("\t")
            def g(col, dflt=""):
                return f[idx[col]] if col in idx and idx[col] < len(f) else dflt
            g5, g3 = norm_gene(g("gene1")), norm_gene(g("gene2"))
            if not g5 or not g3: continue
            p5, p3 = partner(g5, g("gene_id1")), partner(g3, g("gene_id2"))
            bp5 = "{}:{}".format(*parse_bp(g('breakpoint1')))
            bp3 = "{}:{}".format(*parse_bp(g('breakpoint2')))
            sup = 0
            for c in ("split_reads1", "split_reads2", "discordant_mates"):
                try: sup += int(g(c, "0") or 0)
                except Exception: pass
            s5, s3 = g("strand1(gene/fusion)").split("/")[0], g("strand2(gene/fusion)").split("/")[0]
            add(store, p5, p3, bp5, bp3, "arriba", sup,
                {"arriba_conf": g("confidence"), "type": g("type")},
                (s5 if s5 in ("+", "-") else "", s3 if s3 in ("+", "-") else ""))
except FileNotFoundError:
    pass

# ---- STAR-Fusion (tab; header starts with '#FusionName') ----
try:
    with open("starfusion.tsv") as fh:
        rows = [l.rstrip("\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].lstrip("#").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in rows[1:]:
            if line.startswith("#"): continue
            f = line.split("\t")
            def g(col, dflt=""):
                return f[idx[col]] if col in idx and idx[col] < len(f) else dflt
            fname = g("FusionName")
            if "--" in fname:
                g5, g3 = fname.split("--")[0], fname.split("--")[1]
            else:
                g5, g3 = g("LeftGene"), g("RightGene")
            g5, g3 = norm_gene(g5), norm_gene(g3)
            if not g5 or not g3: continue
            # LeftGene / RightGene are "SYMBOL^ENSG..."
            id5 = g("LeftGene").split("^")[1] if "^" in g("LeftGene") else ""
            id3 = g("RightGene").split("^")[1] if "^" in g("RightGene") else ""
            p5, p3 = partner(g5, id5), partner(g3, id3)
            lb, rb = parse_bp(g("LeftBreakpoint")), parse_bp(g("RightBreakpoint"))
            bp5, bp3 = f"{lb[0]}:{lb[1]}", f"{rb[0]}:{rb[1]}"
            sup = 0
            for c in ("JunctionReadCount", "SpanningFragCount"):
                try: sup += int(float(g(c, "0") or 0))
                except Exception: pass
            add(store, p5, p3, bp5, bp3, "starfusion", sup, {},
                (bp_strand(g("LeftBreakpoint")), bp_strand(g("RightBreakpoint"))))
except FileNotFoundError:
    pass

# ---- FusionCatcher (tab; long header) ----
try:
    with open("fusioncatcher.txt") as fh:
        rows = [l.rstrip("\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].split("\t")
        def find(sub):
            for i, h in enumerate(hdr):
                if sub in h: return i
            return -1
        i_g1, i_g2 = find("Gene_1_symbol"), find("Gene_2_symbol")
        i_d1, i_d2 = find("Gene_1_id"), find("Gene_2_id")
        i_b1, i_b2 = find("Fusion_point_for_gene_1"), find("Fusion_point_for_gene_2")
        i_sp, i_su = find("Spanning_pairs"), find("Spanning_unique_reads")
        for line in rows[1:]:
            f = line.split("\t")
            if i_g1 < 0 or i_g2 < 0 or max(i_g1, i_g2) >= len(f): continue
            g5, g3 = norm_gene(f[i_g1]), norm_gene(f[i_g2])
            if not g5 or not g3: continue
            p5 = partner(g5, f[i_d1] if 0 <= i_d1 < len(f) else "")
            p3 = partner(g3, f[i_d2] if 0 <= i_d2 < len(f) else "")
            b1 = parse_bp(f[i_b1]) if 0 <= i_b1 < len(f) else ("", "")
            b2 = parse_bp(f[i_b2]) if 0 <= i_b2 < len(f) else ("", "")
            bp5, bp3 = f"{ref_chrom(b1[0])}:{b1[1]}", f"{ref_chrom(b2[0])}:{b2[1]}"
            sup = 0
            for i_c in (i_sp, i_su):
                if 0 <= i_c < len(f):
                    try: sup += int(float(f[i_c] or 0))
                    except Exception: pass
            add(store, p5, p3, bp5, bp3, "fusioncatcher", sup, {},
                (bp_strand(f[i_b1]) if 0 <= i_b1 < len(f) else "", bp_strand(f[i_b2]) if 0 <= i_b2 < len(f) else ""))
except FileNotFoundError:
    pass

# ---- known oncogenic fusions (order-insensitive pair set) ----
known = set()
if KNOWN_PATH != "NONE":
    try:
        with open(KNOWN_PATH) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = split_pair(line)
                if len(parts) >= 2:
                    known.add("--".join(sorted([norm_gene(parts[0]), norm_gene(parts[1])])))
                else:
                    sys.stderr.write(f"WARN known_fusions: unparsed line (use TAB or '--'): {line}\n")
    except FileNotFoundError:
        pass

def readthrough_like(bp5, bp3, strands):
    # Transcriptional read-through: both genes on the same strand, the 5' partner upstream of
    # the 3' partner on that strand, within readthrough_kb. Pairs in the other orientations
    # (tandem duplication, inversion: FGFR3::TACC3) are genomic rearrangements, not
    # read-through. Without gene strands the pair is not judged here (Arriba's own
    # read-through type still applies).
    try:
        c5, p5 = bp5.split(":")[:2]; c3, p3 = bp3.split(":")[:2]
        s5, s3 = strands
        if c5 == "" or norm_chrom(c5) != norm_chrom(c3) or not s5 or s5 != s3: return False
        d = int(p3) - int(p5) if s5 == "+" else int(p5) - int(p3)
        return 0 < d <= READTHRU_BP
    except Exception:
        return False

# ---- verdicts ----
# Callers that ran (a disabled caller's file is not staged). As in the DNA consensus, the
# threshold is capped at that number, so an Arriba-only run is not left without PASS calls.
import os
n_run = sum(os.path.exists(p) for p in ("arriba.tsv", "starfusion.tsv", "fusioncatcher.txt"))
EFF_MIN = max(1, min(MIN_CALLERS, n_run))
out_rows = []
n_pass = n_known = n_review = 0
for key, d in store.items():
    g5, g3 = d["genes"]
    ncall = len(d["callers"])
    is_known = any("--".join(sorted([a, b])) in known for a in d["alias5"] for b in d["alias3"])
    arriba_type = (d.get("type") or "").lower()
    readthrough = ("read_through" in arriba_type or "read-through" in arriba_type
                   or readthrough_like(d["bp5"], d["bp3"], d["strands"]))
    # A single-caller run has no vote; Arriba's own low-confidence calls are then left for review.
    single_low = (n_run == 1 and d["callers"] == {"arriba"} and (d["arriba_conf"] or "").lower() not in ("high", "medium"))
    if ncall >= EFF_MIN and not single_low:
        verdict = "PASS_CONSENSUS"; n_pass += 1
    elif is_known and ncall >= 1:
        verdict = "PASS_KNOWN"; n_known += 1
    else:
        verdict = "REVIEW"; n_review += 1
    out_rows.append([
        f"{g5}--{g3}", g5, g3, d["bp5"], d["bp3"],
        ",".join(sorted(d["callers"])), str(ncall), str(d["support"]),
        d["arriba_conf"] or ".", "yes" if is_known else "no",
        "yes" if readthrough else "no", verdict
    ])

# sort: PASS_CONSENSUS, PASS_KNOWN, REVIEW; then by caller count / support desc
order = {"PASS_CONSENSUS": 0, "PASS_KNOWN": 1, "REVIEW": 2}
out_rows.sort(key=lambda r: (order.get(r[11], 3), -int(r[6]), -int(r[7])))

with open(f"{SAMPLE}_fusion_consensus.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["FusionName","Gene5","Gene3","Breakpoint5","Breakpoint3",
                "Callers","NumCallers","TotalSupport","Arriba_Confidence",
                "KnownFusion","ReadThrough","Verdict"])
    w.writerows(out_rows)

# FusionInspector candidate list = confident set only (PASS_*), directional "A--B", in the CTAT
# library's symbols (FusionInspector looks genes up in its own GTF). The name table maps those
# names back to the reported names for the annotation step.
ctat_name = {f"{d['genes'][0]}--{d['genes'][1]}": f"{d['ctat'][0]}--{d['ctat'][1]}" for d in store.values()}
with open(f"{SAMPLE}_fusion_candidates.txt", "w") as out, open(f"{SAMPLE}_fusion_names.tsv", "w") as nm:
    nm.write("ctat_name\treported_name\n")
    seen = set()
    for r in out_rows:
        c = ctat_name.get(r[0], r[0])
        if r[11].startswith("PASS") and c not in seen:
            out.write(c + "\n"); nm.write(f"{c}\t{r[0]}\n"); seen.add(c)

with open(f"{SAMPLE}_fusion_consensus_stats.txt", "w") as out:
    out.write(f"sample_id\t{SAMPLE}\n")
    out.write(f"total_unique_fusions\t{len(store)}\n")
    out.write(f"pass_consensus\t{n_pass}\n")
    out.write(f"pass_known_single_caller\t{n_known}\n")
    out.write(f"review_only\t{n_review}\n")
    out.write(f"min_callers\t{MIN_CALLERS}\n")
    out.write(f"callers_run\t{n_run}\n")
    out.write(f"effective_min_callers\t{EFF_MIN}\n")
