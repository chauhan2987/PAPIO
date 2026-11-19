#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Annotate and tier consensus fusions.
#          Called by the PAPIO process ANNOTATE_FUSIONS; moved from an inline heredoc.
# Usage  : annotate_fusions.py <cgc_arg> <known_arg> <sample_id> <min_known>
#------------------------------------------------------------------------------
import sys, re

CGC_PATH   = sys.argv[1]
KNOWN_PATH = sys.argv[2]
SAMPLE     = sys.argv[3]

def norm_gene(g):
    g = str(g or "").strip().upper().split('^')[0]
    return re.sub(r'\(.*?\)', '', g).strip()

def split_pair(line):
    # Separators: TAB, '--', comma, semicolon, space. A single hyphen is a separator
    # ONLY when the line has exactly one ('BCR-ABL1'); symbols such as NKX2-1 / HLA-A
    # contain hyphens, so 'A--B' or TAB is required for those.
    parts = [p for p in re.split(r'[\t,; ]+|--', line) if p]
    if len(parts) == 1 and parts[0].count('-') == 1:
        parts = parts[0].split('-')
    return parts

# CGC gene -> role (col layout matches the DNA CLASSIFY_VARIANTS usage:
# col1=symbol, col15=role_in_cancer). Robust to header presence.
role = {}
if CGC_PATH != "NONE":
    # Columns located by HEADER NAME (CGC layouts differ between releases).
    with open(CGC_PATH) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    if lines:
        hdr_cgc = [re.sub(r"\s+", "_", x.strip().upper()) for x in lines[0]]
        if "GENE_SYMBOL" not in hdr_cgc or "ROLE_IN_CANCER" not in hdr_cgc:
            sys.exit("ERROR: CGC header lacks GENE_SYMBOL / ROLE_IN_CANCER: " + " | ".join(lines[0]))
        gi, ri = hdr_cgc.index("GENE_SYMBOL"), hdr_cgc.index("ROLE_IN_CANCER")
        for f in lines[1:]:
            if len(f) > max(gi, ri) and f[gi]:
                role[norm_gene(f[gi])] = f[ri]

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
    except FileNotFoundError:
        pass

# FusionInspector confirmation set (directional A--B, upper). fi.tsv is staged only
# when FusionInspector ran for this sample.
import os
FI_RAN    = os.path.exists("fi.tsv")
MIN_KNOWN = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
# FusionInspector reports the CTAT library's gene symbols; names.tsv (from FUSION_CONSENSUS)
# maps them to the reported names, which can differ when the RNA GTF is a newer GENCODE release.
to_reported = {}
try:
    with open("names.tsv") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] != "ctat_name":
                to_reported["--".join(norm_gene(x) for x in f[0].split("--"))] = \
                    "--".join(norm_gene(x) for x in f[1].split("--"))
except FileNotFoundError:
    pass

confirmed = {}
try:
    with open("fi.tsv") as fh:
        rows = [l.rstrip("\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].lstrip("#").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in rows[1:]:
            f = line.split("\t")
            name = f[idx["FusionName"]] if "FusionName" in idx and idx["FusionName"] < len(f) else ""
            name = "--".join(norm_gene(x) for x in name.split("--")) if "--" in name else norm_gene(name)
            name = to_reported.get(name, name)
            jrc = f[idx["JunctionReadCount"]] if "JunctionReadCount" in idx and idx["JunctionReadCount"] < len(f) else "0"
            confirmed[name] = jrc
except FileNotFoundError:
    pass

# Immunoglobulin / T-cell receptor segment genes and HLA. Anchored so that ordinary genes
# sharing a prefix (TRAF3, TRDN, IGLON5, IGHMBP2, ...) are not caught.
IG_HLA = re.compile(r'^(IG[HKL]@?|IG[HKL][VDJ]\d.*|IGH(A[12]|D|E|G[1-4]|M)|IGKC|IGLC\d+|'
                    r'TR[ABDG]@?|TR[ABDG][VDJ]\d.*|TR[ABDG]C\d?|HLA-.*)$')

out = []
with open("consensus.tsv") as fh:
    rows = [l.rstrip("\n") for l in fh if l.strip()]
hdr = rows[0].split("\t") if rows else []
cidx = {h: i for i, h in enumerate(hdr)}
def col(f, name, d=""):
    return f[cidx[name]] if name in cidx and cidx[name] < len(f) else d

for line in rows[1:]:
    f = line.split("\t")
    g5, g3 = norm_gene(col(f, "Gene5")), norm_gene(col(f, "Gene3"))
    key = "--".join(sorted([g5, g3]))
    dirname = f"{g5}--{g3}"
    verdict = col(f, "Verdict")
    readthrough = col(f, "ReadThrough") == "yes"
    is_known = key in known or col(f, "KnownFusion") == "yes"
    r5, r3 = role.get(g5, ""), role.get(g3, "")
    in_cgc = bool(r5) or bool(r3)
    fi_reads = confirmed.get(dirname, "")          # directional: FI confirms this 5'->3' product only
    fi_confirmed = "yes" if (fi_reads not in ("", None) and str(fi_reads).strip() not in ("", "0")) else "no"
    ig_hla = bool(IG_HLA.match(g5) or IG_HLA.match(g3))
    # IG/TR partner + a cancer gene (IGH--MYC, IGH--BCL2, IGH--CCND1) is a driver rearrangement,
    # not an artifact. HLA, IG/TR--IG/TR, and IG/TR with a non-cancer gene stay LIKELY_ARTIFACT
    # (abundant immune transcripts in infiltrated tissue).
    ig5, ig3 = bool(IG_HLA.match(g5)), bool(IG_HLA.match(g3))
    partner_cgc = (ig5 != ig3) and bool(role.get(g3 if ig5 else g5, ""))
    hla = g5.startswith("HLA-") or g3.startswith("HLA-")
    ig_artifact = ig_hla and (hla or not partner_cgc)

    # tiering (computational; for pathologist review - not a final AMP/CAP tier)
    # A catalogued pair needs real support: multi-caller consensus, OR a single-caller
    # known rescue with >= MIN_KNOWN reads AND FusionInspector confirmation when FI ran
    # (1-read single-caller hits are typical of index hopping / cross-sample bleed).
    try:
        support = float(col(f, "TotalSupport", "0") or 0)
    except ValueError:
        support = 0.0
    known_supported = (verdict == "PASS_CONSENSUS") or (
        verdict.startswith("PASS") and support >= MIN_KNOWN and (not FI_RAN or fi_confirmed == "yes"))
    if is_known and known_supported:
        tier = "KNOWN_ONCOGENIC"
    elif is_known:
        tier = "KNOWN_LOW_SUPPORT"
    elif readthrough or ig_artifact:
        tier = "LIKELY_ARTIFACT"
    elif in_cgc and verdict == "PASS_CONSENSUS":
        tier = "LIKELY_ONCOGENIC"
    elif verdict == "PASS_CONSENSUS":
        tier = "VUS"
    else:
        tier = "LOW_CONFIDENCE"

    out.append([
        dirname, g5, g3, col(f, "Breakpoint5"), col(f, "Breakpoint3"),
        col(f, "Callers"), col(f, "NumCallers"), col(f, "TotalSupport"),
        col(f, "Arriba_Confidence"), "yes" if is_known else "no",
        "yes" if in_cgc else "no", (r5 or ".") + "|" + (r3 or "."),
        "yes" if readthrough else "no", "yes" if ig_hla else "no",
        fi_confirmed, fi_reads or ".", verdict, tier
    ])

order = {"KNOWN_ONCOGENIC":0,"LIKELY_ONCOGENIC":1,"KNOWN_LOW_SUPPORT":2,
         "VUS":3,"LOW_CONFIDENCE":4,"LIKELY_ARTIFACT":5}
out.sort(key=lambda r: (order.get(r[17], 9), r[0]))

with open(f"{SAMPLE}_fusions_annotated.tsv", "w") as o:
    o.write("# COMPUTATIONAL FUSION PRE-CLASSIFICATION - NOT A FINAL AMP/ASCO/CAP CLINICAL TIER\n")
    o.write("# Final tiering requires a curated knowledge base (OncoKB/CIViC/COSMIC) + expert review.\n")
    o.write("\t".join([
        "FusionName","Gene5","Gene3","Breakpoint5","Breakpoint3","Callers","NumCallers",
        "TotalSupport","Arriba_Confidence","KnownFusion","In_CGC","CGC_Roles(5|3)",
        "ReadThrough","IG_HLA_Locus","FI_Confirmed","FI_JunctionReads","Verdict","TIER"
    ]) + "\n")
    for r in out:
        o.write("\t".join(r) + "\n")
