#!/usr/bin/awk -f
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Sex-aware copy-number calls for the patient report.
#          Called by the PAPIO process CONSOLIDATE_PATIENT_REPORT; moved from an inline heredoc.
# Usage  : awk -v src=sequenza|cnvkit -v sex=.. -v sexlab=.. -v ampmin=.. -f cn_calls.awk <segments>
#------------------------------------------------------------------------------
BEGIN { FS = OFS = "\t" }
function isauto(c) { return c ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ }
function isx(c)    { return c ~ /^(chr)?X$/ }
function isy(c)    { return c ~ /^(chr)?Y$/ }
# pass 1 (sequenza only): tumour ploidy = length-weighted mean autosomal CNt
NR == FNR && src == "sequenza" {
    if (FNR == 1) { for (i = 1; i <= NF; i++) h[$i] = i; next }
    cn = $(h["CNt"]); if (cn == "" || cn == "NA") next
    if (isauto($(h["chromosome"]))) { L = $(h["end.pos"]) - $(h["start.pos"]) + 1; sw += L * cn; sl += L }
    next
}
FNR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    if (src == "sequenza") { ploidy = (sl > 0) ? sw / sl : 2; base = int(ploidy + 0.5)
                             cc = "chromosome"; sc = "start.pos"; ec = "end.pos"; nc = "CNt"
                             note = (sl > 0) ? sprintf("tumour_ploidy=%.2f (length-weighted mean autosomal CNt)", ploidy) \
                                             : "tumour_ploidy=NA (no autosomal CNt; baseline 2 assumed)" }
    else                   { base = 2; cc = "chromosome"; sc = "start"; ec = "end"; nc = "cn"
                             note = "NOT purity/ploidy-adjusted" }
    printf "## cn_source=%s  %s  autosomal_baseline=%d  sex=%s  amp_min_cn=%s\n", src, note, base, sexlab, (ampmin == "" ? "unset" : ampmin)
    print "CHROM", "START", "END", "GENE", "LOG2", "COPY_NUMBER", "CALL"
    next
}
{
    c = $(h[cc]); cn = $(h[nc]); if (cn == "" || cn == "NA") next; cn += 0
    if (isauto(c))    ex = base
    else if (isx(c)) { if (sex == "1") ex = base / 2; else if (sex == "0") ex = base; else next }
    else if (isy(c)) { if (sex == "1") ex = base / 2; else next }
    else next
    if (cn == 0)                              call = "HOMDEL"
    else if (ampmin != "" && cn >= ampmin + 0) call = "AMP"
    else if (cn > ex)                         call = "GAIN"
    else if (cn < ex)                         call = "LOSS"
    else next
    if (src == "sequenza") { dr = $(h["depth.ratio"]); lg = (dr + 0 > 0) ? sprintf("%.3f", log(dr) / log(2)) : "."; g = "." }
    else                   { lg = $(h["log2"]); g = $(h["gene"]) }
    print c, $(h[sc]), $(h[ec]), g, lg, cn, call
}
