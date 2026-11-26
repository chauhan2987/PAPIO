#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: GATK gCNV genotyped segments -> CNV calls for AnnotSV. GATK writes every segment
#          (copy-neutral ones too, GT 0) with ALT <DEL>,<DUP>; only the calls are kept
#          (GT 1 = DEL, GT 2 = DUP, relative to the contig's ploidy), each with its own ALT
#          and INFO/SVTYPE. Called by the PAPIO process ANNOTATE_SV (gcnv source).
# Usage  : gcnv_calls_to_sv.sh <segments.vcf.gz> <out.vcf>
# Output : <out.vcf> (header + non-reference segments)
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
IN=$1; OUT=$2

##### Main Script
zcat -f "$IN" | awk -F'\t' 'BEGIN { OFS = "\t" }
    /^##INFO=<ID=SVTYPE,/ { has = 1 }
    /^##/ { print; next }
    /^#CHROM/ { if (!has) print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"gCNV call from GT: DEL or DUP\">"; print; next }
    { n = split($9, f, ":"); gi = 0; for (i = 1; i <= n; i++) if (f[i] == "GT") gi = i
      split($10, v, ":"); gt = gi ? v[gi] : "."
      if (gt == "1") t = "DEL"; else if (gt == "2") t = "DUP"; else next
      $5 = "<" t ">"; $8 = ($8 == "." ? "" : $8 ";") "SVTYPE=" t; print }' > "$OUT"
