#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: CNVkit integer copy number. With a Sequenza fit of near-diploid ploidy
#          (1.6-2.4) the tumor purity is taken from it (cnvkit.py call -m clonal --purity);
#          otherwise CNVkit's default threshold calling is used: CNVkit centres log2 on the
#          median, so purity rescaling of a polyploid tumor would give wrong integer copy
#          numbers. Called by the PAPIO process CNVKIT_CALL.
# Usage  : cnvkit_call.sh <sample_id> <sample.cns> <sex: female|male|unknown> [<sequenza_confints_CP.txt>]
# Output : <sample>-call.cns
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=$1; CNS=$2; SEX=$3; CONFINTS=${4:-}

##### Main Script
sex_flag=()
case "$SEX" in female|male) sex_flag=(-x "$SEX") ;; esac

purity=""; ploidy=""
if [ -n "$CONFINTS" ] && [ -s "$CONFINTS" ]; then
    ## confints_CP.txt: header + (lower, point, upper) rows; the point estimate is data row 2
    read -r purity ploidy < <(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "cellularity") c = i; if ($i == "ploidy.estimate") p = i }; next }
                         NR == 3 && c { print $c, (p ? $p : "NA") }' "$CONFINTS") || true
fi
if [ -n "$purity" ] && ! awk -v p="$ploidy" 'BEGIN { exit !(p != "NA" && p + 0 >= 1.6 && p + 0 <= 2.4) }'; then
    echo "CNVKIT_CALL: Sequenza ploidy ${ploidy:-NA} outside 1.6-2.4 - purity not applied" >&2
    purity=""
fi

if [ -n "$purity" ]; then
    echo "CNVKIT_CALL: purity $purity from Sequenza (ploidy $ploidy; -m clonal)" >&2
    cnvkit.py call "$CNS" -m clonal --purity "$purity" "${sex_flag[@]}" -o "${SAMPLE}-call.cns"
else
    echo "CNVKIT_CALL: CNVkit default threshold calling (no purity applied)" >&2
    cnvkit.py call "$CNS" "${sex_flag[@]}" -o "${SAMPLE}-call.cns"
fi
