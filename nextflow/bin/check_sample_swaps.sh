#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Sample-identity verdicts from somalier relate pairs.
#          Called by the PAPIO process CHECK_SAMPLE_SWAPS.
# Usage  : check_sample_swaps.sh --pairs somalier.pairs.tsv --analysis somatic|germline
#                                --min-rel F --max-unrel F --max-ibs0 F [--dup-min F]
#          Somatic: samples of one patient must be related and different patients must
#          not be. Low relatedness with an opposite-homozygote (IBS0) fraction <= max-ibs0
#          is tumor LOH in one person, not a swap. Germline: only identity-level matches
#          (relatedness > dup-min, default 0.80) between different sample IDs are flagged;
#          relatives are expected.
#          The patient is the sample ID without its last "_<label>".
# Output : sample_swap_report.tsv
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
PAIRS=""; ANALYSIS="germline"; MIN_REL=""; MAX_UNREL=""; MAX_IBS0=""; DUP_MIN="0.80"
OUT="sample_swap_report.tsv"

##### Functions
## Function to print the 1-based position of a column in the pairs header
col() { echo "$HDR" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1 || true; }

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "check_sample_swaps.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --dup-min)   DUP_MIN=$2 ;;
        --pairs)     PAIRS=$2 ;;
        --analysis)  ANALYSIS=$2 ;;
        --min-rel)   MIN_REL=$2 ;;
        --max-unrel) MAX_UNREL=$2 ;;
        --max-ibs0)  MAX_IBS0=$2 ;;
        *) echo "check_sample_swaps.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done
germ=$([ "$ANALYSIS" = "germline" ] && echo true || echo false)

HDR=$(head -1 "$PAIRS" | sed 's/^#//')
ca=$(col sample_a); cb=$(col sample_b); cr=$(col relatedness)
ci0=$(col ibs0); cn=$(col n)          # optional: absent -> relatedness-only rule
if [ -z "$ca" ] || [ -z "$cb" ] || [ -z "$cr" ]; then
    echo "ERROR: could not find sample_a/sample_b/relatedness columns in $PAIRS" >&2
    head -1 "$PAIRS" >&2; exit 1
fi

awk -F'\t' -v ca="$ca" -v cb="$cb" -v cr="$cr" -v ci0="${ci0:-0}" -v cn="${cn:-0}" \
    -v minrel="$MIN_REL" -v maxunrel="$MAX_UNREL" -v maxibs0="$MAX_IBS0" -v germ="$germ" -v dupmin="$DUP_MIN" '
    function patient(s,   n) { n = s; sub(/_[^_]+$/, "", n); return n }
    BEGIN { OFS = "\t"; print "sample_a", "sample_b", "relatedness", "same_patient", "verdict"
            swap = 0; contam = 0; dup = 0; review = 0 }
    NR > 1 && $ca != "" {
        a = $ca; b = $cb; rel = $cr + 0
        same = (patient(a) == patient(b)) ? "yes" : "no"
        verdict = "PASS"
        if (germ == "true") {
            if (rel > dupmin + 0) { verdict = "FLAG_POSSIBLE_DUPLICATE"; dup++; print a, b, $cr, same, verdict }
        } else {
            if (same == "yes" && rel < minrel) {
                if (ci0 > 0 && cn > 0 && $cn + 0 > 0 && ($ci0 / $cn) <= maxibs0) {
                    verdict = "REVIEW_LOW_RELATEDNESS_IBS0_OK"; review++
                } else { verdict = "FLAG_SWAP_LOW_RELATEDNESS"; swap++ }
            }
            if (same == "no" && rel > maxunrel) { verdict = "FLAG_UNEXPECTED_RELATEDNESS"; contam++ }
            if (same == "yes" || verdict != "PASS") print a, b, $cr, same, verdict
        }
    }
    END {
        if (germ == "true") {
            print "# possible_duplicates=" dup "  (relatedness>" dupmin " between distinct sample IDs)" > "/dev/stderr"
        } else {
            print "# REVIEW_LOW_RELATEDNESS_IBS0_OK = IBS0~0: tumour LOH in one person, OR a parent-child pair (also IBS0~0) - confirm identity"
            print "# same_patient_swaps=" swap "  loh_compatible_reviews=" review "  unexpected_relatedness_pairs=" contam \
                  "  (min_rel=" minrel ", max_unrel=" maxunrel ")" > "/dev/stderr"
        }
    }' "$PAIRS" > "$OUT"

## data rows only: exclude the column header and the '#' legend line
n=$(awk 'NR > 1 && !/^#/' "$OUT" | wc -l)
echo "CHECK_SAMPLE_SWAPS ($ANALYSIS): $n flagged/relevant pair(s)" >&2
