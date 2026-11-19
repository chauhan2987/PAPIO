#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Infer sex from the chrX and chrY / autosome coverage ratios and compare it with
#          the samplesheet. Called by the PAPIO process INFER_SEX.
# Usage  : infer_sex.sh --sample ID --summary MOSDEPTH_SUMMARY --expected 0|1|NA
#                       --type DT|DN --female-min F --male-max F [--female-max-y F --male-min-y F]
#          chrX decides; chrY, when it has data, must agree: female needs chrY <= female-max-y,
#          male needs chrY >= male-min-y (normals; tumors may have lost chrY). Disagreement
#          (X0, XXY, contamination) is AMBIGUOUS.
#          WES (mosdepth --by targets): on-target *_region rows are used, because
#          whole-chromosome means are dominated by off-target reads.
# Output : <sample>_sex.tsv
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; SUMMARY=""; EXPECTED_CODE="NA"; TYPE=""; FEMALE_MIN=""; MALE_MAX=""; FEMALE_MAX_Y="0.10"; MALE_MIN_Y="0.20"

##### Functions
## Function to print the mean depth of the autosomes (chr1-22 or 1-22)
autosome_mean() {
    awk -v s="$1" '{ c = $1; if (s != "") { if (c !~ /_region$/) next; sub(/_region$/, "", c) } else if (c ~ /_region$/) next }
        c ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ { sum += $4; n++ }
        END { if (n > 0) printf "%.4f", sum / n; else print "NA" }' "$2"
}

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "infer_sex.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)     SAMPLE=$2 ;;
        --female-max-y) FEMALE_MAX_Y=$2 ;;
        --male-min-y)   MALE_MIN_Y=$2 ;;
        --summary)    SUMMARY=$2 ;;
        --expected)   EXPECTED_CODE=$2 ;;
        --type)       TYPE=$2 ;;
        --female-min) FEMALE_MIN=$2 ;;
        --male-max)   MALE_MAX=$2 ;;
        *) echo "infer_sex.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

if grep -q '_region' "$SUMMARY"; then suf="_region"; else suf=""; fi
auto=$(autosome_mean "$suf" "$SUMMARY")
xdep=$(awk -v s="$suf" '$1 == ("chrX" s) || $1 == ("X" s) { print $4; exit }' "$SUMMARY")
xdep=${xdep:-NA}
ydep=$(awk -v s="$suf" '$1 == ("chrY" s) || $1 == ("Y" s) { print $4; exit }' "$SUMMARY")
ydep=${ydep:-NA}

ratio="NA"; yratio="NA"; inferred="UNKNOWN"
if [ "$auto" != "NA" ] && [ "$xdep" != "NA" ] && awk -v a="$auto" 'BEGIN { exit !(a > 0) }'; then
    ratio=$(awk -v x="$xdep" -v a="$auto" 'BEGIN { printf "%.3f", x / a }')
    [ "$ydep" = "NA" ] || yratio=$(awk -v y="$ydep" -v a="$auto" 'BEGIN { printf "%.3f", y / a }')
    inferred=$(awk -v r="$ratio" -v y="$yratio" -v fmin="$FEMALE_MIN" -v mmax="$MALE_MAX" \
                   -v fmaxy="$FEMALE_MAX_Y" -v mminy="$MALE_MIN_Y" -v t="$TYPE" 'BEGIN {
        x = (r + 0 >= fmin + 0) ? "female" : (r + 0 <= mmax + 0) ? "male" : "AMBIGUOUS"
        if (y != "NA" && x == "female" && y + 0 > fmaxy + 0) x = "AMBIGUOUS"
        if (y != "NA" && x == "male" && y + 0 < mminy + 0 && t != "DT") x = "AMBIGUOUS"
        print x }')
fi

case "$EXPECTED_CODE" in
    0) expected="female" ;;
    1) expected="male" ;;
    *) expected="unknown" ;;
esac

if [ "$inferred" = "UNKNOWN" ] || [ "$expected" = "unknown" ]; then
    flag="SKIP"
elif [ "$inferred" = "AMBIGUOUS" ]; then
    flag="REVIEW"           # between the cut-offs, or chrX and chrY disagree (X0, XXY, contamination)
elif [ "$inferred" = "$expected" ]; then
    flag="PASS"
elif [ "$TYPE" = "DT" ]; then
    flag="REVIEW_TUMOR"     # chrX loss / loss of Y in tumors can explain it - judge on the normal
else
    flag="MISMATCH"
fi

{   printf 'sample_id\t%s\n' "$SAMPLE"
    printf 'chrX_autosome_ratio\t%s\n' "$ratio"
    printf 'chrY_autosome_ratio\t%s\n' "$yratio"
    printf 'inferred_sex\t%s\n' "$inferred"
    printf 'expected_sex\t%s\n' "$expected"
    printf 'sex_concordance\t%s\n' "$flag"
} > "${SAMPLE}_sex.tsv"
