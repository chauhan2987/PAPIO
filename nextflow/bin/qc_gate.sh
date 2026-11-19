#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Soft QC gate for one DNA sample (report only; never stops a run).
#          Called by the PAPIO process QC_GATE.
# Usage  : qc_gate.sh --sample ID --summary MOSDEPTH_SUMMARY [--selfsm FILE] [--hs FILE]
#                     [--stats FILE] --min-cov N --max-contam F --min-breadth30 F
#                     --min-breadth20 F --gate20 true|false --min-ontarget F
#                     --min-mapping F --max-fold80 F
#          [--global-dist MOSDEPTH_GLOBAL_DIST]  (WGS: breadth rows = share of the autosomes,
#          length-weighted, at >= 30x / 20x; used when there are no HsMetrics)
#          --gate20 true (normals, germline): breadth judged at 20x, the 30x row is INFO;
#          false (tumors): breadth judged at 30x, the 20x row is INFO.
# Output : <sample>_qc_gate.tsv
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; SUMMARY=""; SELFSM=""; HS=""; STATS=""; GDIST=""
MIN_COV=""; MAX_CONTAM=""; MIN_B30=""; MIN_B20=""; GATE20="false"
MIN_ONT=""; MIN_MAP=""; MAX_F80=""

##### Functions
## Function to read one column (by header name) from the first data row of a Picard metrics block
picard_value() {
    awk -F'\t' -v col="$1" '/^BAIT_SET/ { for (i = 1; i <= NF; i++) if ($i == col) c = i; getline; if (c) print $c }' "$2"
}

## Function to flag a value against a threshold (dir ">=" or "<="; missing -> the given word)
flag() {
    awk -v v="$1" -v t="$2" -v d="$3" -v na="$4" \
        'BEGIN { if (v == "NA") print na; else if (d == ">=") print ((v + 0 >= t + 0) ? "PASS" : "WARN"); else print ((v + 0 <= t + 0) ? "PASS" : "WARN") }'
}

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "qc_gate.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)        SAMPLE=$2 ;;
        --summary)       SUMMARY=$2 ;;
        --selfsm)        SELFSM=$2 ;;
        --hs)            HS=$2 ;;
        --global-dist)   GDIST=$2 ;;
        --stats)         STATS=$2 ;;
        --min-cov)       MIN_COV=$2 ;;
        --max-contam)    MAX_CONTAM=$2 ;;
        --min-breadth30) MIN_B30=$2 ;;
        --min-breadth20) MIN_B20=$2 ;;
        --gate20)        GATE20=$2 ;;
        --min-ontarget)  MIN_ONT=$2 ;;
        --min-mapping)   MIN_MAP=$2 ;;
        --max-fold80)    MAX_F80=$2 ;;
        *) echo "qc_gate.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

## mean coverage: on-target (total_region) when mosdepth ran --by, else genome-wide (total)
cov=$(awk -F'\t' '$1 == "total_region" { print $4 }' "$SUMMARY")
[ -n "$cov" ] || cov=$(awk -F'\t' '$1 == "total" { print $4 }' "$SUMMARY")
cov=${cov:-NA}

## contamination (VerifyBamID2 FREEMIX; normals only)
contam="NA"
if [ -n "$SELFSM" ] && [ -s "$SELFSM" ]; then
    contam=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "FREEMIX") c = i } NR == 2 && c { print $c }' "$SELFSM")
    contam=${contam:-NA}
fi

## Picard HsMetrics (WES only): breadth, on-target fraction, fold-80 penalty
breadth="NA"; breadth20="NA"; ontarget="NA"; fold80="NA"
if [ -n "$HS" ] && [ -s "$HS" ]; then
    breadth=$(picard_value PCT_TARGET_BASES_30X "$HS");  breadth=${breadth:-NA}
    breadth20=$(picard_value PCT_TARGET_BASES_20X "$HS"); breadth20=${breadth20:-NA}
    ontarget=$(picard_value PCT_SELECTED_BASES "$HS");   ontarget=${ontarget:-NA}
    fold80=$(picard_value FOLD_80_BASE_PENALTY "$HS");   fold80=${fold80:-NA}
fi

## WGS: genome breadth from the mosdepth global distribution (autosomes, weighted by length)
if [ "$breadth" = "NA" ] && [ -n "$GDIST" ] && [ -s "$GDIST" ]; then
    gbreadth() {
        awk -F'\t' -v d="$1" 'NR == FNR { if ($1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/) len[$1] = $2; next }
            ($1 in len) && $2 + 0 == d + 0 { s += len[$1] * $3; l += len[$1] }
            END { if (l > 0) printf "%.4f", s / l; else print "NA" }' "$SUMMARY" "$GDIST"
    }
    breadth=$(gbreadth 30); breadth20=$(gbreadth 20)
fi

## mapping rate = reads mapped / raw total sequences (samtools stats SN section)
mapping="NA"
if [ -n "$STATS" ] && [ -s "$STATS" ]; then
    mapped=$(awk -F'\t' '/^SN\treads mapped:/ { print $3; exit }' "$STATS")
    total=$(awk -F'\t' '/^SN\traw total sequences:/ { print $3; exit }' "$STATS")
    if [ -n "$mapped" ] && [ -n "$total" ] && [ "$total" != "0" ]; then
        mapping=$(awk -v m="$mapped" -v t="$total" 'BEGIN { printf "%.4f", m / t }')
    fi
fi

cov_flag=$(flag "$cov" "$MIN_COV" ">=" UNKNOWN)
con_flag=$(flag "$contam" "$MAX_CONTAM" "<=" UNKNOWN)
brd_flag=$(flag "$breadth" "$MIN_B30" ">=" SKIP)
b20_flag=$(flag "$breadth20" "$MIN_B20" ">=" SKIP)
if [ "$GATE20" = "true" ]; then [ "$brd_flag" = "SKIP" ] || brd_flag="INFO"
else [ "$b20_flag" = "SKIP" ] || b20_flag="INFO"; fi
ont_flag=$(flag "$ontarget" "$MIN_ONT" ">=" SKIP)
map_flag=$(flag "$mapping" "$MIN_MAP" ">=" UNKNOWN)
f80_flag=$(flag "$fold80" "$MAX_F80" "<=" SKIP)

{   printf 'sample_id\t%s\n' "$SAMPLE"
    printf 'metric\tvalue\tthreshold\tflag\n'
    printf 'mean_coverage\t%s\t>=%s\t%s\n'         "$cov"       "$MIN_COV"    "$cov_flag"
    printf 'contamination_freemix\t%s\t<=%s\t%s\n' "$contam"    "$MAX_CONTAM" "$con_flag"
    printf 'pct_target_30x\t%s\t>=%s\t%s\n'        "$breadth"   "$MIN_B30"    "$brd_flag"
    printf 'pct_target_20x\t%s\t>=%s\t%s\n'        "$breadth20" "$MIN_B20"    "$b20_flag"
    printf 'pct_on_target\t%s\t>=%s\t%s\n'         "$ontarget"  "$MIN_ONT"    "$ont_flag"
    printf 'mapping_rate\t%s\t>=%s\t%s\n'          "$mapping"   "$MIN_MAP"    "$map_flag"
    printf 'fold_80_penalty\t%s\t<=%s\t%s\n'       "$fold80"    "$MAX_F80"    "$f80_flag"
    printf 'note\tT/N_concordance_assessed_in_somalier.pairs.tsv\t-\t-\n'
} > "${SAMPLE}_qc_gate.tsv"
