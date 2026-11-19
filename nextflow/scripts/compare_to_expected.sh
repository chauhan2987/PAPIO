#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Date   : 2026-11-12
# Purpose: Record / check the expected outputs of a PAPIO test run.
#
#   record <Result dir> <expected dir>   developer: writes expected_metrics.tsv
#   check  <Result dir> <expected dir>   user: compares a run against it
#
#   Metrics (per key output file, keyed by path relative to Result/):
#     - file present
#     - VCF: records and PASS records
#     - TSV/TXT: data rows (lines not starting with '#')
#     - fusion consensus: the exact set of PASS fusions (name|verdict)
#   Counts pass if within TOL_PCT percent or TOL_ABS records of expected
#   (multi-threaded callers are not bit-identical across machines);
#   fusion sets must match exactly.
#
# Exit   : 0 = all metrics pass, 1 = any FAIL
#------------------------------------------------------------------------------
set -uo pipefail

###### GLOBAL VARIABLES
MODE=${1:?usage: compare_to_expected.sh record|check <Result dir> <expected dir>}
RESULT=${2:?missing Result dir}
EXPECTED=${3:?missing expected dir}
TOL_PCT=${TOL_PCT:-5}
TOL_ABS=${TOL_ABS:-2}
## key output files (find -name patterns)
PATTERNS="*_consensus.vcf.gz *_filtered.vcf.gz *_vep.vcf.gz *_consensus_stats.txt *_preclassification.tsv
*consolidated*.tsv *_fusion_consensus.tsv *_fusions_annotated.tsv *_fusion_report.tsv *_msi"

#------------------------------------------------------------------------------
##### Functions

## Function to compute "key<TAB>metric<TAB>value" lines for a Result dir
metrics() {
    local root=$1 pat f rel
    for pat in $PATTERNS; do
        find "$root" -path '*/work/*' -prune -o -type f -name "$pat" -print
    done | LC_ALL=C sort -u | while read -r f; do
        rel=${f#"$root"/}
        printf '%s\tpresent\t1\n' "$rel"
        case "$f" in
            *.vcf.gz)
                zcat "$f" | awk -F'\t' -v k="$rel" '!/^#/ {n++; if ($7=="PASS") p++}
                    END {printf "%s\trecords\t%d\n%s\tpass_records\t%d\n", k, n, k, p}' ;;
            *_fusion_consensus.tsv)
                awk -F'\t' -v k="$rel" 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
                    {n++; if ($h["Verdict"] ~ /^PASS/) s[$h["FusionName"] "|" $h["Verdict"]]=1}
                    END {printf "%s\trows\t%d\n", k, n; m=""; for (x in s) m = m (m==""?"":",") x
                         printf "%s\tpass_fusions\t%s\n", k, (m==""?"none":m)}' "$f" ;;
            *)
                awk -v k="$rel" '!/^#/ && NF {n++} END {printf "%s\trows\t%d\n", k, n}' "$f" ;;
        esac
    done
}

## Function to sort the comma list inside pass_fusions so order never matters
## (portable: no gawk-only asort)
canon() {
    local k m v
    while IFS=$'\t' read -r k m v; do
        [ "$m" = "pass_fusions" ] && v=$(echo "$v" | tr ',' '\n' | LC_ALL=C sort | paste -sd, -)
        printf '%s\t%s\t%s\n' "$k" "$m" "$v"
    done
}

#------------------------------------------------------------------------------
##### Main Script

[ -d "$RESULT" ] || { echo "ERROR: Result dir not found: $RESULT" >&2; exit 1; }

if [ "$MODE" = "record" ]; then
    mkdir -p "$EXPECTED"
    metrics "$RESULT" | canon > "${EXPECTED}/expected_metrics.tsv"
    echo "Recorded $(cut -f1 "${EXPECTED}/expected_metrics.tsv" | sort -u | wc -l) files, $(wc -l < "${EXPECTED}/expected_metrics.tsv") metrics -> ${EXPECTED}/expected_metrics.tsv"
    exit 0
fi

[ "$MODE" = "check" ] || { echo "ERROR: mode must be record or check" >&2; exit 1; }
[ -s "${EXPECTED}/expected_metrics.tsv" ] || { echo "ERROR: ${EXPECTED}/expected_metrics.tsv missing" >&2; exit 1; }
metrics "$RESULT" | canon > "${RESULT}/.observed_metrics.tsv"

awk -F'\t' -v tp="$TOL_PCT" -v ta="$TOL_ABS" '
    NR == FNR { obs[$1 FS $2] = $3; next }
    {
        key = $1 FS $2; e = $3
        if (!(key in obs)) { printf "  FAIL  %-60s %-13s missing\n", $1, $2; bad++; next }
        o = obs[key]
        if ($2 == "pass_fusions") {
            if (o == e) { printf "  PASS  %-60s %-13s identical\n", $1, $2; ok++ }
            else { printf "  FAIL  %-60s %-13s expected [%s] got [%s]\n", $1, $2, e, o; bad++ }
            next
        }
        d = o - e; if (d < 0) d = -d
        lim = e * tp / 100; if (lim < ta) lim = ta
        if (d <= lim) { printf "  PASS  %-60s %-13s %s (expected %s)\n", $1, $2, o, e; ok++ }
        else          { printf "  FAIL  %-60s %-13s %s (expected %s, tolerance %d)\n", $1, $2, o, e, lim; bad++ }
    }
    END { printf "\n%d PASS, %d FAIL\n", ok, bad; exit (bad > 0) }' \
    "${RESULT}/.observed_metrics.tsv" "${EXPECTED}/expected_metrics.tsv"
