#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Cohort QC table from the per-sample QC gate reports (one row per sample).
#          Called by the PAPIO process QC_GATE_SUMMARY.
# Usage  : qc_gate_summary.sh <sample>_qc_gate.tsv [...]
# Output : cohort_qc_summary.tsv (OVERALL = WARN when any gated metric is WARN)
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
METRICS="mean_coverage contamination_freemix pct_target_30x pct_target_20x pct_on_target mapping_rate fold_80_penalty"
OUT="cohort_qc_summary.tsv"

##### Functions
## Function to write one sample's row
sample_row() {
    awk -F'\t' -v metrics="$METRICS" '
        $1 == "sample_id" { sid = $2; next }
        { val[$1] = $2; flg[$1] = $4 }
        END {
            n = split(metrics, m, " "); overall = "PASS"; line = sid
            for (i = 1; i <= n; i++) {
                v = (m[i] in val && val[m[i]] != "") ? val[m[i]] : "NA"
                if (flg[m[i]] == "WARN") overall = "WARN"
                line = line "\t" v
            }
            print line "\t" overall
        }' "$1"
}

##### Main Script
{   printf 'SAMPLE'; for m in $METRICS; do printf '\t%s' "$m"; done; printf '\tOVERALL\n'
    for f in "$@"; do
        [ -s "$f" ] || continue
        sample_row "$f"
    done | LC_ALL=C sort -t$'\t' -k1,1
} > "$OUT"
