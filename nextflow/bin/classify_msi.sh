#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: MSI status from the MSIsensor-pro score. Called by the PAPIO process CLASSIFY_MSI.
# Usage  : classify_msi.sh --sample ID --msi <sample>_msi --high F --low F --min-sites N
#          score >= high -> MSI-H; score < low -> MSS; in between -> MSI-INDETERMINATE;
#          fewer than min-sites evaluable microsatellites -> NA_LOW_SITES.
# Output : <sample>_msi_status.txt
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; MSI=""; HIGH=""; LOW=""; MIN_SITES=50

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "classify_msi.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)    SAMPLE=$2 ;;
        --msi)       MSI=$2 ;;
        --high)      HIGH=$2 ;;
        --low)       LOW=$2 ;;
        --min-sites) MIN_SITES=$2 ;;
        *) echo "classify_msi.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

if [ ! -s "$MSI" ]; then
    echo "ERROR: MSIsensor-pro main output '$MSI' not found among staged files" >&2
    ls -l >&2
    exit 1
fi

score=$(awk 'NR == 2 && $3 != "" { print $3 }' "$MSI")
if [ -z "$score" ]; then
    score="NA"; status="NA"
else
    status=$(awk -v s="$score" -v hi="$HIGH" -v lo="$LOW" \
        'BEGIN { if (s + 0 >= hi + 0) print "MSI-H"; else if (s + 0 < lo + 0) print "MSS"; else print "MSI-INDETERMINATE" }')
fi
total=$(awk 'NR == 2 { print $1 }' "$MSI")
somatic=$(awk 'NR == 2 { print $2 }' "$MSI")
## too few evaluable microsatellites: the score is not interpretable
if [ -n "$total" ] && [ "$status" != "NA" ] && [ "${total%.*}" -lt "$MIN_SITES" ]; then
    status="NA_LOW_SITES"
fi

{   printf 'sample_id\t%s\n' "$SAMPLE"
    printf 'total_sites\t%s\n' "${total:-NA}"
    printf 'somatic_unstable_sites\t%s\n' "${somatic:-NA}"
    printf 'msi_score_pct\t%s\n' "$score"
    printf 'msi_cutoff_pct\t%s-%s\n' "$LOW" "$HIGH"
    printf 'msi_status\t%s\n' "$status"
} > "${SAMPLE}_msi_status.txt"
