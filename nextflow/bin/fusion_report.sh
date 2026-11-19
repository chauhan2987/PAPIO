#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Consolidated RNA fusion report for one patient.
#          Called by the PAPIO process CONSOLIDATE_FUSION_REPORT.
# Usage  : fusion_report.sh <patient_id> <rna_sample_id> <consensus_stats.txt> <annotated.tsv>
# Output : <patient_id>_fusion_report.tsv
#------------------------------------------------------------------------------
set -uo pipefail

###### GLOBAL VARIABLES
PATIENT=$1; SAMPLE=$2; STATS=$3; ANNOT=$4
REPORT="${PATIENT}_fusion_report.tsv"

##### Main Script
{
    echo "## ============================================================"
    echo "## CONSOLIDATED RNA FUSION REPORT - PATIENT: ${PATIENT}   (RNA sample: ${SAMPLE})"
    echo "## Generated: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "## NOTE: computational summary for review; not a final clinical tier."
    echo "## ============================================================"

    printf '\n## [1] CALLING SUMMARY\n'
    if [ -s "$STATS" ]; then cat "$STATS"; else echo "(no stats)"; fi

    printf '\n## [2] PRIORITISED FUSIONS (KNOWN_ONCOGENIC / LIKELY_ONCOGENIC first)\n'
    if [ -s "$ANNOT" ]; then
        grep -E '^FusionName' "$ANNOT" || true
        grep -E 'KNOWN_ONCOGENIC|LIKELY_ONCOGENIC' "$ANNOT" || true
    else
        echo "(no fusion calls)"
    fi

    printf '\n## [3] ALL CONSENSUS FUSIONS (full annotated table)\n'
    if [ -s "$ANNOT" ]; then grep -v '^#' "$ANNOT"; else echo "(none)"; fi
} > "$REPORT"
exit 0
