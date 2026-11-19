#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Fusion consensus across callers. Builds the gene-ID tables and the reference
#          chromosome naming that fusion_consensus.py uses, then runs it.
#          Called by the PAPIO process FUSION_CONSENSUS.
# Usage  : fusion_consensus.sh <sample_id> <min_callers> <readthrough_kb> <known_fusions|NONE>
#                              <rna_gtf> [<ctat_ref_annot.gtf>]
#          Caller files are read from the working directory (arriba.tsv, starfusion.tsv,
#          fusioncatcher.txt; a disabled caller's file is absent).
# Output : <sample>_fusion_consensus.tsv, _fusion_consensus_stats.txt, _fusion_candidates.txt,
#          _fusion_names.tsv
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=$1; MIN_CALLERS=$2; READTHRU_KB=$3; KNOWN=$4; GTF=$5; CTAT_GTF=${6:-}

##### Functions
## Function to write gene_id <TAB> gene_name for every gene of a GTF (plain or gzipped)
gene_table() {
    zcat -f "$1" | awk -F'\t' '$3 == "gene" || $3 == "transcript" || $3 == "exon" {
        if (match($9, /gene_id "[^"]+"/)) { id = substr($9, RSTART + 9, RLENGTH - 10)
            if (match($9, /gene_name "[^"]+"/)) print id "\t" substr($9, RSTART + 11, RLENGTH - 12) } }' \
        | LC_ALL=C sort -u
}

##### Main Script
gene_table "$GTF" > genes_main.tsv
if [ -n "$CTAT_GTF" ] && [ -s "$CTAT_GTF" ]; then gene_table "$CTAT_GTF" > genes_ctat.tsv; else cp genes_main.tsv genes_ctat.tsv; fi

## reference chromosome naming (chr1 or 1); FusionCatcher breakpoints are written the same way.
## awk stops at the first record: '|| true' keeps zcat's SIGPIPE from failing under pipefail.
chr_style=$({ zcat -f "$GTF" || true; } | awk -F'\t' '!/^#/ { print (($1 ~ /^chr/) ? "chr" : "plain"); exit }')

PY=$(command -v python3 || echo /opt/conda/bin/python3)
"$PY" "$(command -v fusion_consensus.py)" "$MIN_CALLERS" "$READTHRU_KB" "$KNOWN" "$SAMPLE" \
    genes_main.tsv genes_ctat.tsv "${chr_style:-asis}"
## never hand FusionInspector a missing candidate list
[ -e "${SAMPLE}_fusion_candidates.txt" ] || : > "${SAMPLE}_fusion_candidates.txt"
