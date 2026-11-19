#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Keep FusionInspector candidates whose two genes are in the CTAT annotation,
#          written in the GTF's own case (FusionInspector looks names up exactly).
#          Called by the PAPIO process FILTER_FUSION_CANDIDATES.
# Usage  : filter_fusion_candidates.sh <candidates.txt> <ctat_genome_lib_dir> <out.txt>
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
CANDIDATES=$1; CTAT=$2; OUT=$3

##### Main Script
if [ ! -s "$CANDIDATES" ]; then : > "$OUT"; exit 0; fi

gtf=$(find -L "$CTAT" -name ref_annot.gtf -print -quit)
[ -n "$gtf" ] || { echo "ERROR: ref_annot.gtf not found in $CTAT" >&2; exit 1; }
awk -F'\t' '$3 == "gene" || $3 == "exon" { print $9 }' "$gtf" \
    | grep -oE '(gene_id|gene_name|gene_symbol) "[^"]+"' \
    | cut -d'"' -f2 | sort -u > valid_genes.txt

## case-insensitive match: FUSION_CONSENSUS upper-cases symbols, the GTF keeps mixed case
awk -F'\t' -v gtf="valid_genes.txt" '
    BEGIN { while ((getline line < gtf) > 0) { u = toupper(line); if (!(u in valid)) valid[u] = line } }
    /^#/ { print; next }
    {
        n = split($1, genes, "--")
        if (n > 1) {
            a = toupper(genes[1]); b = toupper(genes[2])
            if ((a in valid) && (b in valid)) print valid[a] "--" valid[b]
            else print "WARN: dropped (gene not in CTAT GTF): " $1 > "/dev/stderr"
        } else if (toupper($1) in valid) {
            print valid[toupper($1)]
        }
    }' "$CANDIDATES" > "$OUT"
