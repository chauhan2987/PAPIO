#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Per-gene coverage of reportable genes (length-weighted mean depth and the
#          fraction of bases at 10/20/30/50/100x), MAPQ-filtered.
#          Called by the PAPIO process MOSDEPTH when coverage_report_bed is set.
# Usage  : gene_coverage.sh --sample ID --bam BAM --bed REPORT_BED --mapq N --threads N
#          REPORT_BED column 4 = gene. When no region matches the BAM (chr1 vs 1
#          naming), the BED is retried once with toggled contig names.
# Output : <sample>.gene_coverage.tsv (+ mosdepth <sample>.genes.* files)
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; BAM=""; BED=""; MAPQ=20; THREADS=1

##### Functions
## Function to run mosdepth over the report regions
gene_md() {
    mosdepth --threads "$THREADS" --no-per-base --mapq "$MAPQ" \
        --by "$1" --thresholds 10,20,30,50,100 "${SAMPLE}.genes" "$BAM"
}

## Function to test whether mosdepth wrote any region (the pipe is read in full: pipefail-safe)
has_regions() { [ "$(zcat "${SAMPLE}.genes.regions.bed.gz" | wc -l)" -gt 0 ]; }

## Function to write a copy of a BED with chr1 <-> 1 (and chrM <-> MT) contig names
toggle_chr() {
    awk -F'\t' 'BEGIN { OFS = "\t" } /^(#|track|browser)/ { next }
        { c = $1; d = (c ~ /^chr/) ? substr(c, 4) : "chr" c
          if (d == "chrMT") d = "chrM"; else if (d == "M") d = "MT"; $1 = d; print }' "$1"
}

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "gene_coverage.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)  SAMPLE=$2 ;;
        --bam)     BAM=$2 ;;
        --bed)     BED=$2 ;;
        --mapq)    MAPQ=$2 ;;
        --threads) THREADS=$2 ;;
        *) echo "gene_coverage.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done
OUT="${SAMPLE}.gene_coverage.tsv"

gene_md "$BED"
if ! has_regions; then
    echo "WARN: coverage_report_bed matched no BAM contig - retrying with toggled chr naming" >&2
    toggle_chr "$BED" > report.toggled.bed
    gene_md report.toggled.bed
fi
if ! has_regions; then
    echo "WARN: coverage_report_bed shares no contig with the BAM - gene coverage NOT computed" >&2
    printf '# ERROR: coverage_report_bed shares no contig with the BAM (check chr naming / build)\n' > "$OUT"
    exit 0
fi

## length-weighted mean depth per gene, and threshold base counts per gene
zcat "${SAMPLE}.genes.regions.bed.gz" \
    | awk -F'\t' '{ L = $3 - $2; bp[$4] += L; dsum[$4] += L * $5 }
                  END { for (g in bp) printf "%s\t%d\t%.2f\n", g, bp[g], dsum[g] / bp[g] }' \
    | LC_ALL=C sort -k1,1 > depth.tmp
zcat "${SAMPLE}.genes.thresholds.bed.gz" \
    | awk -F'\t' '!/^#/ { g = $4; for (i = 5; i <= 9; i++) t[g, i] += $i }
                  END { for (k in t) { split(k, a, SUBSEP); seen[a[1]] = 1 }
                        for (g in seen) printf "%s\t%d\t%d\t%d\t%d\t%d\n", g, t[g,5], t[g,6], t[g,7], t[g,8], t[g,9] }' \
    | LC_ALL=C sort -k1,1 > thr.tmp
{   printf 'gene\ttarget_bp\tmean_depth\tpct_10x\tpct_20x\tpct_30x\tpct_50x\tpct_100x\n'
    LC_ALL=C join -t $'\t' depth.tmp thr.tmp \
        | awk -F'\t' 'BEGIN { OFS = "\t" }
            { L = $2; printf "%s\t%d\t%s", $1, L, $3
              for (i = 4; i <= 8; i++) printf "\t%.4f", (L > 0 ? $i / L : 0); printf "\n" }'
} > "$OUT"
rm -f depth.tmp thr.tmp
