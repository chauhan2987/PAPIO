#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Sequenza purity/ploidy/allele-specific copy number for one tumor-normal pair
#          (bam2seqz per chromosome, merge, binning, R fit, fit-ambiguity QC).
#          Called by the PAPIO process SEQUENZA.
# Usage  : sequenza.sh --sample ID --tumor BAM --normal BAM --fasta REF --fai REF.fai
#                      --gc-wiggle FILE --parallel N --window N --gamma N --kmin N
#                      --female TRUE|FALSE|NA [--female-min-xratio F --male-max-xratio F]
#                      --segment-filter BP --assembly hg38|hg19 [--chroms "chr1 chr2 ..."]
#          Without --chroms the primary chromosomes (1-22, X, Y) are read from the .fai.
#          --female NA: sex from the normal's chrX / autosome depth in the seqz data
#          (>= female-min-xratio female, <= male-max-xratio male); when still unclear the
#          male model is used, which keeps chrX/chrY out of the purity/ploidy fit.
# Output : <sample>_segments.txt, _confints_CP.txt, _alternative_solutions.txt, plots,
#          RData, _mutations.txt, _sequenza_log.txt, _fit_qc.tsv
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; TUMOR=""; NORMAL=""; FASTA=""; FAI=""; GC=""; SLOTS=4; WINDOW=50
GAMMA=80; KMIN=10; FEMALE=FALSE; SEG_FILT=3000000; ASSEMBLY=hg38; CHROMS=""
FEMALE_MIN=0.80; MALE_MAX=0.65; SEX_NOTE=""

##### Functions
## Function to remove intermediates and plots that are not published
cleanup() {
    rm -f "${SAMPLE}".seqz.gz* "${SAMPLE}"_*.seqz.gz "${SAMPLE}".merged.seqz.gz "${SAMPLE}".binned.seqz.gz \
          "${SAMPLE}_CN_bars.pdf" "${SAMPLE}_chromosome_depths.pdf" "${SAMPLE}_gc_plots.pdf" \
          "${SAMPLE}_model_fit.pdf" "${SAMPLE}_alternative_fit.pdf" fit_top2.tmp
}

## Function to merge per-chromosome seqz files in chromosome order under one header
merge_seqz() {
    local first=1 chr f
    for chr in $CHROMS; do
        f="${SAMPLE}_${chr}.seqz.gz"
        if [ ! -s "$f" ]; then echo "WARN: no seqz for ${chr} ($f) - skipping" >&2; continue; fi
        if [ "$first" = 1 ]; then zcat "$f"; first=0
        else zcat "$f" | awk 'NR > 1 || $1 != "chromosome"'; fi
    done | bgzip > "${SAMPLE}.merged.seqz.gz"
}

## Function to set FEMALE from the normal depth of chrX vs the autosomes in a seqz file
infer_female() {
    local ratio
    ratio=$(zcat "$1" | awk -F'\t' 'NR > 1 && $4 > 0 {
                if ($1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/) { a += $4; na++ }
                else if ($1 ~ /^(chr)?X$/)               { x += $4; nx++ } }
            END { if (na > 0 && nx > 0 && a > 0) printf "%.3f", (x / nx) / (a / na); else print "NA" }')
    if   [ "$ratio" = "NA" ]; then FEMALE=FALSE; SEX_NOTE="sex NA, no chrX data: male model (chrX/chrY not fitted)"
    elif awk -v r="$ratio" -v t="$FEMALE_MIN" 'BEGIN { exit !(r + 0 >= t + 0) }'; then FEMALE=TRUE; SEX_NOTE="sex NA, inferred female (normal chrX/autosome depth $ratio)"
    elif awk -v r="$ratio" -v t="$MALE_MAX"   'BEGIN { exit !(r + 0 <= t + 0) }'; then FEMALE=FALSE; SEX_NOTE="sex NA, inferred male (normal chrX/autosome depth $ratio)"
    else FEMALE=FALSE; SEX_NOTE="sex NA, ambiguous chrX/autosome depth $ratio: male model (chrX/chrY not fitted)"; fi
    echo "SEQUENZA: $SEX_NOTE" >&2
}

## Function to write the top-2 purity/ploidy solutions by SLPP and their difference
## (reported, not judged: close solutions are common in hypodiploid or low-purity tumors)
fit_qc() {
    awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i
                          if (!("cellularity" in h) || !("ploidy" in h) || !("SLPP" in h)) exit
                          next }
                { print $(h["cellularity"]) "\t" $(h["ploidy"]) "\t" $(h["SLPP"]) }' \
        "${SAMPLE}_alternative_solutions.txt" | LC_ALL=C sort -t$'\t' -k3,3gr | awk 'NR <= 2' > fit_top2.tmp
    {   printf 'rank\tcellularity\tploidy\tSLPP\n'
        awk -F'\t' '{ print NR "\t" $0 }' fit_top2.tmp
        awk -F'\t' 'NR == 1 { c1 = $1; p1 = $2; l1 = $3 } NR == 2 { c2 = $1; p2 = $2; l2 = $3 }
                    END { if (NR >= 2) printf "delta_1_minus_2\t%.3f\t%.2f\t%.3f\n", c1 - c2, p1 - p2, l1 - l2
                          else print "delta_1_minus_2\tNA\tNA\tNA" }' fit_top2.tmp
    } > "${SAMPLE}_fit_qc.tsv"
}

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "sequenza.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)         SAMPLE=$2 ;;
        --tumor)          TUMOR=$2 ;;
        --normal)         NORMAL=$2 ;;
        --fasta)          FASTA=$2 ;;
        --fai)            FAI=$2 ;;
        --gc-wiggle)      GC=$2 ;;
        --parallel)       SLOTS=$2 ;;
        --window)         WINDOW=$2 ;;
        --gamma)          GAMMA=$2 ;;
        --kmin)           KMIN=$2 ;;
        --female)         FEMALE=$2 ;;
        --female-min-xratio) FEMALE_MIN=$2 ;;
        --male-max-xratio)   MALE_MAX=$2 ;;
        --segment-filter) SEG_FILT=$2 ;;
        --assembly)       ASSEMBLY=$2 ;;
        --chroms)         CHROMS=$2 ;;
        *) echo "sequenza.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done
trap cleanup EXIT

if [ -z "$CHROMS" ]; then
    CHROMS=$(awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/ { printf "%s ", $1 }' "$FAI")
fi
[ -n "$CHROMS" ] || { echo "ERROR: no primary chromosomes (1-22, X, Y) in $FAI" >&2; exit 1; }

## 1. raw seqz, one file per chromosome (<prefix>_<chrom>.seqz.gz); -C skips decoy/viral/HLA contigs
# shellcheck disable=SC2086
sequenza-utils bam2seqz --normal "$NORMAL" --tumor "$TUMOR" --fasta "$FASTA" -gc "$GC" \
    -C $CHROMS --parallel "$SLOTS" --output "${SAMPLE}.seqz.gz" -q 20 -N 30

## 2. merge; an empty merge would end in a ZeroDivisionError during binning
merge_seqz
if [ "$(zcat "${SAMPLE}.merged.seqz.gz" | awk 'NR <= 2' | wc -l)" -lt 2 ]; then
    echo "ERROR: merged seqz empty - no per-chromosome files matched ${SAMPLE}_<chr>.seqz.gz" >&2
    ls -la "${SAMPLE}"*.seqz.gz >&2 || true
    exit 1
fi

## 3. bin to cut the memory of the R fit
sequenza-utils seqz_binning --seqz "${SAMPLE}.merged.seqz.gz" --window "$WINDOW" -o "${SAMPLE}.binned.seqz.gz"

## sex for the model when the samplesheet has none
[ "$FEMALE" = "NA" ] && infer_female "${SAMPLE}.binned.seqz.gz"

## 4. fit (bin/run_sequenza.R) and fit-ambiguity QC
Rscript --vanilla "$(command -v run_sequenza.R)" "${SAMPLE}.binned.seqz.gz" "$SAMPLE" \
    "$ASSEMBLY" "$GAMMA" "$KMIN" "$FEMALE" "$SEG_FILT"
fit_qc
[ -z "$SEX_NOTE" ] || echo "sex_model: $SEX_NOTE" >> "${SAMPLE}_sequenza_log.txt"
