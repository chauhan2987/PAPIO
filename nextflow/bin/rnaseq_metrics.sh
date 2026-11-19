#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: RNA-seq library strandedness and Picard CollectRnaSeqMetrics.
#          Called by the PAPIO process RNASEQ_METRICS.
# Usage  : rnaseq_metrics.sh --sample ID --bam BAM --strand auto|NONE|FIRST_READ_TRANSCRIPTION_STRAND|SECOND_READ_TRANSCRIPTION_STRAND
#                            --basis TEXT [--gene-counts ReadsPerGene.out.tab] [--refflat FILE]
#                            [--rrna FILE] --qc-filter "-f 2"|"-F 0x900" --threads N --java-opts "..."
#          --strand auto: inferred from the STAR gene counts (>= 80% of gene-assigned reads
#          on one strand), NONE when no counts are given. Without --refflat only the
#          strandedness is written.
# Output : <sample>.strandedness.txt, <sample>.rna_metrics
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; BAM=""; STRAND="NONE"; BASIS=""; COUNTS=""; REFFLAT=""; RRNA=""
QC_FILTER="-f 2"; THREADS=1; JAVA_OPTS=""

##### Functions
## Function to infer strandedness from STAR ReadsPerGene (col3 = read-1 sense, col4 = antisense)
infer_strand() {
    local fwd rev tot
    read -r fwd rev < <(awk -F'\t' '!/^N_/ { f += $3; r += $4 } END { print f + 0, r + 0 }' "$1")
    tot=$((fwd + rev))
    if   [ "$tot" -gt 0 ] && [ $((100 * fwd / tot)) -ge 80 ]; then STRAND="FIRST_READ_TRANSCRIPTION_STRAND"
    elif [ "$tot" -gt 0 ] && [ $((100 * rev / tot)) -ge 80 ]; then STRAND="SECOND_READ_TRANSCRIPTION_STRAND"
    else STRAND="NONE"; fi
    BASIS="inferred from STAR gene counts: sense=$fwd antisense=$rev"
}

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "rnaseq_metrics.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)      SAMPLE=$2 ;;
        --bam)         BAM=$2 ;;
        --strand)      STRAND=$2 ;;
        --basis)       BASIS=$2 ;;
        --gene-counts) COUNTS=$2 ;;
        --refflat)     REFFLAT=$2 ;;
        --rrna)        RRNA=$2 ;;
        --qc-filter)   QC_FILTER=$2 ;;
        --threads)     THREADS=$2 ;;
        --java-opts)   JAVA_OPTS=$2 ;;
        *) echo "rnaseq_metrics.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

if [ "$STRAND" = "auto" ]; then
    if [ -n "$COUNTS" ] && [ -s "$COUNTS" ]; then infer_strand "$COUNTS"; else STRAND="NONE"; fi
fi
printf 'sample_id\t%s\nstrand_specificity\t%s\nbasis\t%s\n' "$SAMPLE" "$STRAND" "$BASIS" > "${SAMPLE}.strandedness.txt"

if [ -z "$REFFLAT" ]; then
    echo "RNASEQ_METRICS: no refflat - CollectRnaSeqMetrics skipped (strandedness still reported)" >&2
    exit 0
fi

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/${SAMPLE}_rnaseqc_tmp.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT
## QC-only copy: properly paired reads for paired-end (Picard stops on interchromosomal
## discordant pairs), all primary reads for single-end. The callers use the raw BAM.
# shellcheck disable=SC2086
samtools view -@ "$THREADS" $QC_FILTER -b "$BAM" -o qc_input.bam
samtools index -@ "$THREADS" qc_input.bam
rrna_flag=()
[ -n "$RRNA" ] && rrna_flag=(--RIBOSOMAL_INTERVALS "$RRNA")
gatk --java-options "$JAVA_OPTS" CollectRnaSeqMetrics \
    -I qc_input.bam \
    -O "${SAMPLE}.rna_metrics" \
    --REF_FLAT "$REFFLAT" \
    --STRAND_SPECIFICITY "$STRAND" \
    "${rrna_flag[@]}" \
    --VALIDATION_STRINGENCY SILENT --TMP_DIR "$tmp_dir"
