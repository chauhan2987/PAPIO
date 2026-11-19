#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: CNVkit copy number for one tumor. Called by the PAPIO process CNVKIT.
# Usage  : run_cnvkit.sh --sample ID --tumor BAM --tumor-bai BAI --assay wes|wgs --threads N
#                        --fasta REF [--normal BAM --normal-bai BAI] [--pon reference.cnn] [--intervals BED]
#          With a panel of normals only --reference is used (it already holds the bins);
#          without one, a reference is built from the matched normal: WES bins from the
#          target BED (--targets), WGS bins genome-wide (--method wgs, --access = BED).
#          Tumor-only needs a panel of normals; a WGS run refuses an exome-style panel.
# Output : <sample>.cnr, <sample>.cns (+ CNVkit's other files)
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; TUMOR=""; TUMOR_BAI=""; NORMAL=""; NORMAL_BAI=""; PON=""; INTERVALS=""; ASSAY="wes"; THREADS=1; FASTA=""

##### Functions
## Function to test whether a CNVkit reference has antitarget bins (exome-style).
## awk stops at the first hit; '|| true' keeps the reader's SIGPIPE from turning a hit into a miss.
has_antitargets() {
    local reader=cat
    case "$1" in *.gz) reader="zcat" ;; esac
    { "$reader" "$1" || true; } | awk 'NR > 1 && $4 ~ /[Aa]ntitarget/ { f = 1; exit } END { exit !f }'
}

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "run_cnvkit.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)    SAMPLE=$2 ;;
        --tumor)     TUMOR=$2 ;;
        --tumor-bai) TUMOR_BAI=$2 ;;
        --normal)    NORMAL=$2 ;;
        --normal-bai) NORMAL_BAI=$2 ;;
        --pon)       PON=$2 ;;
        --intervals) INTERVALS=$2 ;;
        --assay)     ASSAY=$2 ;;
        --threads)   THREADS=$2 ;;
        --fasta)     FASTA=$2 ;;
        *) echo "run_cnvkit.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

## guards: fail loud rather than emit silently wrong copy number
if [ -z "$NORMAL" ] && [ -z "$PON" ]; then
    echo "ERROR: CNVkit tumor-only requires a panel-of-normals (cnvkit_pon) as --reference." >&2; exit 1
fi
if [ -z "$PON" ] && [ "$ASSAY" != "wgs" ] && { [ -z "$INTERVALS" ] || [ ! -s "$INTERVALS" ]; }; then
    echo "ERROR: assay=wes without a PoN requires a target/interval BED (--targets)." >&2; exit 1
fi
if [ "$ASSAY" = "wgs" ] && [ -n "$PON" ] && has_antitargets "$PON"; then
    echo "ERROR: cnvkit_pon has Antitarget bins (exome-style) but assay=wgs." >&2
    echo "       Build a WGS PoN (cnvkit.py batch --method wgs) or run reference-free." >&2
    exit 1
fi

## fixed input names: outputs do not depend on dedup_mode (_markdup / _consensus / _bqsr)
ln -sf "$TUMOR" tumor.bam
ln -sf "$TUMOR_BAI" tumor.bam.bai

if [ -n "$PON" ]; then
    args=(--reference "$PON")
else
    args=()
    if [ -n "$NORMAL" ]; then
        ln -sf "$NORMAL" normal.bam
        ln -sf "$NORMAL_BAI" normal.bam.bai
        args+=(--normal normal.bam)
    fi
    if [ "$ASSAY" = "wgs" ]; then
        args+=(--method wgs)
        [ -n "$INTERVALS" ] && args+=(--access "$INTERVALS")
    else
        args+=(--targets "$INTERVALS")
    fi
    args+=(--fasta "$FASTA")
fi

## --drop-low-coverage: CNVkit's advice for tumors (no false deep deletions in poorly covered bins)
cnvkit.py batch tumor.bam "${args[@]}" --drop-low-coverage --processes "$THREADS" --output-dir .
mv tumor.cns "${SAMPLE}.cns"
mv tumor.cnr "${SAMPLE}.cnr"
