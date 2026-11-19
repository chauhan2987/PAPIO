#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Somatic variant filters: tumor VAF (from AD), tumor alt reads, tumor depth, and alt
#          evidence in the matched normal. Called by the PAPIO process FILTER_SOMATIC_VCF.
# Usage  : filter_somatic_vcf.sh --vcf VCF --tumor ID --normal ID --mode tag|remove
#                                --min-tumor-vaf F --min-tumor-ad N --min-tumor-dp N
#                                --max-normal-vaf F --min-normal-alt N
#          A record fails when the tumor is below any tumor limit, or the normal has
#          VAF > max-normal-vaf AND >= min-normal-alt alt reads. tag: FILTER LOWSUPPORT;
#          remove: failing records are dropped. Samples are found by name (column order
#          differs between callers).
# Output : <tumor>_filtered.vcf.gz (+ .tbi), <tumor>_filter_stats.txt
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
VCF=""; TUMOR=""; NORMAL=""; MODE="tag"
MIN_T_VAF=""; MIN_T_AD=""; MIN_T_DP=""; MAX_N_VAF=""; MIN_N_ALT=""

##### Functions
## Function to print the 0-based column index of a sample in the VCF
sample_index() { bcftools query -l "$VCF" | awk -v s="$1" '$0 == s { print NR - 1; exit }'; }

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "filter_somatic_vcf.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --vcf)            VCF=$2 ;;
        --tumor)          TUMOR=$2 ;;
        --normal)         NORMAL=$2 ;;
        --mode)           MODE=$2 ;;
        --min-tumor-vaf)  MIN_T_VAF=$2 ;;
        --min-tumor-ad)   MIN_T_AD=$2 ;;
        --min-tumor-dp)   MIN_T_DP=$2 ;;
        --max-normal-vaf) MAX_N_VAF=$2 ;;
        --min-normal-alt) MIN_N_ALT=$2 ;;
        *) echo "filter_somatic_vcf.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done
OUT="${TUMOR}_filtered.vcf.gz"

ti=$(sample_index "$TUMOR"); ni=$(sample_index "$NORMAL")
if [ -z "$ti" ] || [ -z "$ni" ]; then
    echo "ERROR: could not find tumor '$TUMOR' / normal '$NORMAL' columns in $VCF" >&2
    echo "samples found: $(bcftools query -l "$VCF" | tr '\n' ' ')" >&2
    exit 1
fi

## VAF = alt / (ref + alt) from FORMAT/AD for every caller (Mutect2's FORMAT/AF is a model
## estimate, not a read fraction); every term is guarded against missing values ('.') and a
## zero denominator
t_vaf="FMT/AD[${ti}:1] / (FMT/AD[${ti}:0] + FMT/AD[${ti}:1])"
n_vaf="FMT/AD[${ni}:1] / (FMT/AD[${ni}:0] + FMT/AD[${ni}:1])"
t_ok="FMT/AD[${ti}:0]!='.' && FMT/AD[${ti}:1]!='.' && FMT/AD[${ti}:0] + FMT/AD[${ti}:1] > 0"
n_ok="FMT/AD[${ni}:0]!='.' && FMT/AD[${ni}:1]!='.' && FMT/AD[${ni}:0] + FMT/AD[${ni}:1] > 0"
EXPR="(${t_ok} && ${t_vaf} < ${MIN_T_VAF}) || (FMT/AD[${ti}:1]!='.' && FMT/AD[${ti}:1] < ${MIN_T_AD}) || (FMT/DP[${ti}]!='.' && FMT/DP[${ti}] < ${MIN_T_DP}) || (${n_ok} && ${n_vaf} > ${MAX_N_VAF} && FMT/AD[${ni}:1] >= ${MIN_N_ALT})"
if [ "$MODE" = "remove" ]; then
    bcftools view -e "$EXPR" -O z -o "$OUT" "$VCF"
else
    bcftools filter -m + -s LOWSUPPORT -e "$EXPR" -O z -o "$OUT" "$VCF"
fi
bcftools index -t "$OUT"

{   echo "FILTER SUMMARY: $TUMOR"
    echo "mode: $MODE"
    echo "thresholds: tumor_vaf>=${MIN_T_VAF} tumor_alt_ad>=${MIN_T_AD} tumor_dp>=${MIN_T_DP} normal_vaf<=${MAX_N_VAF} (fail only with normal_alt_ad>=${MIN_N_ALT})"
    echo "tumor_col_idx: $ti  normal_col_idx: $ni"
    echo "input variants:       $(bcftools view -H "$VCF" | wc -l)"
    echo "output total records: $(bcftools view -H "$OUT" | wc -l)"
    echo "output PASS records:  $(bcftools view -H -f PASS "$OUT" | wc -l)"
} > "${TUMOR}_filter_stats.txt"
