#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Split VarScan2 multi-allelic records and write standard AD (Number=R) / AF
#          fields. Called by the PAPIO process TAG_VARSCAN_AF.
# Usage  : varscan_tag_af.sh <sample_id> <ref.fa> <varscan_snv.vcf.gz> <varscan_indel.vcf.gz>
# Output : <sample>_varscan_snv.af.vcf.gz, <sample>_varscan_indel.af.vcf.gz (+ .tbi)
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=$1; REF=$2; SNV=$3; INDEL=$4

##### Functions
## Function to tag one VarScan VCF
tag() {
    local kind=$1 in=$2
    ## VarScan's single-value RD/AD cannot represent several ALT alleles: split first
    bcftools norm -m -any -f "$REF" -O v "$in" > "${kind}.split.vcf"
    ## declare AD as Number=R (pysam cannot re-add an existing FORMAT ID)
    sed 's/##FORMAT=<ID=AD,Number=[0-9A-Za-z.]*,Type=Integer/##FORMAT=<ID=AD,Number=R,Type=Integer/' \
        "${kind}.split.vcf" > "${kind}.hdrfix.vcf"
    varscan_add_af.py "${kind}.hdrfix.vcf" "${kind}.af.vcf"
    bgzip -c "${kind}.af.vcf" > "${SAMPLE}_varscan_${kind}.af.vcf.gz"
    bcftools index -t "${SAMPLE}_varscan_${kind}.af.vcf.gz"
}

##### Main Script
tag snv "$SNV"
tag indel "$INDEL"
