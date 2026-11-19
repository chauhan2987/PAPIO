#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Date   : 2026-11-17
# Purpose: Build the slim gnomAD v4.1 joint (exome + genome) VCF used by PAPIO
#          as the population frequency for variant classification
#          (params.gnomad_faf_vcf). Keeps one INFO field, fafmax_faf95_max_joint
#          (grpmax filtering allele frequency, the field ClinGen recommends for
#          ACMG/AMP BA1/BS1), on PASS records only.
#          Run once per gnomAD release; the output is reused across runs.
#
# Usage  : bash prepare_gnomad_faf.sh <out_dir> [threads]
# Output : <out_dir>/gnomad.joint.v4.1.faf95.vcf.gz (+ .tbi)
#
# Note   : downloads the large per-chromosome VCFs one at a time (each is
#          deleted after slimming). Needs bcftools >= 1.15, curl, tabix.
#          Set GNOMAD_SRC to a local copy of the joint VCFs to skip downloading.
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
OUT_DIR=${1:?usage: prepare_gnomad_faf.sh <out_dir> [threads]}
THREADS=${2:-4}
RELEASE="4.1"
URL_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/${RELEASE}/vcf/joint"
GNOMAD_SRC=${GNOMAD_SRC:-}
FIELD="fafmax_faf95_max_joint"
CHROMS=$(echo chr{1..22} chrX chrY)
OUT="${OUT_DIR}/gnomad.joint.v${RELEASE}.faf95.vcf.gz"

#------------------------------------------------------------------------------
##### Functions

## Function to get one chromosome's joint VCF (local copy or download)
fetch_chrom() {
    local chr=$1 name="gnomad.joint.v${RELEASE}.sites.${1}.vcf.bgz"
    if [ -n "$GNOMAD_SRC" ]; then
        echo "${GNOMAD_SRC}/${name}"
    else
        curl -fsSL --retry 3 -o "${OUT_DIR}/tmp/${name}" "${URL_BASE}/${name}"
        echo "${OUT_DIR}/tmp/${name}"
    fi
}

## Function to keep PASS records and the FAF field only
slim_chrom() {
    local src=$1 chr=$2
    bcftools view --threads "$THREADS" -f PASS -Ou "$src" \
    | bcftools annotate --threads "$THREADS" -x "^INFO/${FIELD}" -Ob \
        -o "${OUT_DIR}/tmp/${chr}.bcf"
    bcftools index -f "${OUT_DIR}/tmp/${chr}.bcf"
}

#------------------------------------------------------------------------------
##### Main Script
mkdir -p "${OUT_DIR}/tmp"
for chr in $CHROMS; do
    echo "[$(date '+%F %T')] ${chr}"
    src=$(fetch_chrom "$chr")
    slim_chrom "$src" "$chr"
    if [ -z "$GNOMAD_SRC" ]; then rm -f "$src"; fi
done

bcftools concat --threads "$THREADS" -Oz -o "$OUT" \
    $(for chr in $CHROMS; do echo "${OUT_DIR}/tmp/${chr}.bcf"; done)
tabix -p vcf "$OUT"
rm -rf "${OUT_DIR}/tmp"

n=$(bcftools index -n "$OUT")
echo "Done: ${OUT} (${n} records). Set gnomad_faf_vcf: ${OUT}"
