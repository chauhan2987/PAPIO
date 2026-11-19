#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: COSMIC signature fitting of one sample's PASS SNVs; skipped below a minimum
#          SNV count, where a fit is not meaningful. Called by the PAPIO process
#          MUTATIONAL_SIGNATURES.
# Usage  : mutational_signatures.sh <sample_id> <vcf> <genome_build> <cosmic_version> <exome: True|False> <min_snv>
# Output : <sample>_signatures/ (or <sample>_signatures/INSUFFICIENT_MUTATIONS.txt)
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=$1; VCF=$2; GENOME=$3; COSMIC=$4; EXOME=$5; MIN_SNV=$6
OUT="${SAMPLE}_signatures"

##### Main Script
mkdir -p sig_input "$OUT"
## PASS only (the upstream filters run in 'tag' mode)
bcftools view -f PASS,. "$VCF" -O v -o "sig_input/${SAMPLE}.vcf"
n_snv=$(bcftools view -H -v snps "sig_input/${SAMPLE}.vcf" | wc -l)
echo "signatures: $n_snv SNVs (floor $MIN_SNV)"
if [ "$n_snv" -lt "$MIN_SNV" ]; then
    printf 'SKIPPED: %s SNVs (< %s) - signature fit not meaningful.\nn_snv\t%s\n' \
        "$n_snv" "$MIN_SNV" "$n_snv" > "${OUT}/INSUFFICIENT_MUTATIONS.txt"
    exit 0
fi
python3 "$(command -v fit_signatures.py)" sig_input "$OUT" "$GENOME" "$COSMIC" "$EXOME"
