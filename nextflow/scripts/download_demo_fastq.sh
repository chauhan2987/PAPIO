#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Date   : 2026-11-12
# Purpose: Download the PAPIO demonstration FASTQs from ENA, verify their MD5
#          checksums, and write ready-to-use samplesheets.
#
#   DNA (SEQC2 WES, Agilent SureSelect V6+UTR):
#     SRR7890883  HCC1395   tumour  -> sample_type DT
#     SRR7890874  HCC1395BL normal  -> sample_type DN
#   RNA (bulk RNA-seq):
#     SRR9134735  HCC1395            -> sample_type RN
#
# Usage  : bash download_demo_fastq.sh <out_dir> [patient_id] [rna_patient_id]
#          patient_id defaults to 'valid1' (tumour and normal share it, which is
#          how PAPIO pairs them: <patient_id>_T / _N). The RNA is a separate study
#          and gets its own patient_id, 'rna' by default (<rna_patient_id>_R).
#
# Output : <out_dir>/fastq/                   FASTQs (set params.fastq_root here)
#          <out_dir>/samplesheets/dna_somatic.csv   tumour + normal
#          <out_dir>/samplesheets/dna_germline.csv  normal only
#          <out_dir>/samplesheets/rna_fusion.csv    RNA
#
# Needs  : curl, wget, md5sum. Re-running is safe: finished files whose MD5
#          matches are skipped, partial files are resumed.
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
OUT_DIR=${1:?usage: download_demo_fastq.sh <out_dir> [patient_id]}
PATIENT_ID=${2:-valid1}
RNA_PATIENT_ID=${3:-rna}
FASTQ_DIR="${OUT_DIR}/fastq"
SHEET_DIR="${OUT_DIR}/samplesheets"
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
HEADER="patient_id,sex,tumor_type,batch,sample_type,run_num,flowcell,lane_num,chunk,read1,read2"   # as README "Samplesheets"

## run_accession  sample_type  sex (0 female)  tumor_type  batch
RUNS="SRR7890883 DT 0 BRCA IL_1
SRR7890874 DN 0 BRCA IL_1
SRR9134735 RN 0 BRCA IL_1"

mkdir -p "$FASTQ_DIR" "$SHEET_DIR"

#------------------------------------------------------------------------------
##### Functions

## Function to fetch "url<TAB>md5" lines (one per FASTQ) for a run from ENA
ena_files() {
    curl -sS --fail --retry 3 \
        "${ENA_API}?accession=$1&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv" \
    | awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
                  { n=split($h["fastq_ftp"],u,";"); split($h["fastq_md5"],m,";");
                    for (i=1;i<=n;i++) print u[i] "\t" m[i] }'
}

## Function to download one file and verify its MD5
fetch_verified() {
    local url=$1 md5=$2 dest
    dest="${FASTQ_DIR}/$(basename "$1")"
    if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then
        echo "  ok (already present): $(basename "$dest")"
        return 0
    fi
    wget -c -q --show-progress -O "$dest" "https://${url}"
    if [ "$(md5sum "$dest" | cut -d' ' -f1)" != "$md5" ]; then
        echo "ERROR: MD5 mismatch for $dest - delete it and re-run." >&2
        exit 1
    fi
    echo "  ok (downloaded + MD5 verified): $(basename "$dest")"
}

#------------------------------------------------------------------------------
##### Main Script

echo "$HEADER" > "${SHEET_DIR}/dna_somatic.csv"
echo "$HEADER" > "${SHEET_DIR}/dna_germline.csv"
echo "$HEADER" > "${SHEET_DIR}/rna_fusion.csv"

while read -r RUN TYPE SEX TUMOR BATCH; do
    echo "== ${RUN} (${TYPE})"
    FILES=$(ena_files "$RUN")
    [ -n "$FILES" ] || { echo "ERROR: ENA returned no FASTQ for ${RUN}" >&2; exit 1; }

    R1=""; R2=""
    while IFS=$'\t' read -r URL MD5; do
        fetch_verified "$URL" "$MD5"
        case "$URL" in
            *_1.fastq.gz) R1=$(basename "$URL") ;;
            *_2.fastq.gz) R2=$(basename "$URL") ;;
            *.fastq.gz)   [ -z "$R1" ] && R1=$(basename "$URL") ;;   # single-end
        esac
    done <<< "$FILES"

    ## DNA arm requires paired-end reads; RNA accepts single-end
    if [ "$TYPE" != "RN" ] && [ -z "$R2" ]; then
        echo "ERROR: ${RUN} is not paired-end; the DNA arm needs read1 + read2." >&2
        exit 1
    fi

    ## flowcell/run/lane are not in SRA metadata: the run accession is used as the
    ## flowcell so read groups stay unique; one lane, one chunk per run.
    ROW="${PATIENT_ID},${SEX},${TUMOR},${BATCH},${TYPE},${RUN},${RUN},1,001,${R1},${R2}"
    case "$TYPE" in
        DT) echo "$ROW" >> "${SHEET_DIR}/dna_somatic.csv" ;;
        DN) echo "$ROW" >> "${SHEET_DIR}/dna_somatic.csv"
            echo "$ROW" >> "${SHEET_DIR}/dna_germline.csv" ;;
        RN) echo "${RNA_PATIENT_ID}${ROW#"${PATIENT_ID}"}" >> "${SHEET_DIR}/rna_fusion.csv" ;;
    esac
done <<< "$RUNS"

echo
echo "Done."
echo "  params.fastq_root      : $(CDPATH= builtin cd "$FASTQ_DIR" >/dev/null && builtin pwd -P)"
echo "  params.dna_samplesheet : $(CDPATH= builtin cd "$SHEET_DIR" >/dev/null && builtin pwd -P)/dna_somatic.csv   (germline: dna_germline.csv)"
echo "  params.rna_samplesheet : $(CDPATH= builtin cd "$SHEET_DIR" >/dev/null && builtin pwd -P)/rna_fusion.csv"
