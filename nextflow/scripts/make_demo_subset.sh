#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Date   : 2026-11-12
# Purpose: Build the PAPIO demonstration dataset from the full benchmark run.
#
#   The demo must be small enough for Zenodo or a GitHub release (a few tens of MB).
#   DNA : picks the SMALLEST chromosome (by capture-target bp) that still has
#         >= MIN_SNV truth SNVs and >= MIN_INDEL truth indels inside the
#         evaluation regions (targets x SEQC2 high-confidence). Only the capture
#         targets that hold one of those truth variants, or one of N_MSI exonic
#         microsatellites (so MSI still runs), are kept (+/- WINDOW_PAD bp), and
#         every read pair (tumour + normal) touching them.
#   RNA : read pairs touching the partner genes (+/- FLANK_BP) of the reference
#         fusion list: every chimeric pair (SA tag, mate on another contig) is kept,
#         the others are subsampled to MAX_RNA_PAIRS (fixed seed).
#
#   Read names come from the full run's BAMs; the reads themselves are taken
#   from the ORIGINAL FASTQs (untrimmed, both mates, original order), so the
#   demo re-runs the pipeline from 'mapping' exactly like the full data.
#
# Usage  : edit GLOBAL VARIABLES, then  bash make_demo_subset.sh
# Needs  : samtools, bcftools, awk, gzip  (e.g. inside the align container)
# Output : ${DEMO_DIR}/fastq/*.fastq.gz
#          ${DEMO_DIR}/demo_targets.bed     kept capture targets (test-bundle targets)
#          ${DEMO_DIR}/demo_windows.bed     kept targets +/- WINDOW_PAD, merged (resources are cut to these)
#          ${DEMO_DIR}/samplesheets/{dna_somatic,dna_germline,rna_fusion}.csv
#          ${DEMO_DIR}/demo_selection.tsv   per-chromosome table + chosen chromosome
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
DEMO_DIR=/path/to/demo_subset
PATIENT_ID=valid1

## full-run inputs (DNA)
TARGET_BED=/path/to/papio_resources/capture/S07604624/S07604624_Covered.bed
HIGHCONF_BED=/path/to/seqc2/High-Confidence_Regions_v1.2.bed
TRUTH_SNV=/path/to/seqc2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
TRUTH_INDEL=/path/to/seqc2/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz
TUMOR_BAM=/path/to/run/Result/02_markdup/valid1_T_markdup.bam      # or .cram (set REF_FASTA)
NORMAL_BAM=/path/to/run/Result/02_markdup/valid1_N_markdup.bam
TUMOR_R1=/path/to/demo/fastq/SRR7890945_1.fastq.gz
TUMOR_R2=/path/to/demo/fastq/SRR7890945_2.fastq.gz
NORMAL_R1=/path/to/demo/fastq/SRR7890874_1.fastq.gz
NORMAL_R2=/path/to/demo/fastq/SRR7890874_2.fastq.gz
REF_FASTA=/path/to/papio_resources/genome/GRCh38.d1.vd1.fa          # only needed for CRAM
MSI_LIST=/path/to/papio_resources/msi/GRCh38.d1.vd1.msisensor.list   # MSIsensor-pro scan of the genome

## full-run inputs (RNA)
RNA_BAM=/path/to/run_rna/Result/01_mapping/valid1_R.Aligned.sortedByCoord.out.bam
RNA_R1=/path/to/demo/fastq/SRR9134735_1.fastq.gz
RNA_R2=/path/to/demo/fastq/SRR9134735_2.fastq.gz
CTAT_GTF=/path/to/ctat/ctat_genome_lib_build_dir/ref_annot.gtf
REF_FUSIONS=/path/to/benchmark/reference_fusions.tsv   # 2 columns: GENE5<TAB>GENE3 (header lines starting with # ignored)

## selection thresholds
CANDIDATE_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
MIN_SNV=50          # truth SNVs needed on the demo chromosome
MIN_INDEL=5         # truth indels needed on the demo chromosome
FLANK_BP=10000      # RNA: padding around each fusion partner gene
FORCE_CHROM=""      # set (e.g. chr22) to skip the automatic choice
N_MSI=100           # exonic microsatellites kept (MSIsensor-pro needs sites in the targets)
WINDOW_PAD=200      # bp around each kept target
MAX_RNA_PAIRS=100000   # non-chimeric RNA read pairs kept (all chimeric pairs are kept)

mkdir -p "${DEMO_DIR}/fastq" "${DEMO_DIR}/samplesheets" "${DEMO_DIR}/tmp"
TMP="${DEMO_DIR}/tmp"

#------------------------------------------------------------------------------
##### Functions

## Function to collect unique read names with >= 1 primary record in regions ($2 = region args)
read_names() {
    local bam=$1; shift
    samtools view -F 0x900 ${REF_FASTA:+-T "$REF_FASTA"} "$bam" "$@" | cut -f1 | LC_ALL=C sort -u
}

## Function to keep FASTQ records whose name is in a name list (works per mate; order preserved)
subset_fastq() {
    local names=$1 in_fq=$2 out_fq=$3
    zcat "$in_fq" \
    | awk -v NAMES="$names" '
        BEGIN { while ((getline l < NAMES) > 0) keep[l] = 1 }
        NR % 4 == 1 { h = substr($1, 2); sub(/\/[12]$/, "", h); p = (h in keep) }
        p' \
    | gzip > "$out_fq"
}

## Function to count FASTQ records
n_reads() { echo $(( $(zcat "$1" | wc -l) / 4 )); }

#------------------------------------------------------------------------------
##### Main Script

## ---- 1. per-chromosome table: target bp + truth variants in evaluation regions ----
echo "== Scoring chromosomes"
awk '{bp[$1] += $3 - $2} END {for (c in bp) print c "\t" bp[c]}' "$TARGET_BED" | LC_ALL=C sort -k1,1 > "${TMP}/target_bp.tsv"
TOTAL_BP=$(awk '{s += $2} END {print s}' "${TMP}/target_bp.tsv")

## truth -> inside targets -> inside high-confidence regions (streamed, no bedtools)
for V in SNV INDEL; do
    [ "$V" = "SNV" ] && F=$TRUTH_SNV || F=$TRUTH_INDEL
    bcftools view -T "$TARGET_BED" -Ou "$F" \
    | bcftools view -H -T "$HIGHCONF_BED" - \
    | cut -f1 | LC_ALL=C sort | uniq -c | awk '{print $2 "\t" $1}' > "${TMP}/truth_${V}.tsv"
done

printf 'chrom\ttarget_bp\tpct_of_exome\ttruth_snv\ttruth_indel\n' > "${DEMO_DIR}/demo_selection.tsv"
for C in $CANDIDATE_CHROMS; do
    BP=$(awk -v c="$C" '$1==c {print $2}' "${TMP}/target_bp.tsv"); BP=${BP:-0}
    S=$(awk -v c="$C" '$1==c {print $2}' "${TMP}/truth_SNV.tsv");   S=${S:-0}
    I=$(awk -v c="$C" '$1==c {print $2}' "${TMP}/truth_INDEL.tsv"); I=${I:-0}
    printf '%s\t%s\t%.2f\t%s\t%s\n' "$C" "$BP" "$(awk -v a="$BP" -v b="$TOTAL_BP" 'BEGIN{print 100*a/b}')" "$S" "$I" >> "${DEMO_DIR}/demo_selection.tsv"
done

## smallest target footprint that meets both thresholds
if [ -n "$FORCE_CHROM" ]; then
    CHROM=$FORCE_CHROM
else
    CHROM=$(tail -n +2 "${DEMO_DIR}/demo_selection.tsv" | sort -t$'\t' -k2,2n \
            | awk -F'\t' -v s="$MIN_SNV" -v i="$MIN_INDEL" '$4 >= s && $5 >= i {print $1; exit}')
fi
if [ -z "$CHROM" ]; then
    echo "No single chromosome has >= ${MIN_SNV} SNVs and >= ${MIN_INDEL} indels. Table:" >&2
    cat "${DEMO_DIR}/demo_selection.tsv" >&2
    echo "Lower MIN_SNV / MIN_INDEL or set FORCE_CHROM." >&2
    exit 1
fi
printf '# chosen\t%s\n' "$CHROM" >> "${DEMO_DIR}/demo_selection.tsv"
echo "  chosen: ${CHROM}"; cat "${DEMO_DIR}/demo_selection.tsv"

## ---- 2. DNA subset: targets holding a truth variant or an exonic microsatellite ----
echo "== DNA subset (${CHROM}: targets with truth variants or microsatellites)"
awk -F'\t' -v c="$CHROM" 'BEGIN { OFS = "\t" } $1 == c && !/^(#|track|browser)/ { print $1, $2, $3 }' "$TARGET_BED" \
    | LC_ALL=C sort -k1,1 -k2,2n > "${TMP}/targets_chrom.bed"
## points (0-based start, end): truth SNVs/indels in targets x high-confidence, and microsatellites in targets
for F in "$TRUTH_SNV" "$TRUTH_INDEL"; do
    bcftools view -r "$CHROM" -T "${TMP}/targets_chrom.bed" -Ou "$F" | bcftools view -H -T "$HIGHCONF_BED" - \
        | awk -F'\t' 'BEGIN { OFS = "\t" } { print $1, $2 - 1, $2 }'
done > "${TMP}/points.bed"
awk -F'\t' -v c="$CHROM" 'NR > 1 && $1 == c { print $1 "\t" $2 - 1 "\t" $2 }' "$MSI_LIST" | LC_ALL=C sort -k2,2n > "${TMP}/msi_all.bed"
## microsatellites inside targets (sorted two-pointer sweep), then N_MSI of them spread along the chromosome
awk -F'\t' 'NR == FNR { n++; ts[n] = $2; te[n] = $3; next } { while (i < n && te[i + 1] <= $2) i++; if (i < n && ts[i + 1] <= $2 && $3 <= te[i + 1]) print }' \
    "${TMP}/targets_chrom.bed" "${TMP}/msi_all.bed" > "${TMP}/msi_in_targets.bed"
N_IN=$(awk 'END { print NR }' "${TMP}/msi_in_targets.bed")
awk -v n="$N_IN" -v k="$N_MSI" 'n <= k || (NR - 1) % int(n / k) == 0' "${TMP}/msi_in_targets.bed" | awk -v k="$N_MSI" 'NR <= k' >> "${TMP}/points.bed"
LC_ALL=C sort -k1,1 -k2,2n -o "${TMP}/points.bed" "${TMP}/points.bed"
## targets holding at least one point
awk -F'\t' 'NR == FNR { n++; ps[n] = $2; pe[n] = $3; next }
    { while (j < n && pe[j + 1] <= $2) j++
      if (j < n && ps[j + 1] < $3) print }' "${TMP}/points.bed" "${TMP}/targets_chrom.bed" > "${DEMO_DIR}/demo_targets.bed"
awk -v p="$WINDOW_PAD" 'BEGIN { OFS = "\t" } { s = $2 - p; if (s < 0) s = 0; print $1, s, $3 + p }' "${DEMO_DIR}/demo_targets.bed" \
    | awk 'BEGIN { OFS = "\t" } { if ($1 == c && $2 <= e) { if ($3 > e) e = $3 } else { if (c != "") print c, s, e; c = $1; s = $2; e = $3 } }
           END { if (c != "") print c, s, e }' > "${DEMO_DIR}/demo_windows.bed"
echo "  $(awk 'END { print NR }' "${TMP}/points.bed") points ($(awk 'END { print NR }' "${TMP}/msi_in_targets.bed") microsatellites in targets), $(awk 'END { print NR }' "${DEMO_DIR}/demo_targets.bed") targets, windows $(awk '{ s += $3 - $2 } END { print s }' "${DEMO_DIR}/demo_windows.bed") bp"
[ -s "${DEMO_DIR}/demo_targets.bed" ] || { echo "ERROR: no target holds a truth variant or microsatellite" >&2; exit 1; }
read_names "$TUMOR_BAM"  -L "${DEMO_DIR}/demo_windows.bed" > "${TMP}/tumor.names"
read_names "$NORMAL_BAM" -L "${DEMO_DIR}/demo_windows.bed" > "${TMP}/normal.names"
subset_fastq "${TMP}/tumor.names"  "$TUMOR_R1"  "${DEMO_DIR}/fastq/demo_T_${CHROM}_1.fastq.gz"
subset_fastq "${TMP}/tumor.names"  "$TUMOR_R2"  "${DEMO_DIR}/fastq/demo_T_${CHROM}_2.fastq.gz"
subset_fastq "${TMP}/normal.names" "$NORMAL_R1" "${DEMO_DIR}/fastq/demo_N_${CHROM}_1.fastq.gz"
subset_fastq "${TMP}/normal.names" "$NORMAL_R2" "${DEMO_DIR}/fastq/demo_N_${CHROM}_2.fastq.gz"

## ---- 3. RNA subset: pairs touching any reference-fusion partner gene (+/- flank) ----
echo "== RNA subset (reference fusion partner genes)"
grep -v '^#' "$REF_FUSIONS" | awk -F'\t' 'NF>=2 {print $1; print $2}' | LC_ALL=C sort -u > "${TMP}/fusion_genes.txt"
awk -F'\t' -v F="$FLANK_BP" '
    NR == FNR { g[$1] = 1; next }
    $3 == "gene" && match($9, /gene_name "[^"]+"/) {
        n = substr($9, RSTART + 11, RLENGTH - 12)
        if (n in g) { s = $4 - F; if (s < 1) s = 1; print $1 "\t" s - 1 "\t" $5 + F "\t" n; found[n] = 1 }
    }
    END { for (x in g) if (!(x in found)) print "WARN: gene not in GTF: " x > "/dev/stderr" }' \
    "${TMP}/fusion_genes.txt" "$CTAT_GTF" | LC_ALL=C sort -k1,1 -k2,2n > "${TMP}/fusion_genes.bed"
## chimeric pairs (supplementary/SA alignment, or mate on another contig): all kept; the others subsampled
samtools view -F 0x100 ${REF_FASTA:+-T "$REF_FASTA"} -L "${TMP}/fusion_genes.bed" "$RNA_BAM" \
    | awk -F'\t' '{ chim = (int($2 / 2048) % 2 == 1 || ($7 != "=" && $7 != "*") || $0 ~ /\tSA:Z:/); print $1 "\t" chim }' \
    | LC_ALL=C sort -u > "${TMP}/rna.flags"
awk -F'\t' '$2 == 1 { print $1 }' "${TMP}/rna.flags" | LC_ALL=C sort -u > "${TMP}/rna.chim"
awk -F'\t' '{ print $1 }' "${TMP}/rna.flags" | LC_ALL=C sort -u | LC_ALL=C comm -23 - "${TMP}/rna.chim" > "${TMP}/rna.other"
shuf -n "$MAX_RNA_PAIRS" --random-source=<(yes) "${TMP}/rna.other" > "${TMP}/rna.other.kept" || true
LC_ALL=C sort -u "${TMP}/rna.chim" "${TMP}/rna.other.kept" > "${TMP}/rna.names"
echo "  RNA: $(awk 'END { print NR }' "${TMP}/rna.chim") chimeric + $(awk 'END { print NR }' "${TMP}/rna.other.kept") of $(awk 'END { print NR }' "${TMP}/rna.other") other pairs"
subset_fastq "${TMP}/rna.names" "$RNA_R1" "${DEMO_DIR}/fastq/demo_R_fusiongenes_1.fastq.gz"
subset_fastq "${TMP}/rna.names" "$RNA_R2" "${DEMO_DIR}/fastq/demo_R_fusiongenes_2.fastq.gz"

## ---- 4. samplesheets (same schema as download_demo_fastq.sh) ----------------
HEADER="patient_id,sample_type,sex,batch,flowcell,run_num,lane_num,chunk,read1,read2"
T_ROW="${PATIENT_ID},DT,F,SEQC2_WES,DEMO_T,1,1,1,demo_T_${CHROM}_1.fastq.gz,demo_T_${CHROM}_2.fastq.gz"
N_ROW="${PATIENT_ID},DN,F,SEQC2_WES,DEMO_N,1,1,1,demo_N_${CHROM}_1.fastq.gz,demo_N_${CHROM}_2.fastq.gz"
R_ROW="${PATIENT_ID},RN,F,SEQC2_RNA,DEMO_R,1,1,1,demo_R_fusiongenes_1.fastq.gz,demo_R_fusiongenes_2.fastq.gz"
printf '%s\n%s\n%s\n' "$HEADER" "$T_ROW" "$N_ROW" > "${DEMO_DIR}/samplesheets/dna_somatic.csv"
printf '%s\n%s\n'     "$HEADER" "$N_ROW"          > "${DEMO_DIR}/samplesheets/dna_germline.csv"
printf '%s\n%s\n'     "$HEADER" "$R_ROW"          > "${DEMO_DIR}/samplesheets/rna_fusion.csv"
cp "${TMP}/fusion_genes.bed" "${DEMO_DIR}/demo_rna_fusion_genes.bed"

## ---- 5. sizes ---------------------------------------------------------------
echo "== Demo size"
{   printf 'file\tread_pairs\tbytes\n'
    for f in "${DEMO_DIR}"/fastq/*_1.fastq.gz; do
        printf '%s\t%s\t%s\n' "$(basename "$f" _1.fastq.gz)" "$(n_reads "$f")" "$(( $(stat -c %s "$f") + $(stat -c %s "${f%_1.fastq.gz}_2.fastq.gz") ))"
    done
} | tee "${DEMO_DIR}/demo_sizes.tsv"
rm -rf "$TMP"
echo "Done. fastq_root = ${DEMO_DIR}/fastq ; demo chromosome = ${CHROM}"
