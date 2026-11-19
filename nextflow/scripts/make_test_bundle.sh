#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Date   : 2026-11-12
# Purpose: Build the PAPIO TEST BUNDLE that reviewers/users download to check
#          that every mode runs end to end on a normal Linux server, without
#          the full references. Run ONCE by the developer; upload the tarball
#          to Zenodo; setup.sh --test downloads it.
#
#   The bundle must stay small (tens of MB) for Zenodo or a GitHub release:
#   indexes are NOT shipped (setup.sh --test builds the bwa-mem2 and STAR
#   indexes), and every resource is cut to the demo windows of make_demo_subset.sh.
#
#   reference/  demo chromosome: gzipped FASTA, .fai, .dict
#   resources/  GATK/gnomAD/ClinVar/somalier VCFs cut to the demo windows,
#               header contigs rewritten to match the cut reference
#   targets/    the demo targets (demo_targets.bed): BED, padded call regions,
#               interval_lists, 4 scatter shards; gene BED of the chromosome
#   vep_cache/  VEP 115 cache: the 1 Mb chunks and known variants of the windows
#   msi/        MSIsensor-pro sites inside the windows
#   ctat_mini/  CTAT genome library built from the reference-fusion partner
#               genes (+/- flank) only, without its STAR index and the Pfam/Dfam
#               databases used to build it
#   fastq/ samplesheets/      demo reads from make_demo_subset.sh
#   params_templates/         filled in by setup.sh (@BUNDLE@ etc.)
#
#   NOT in the bundle (licence or size): COSMIC, AnnotSV annotations,
#   FusionCatcher DB, Arriba blacklist (coordinates would not match the mini
#   RNA genome). The test params switch the dependent steps off.
#
# Usage  : edit GLOBAL VARIABLES, then  bash make_test_bundle.sh
# Needs  : singularity, the PAPIO containers, samtools/bcftools/tabix/bgzip
#          on PATH or via SIF_ALIGN, internet for CTAT Pfam/Dfam download.
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
BUNDLE_VERSION=1.0.0
BUNDLE_DIR=/path/to/papio_test_bundle_v${BUNDLE_VERSION}
DEMO_DIR=/path/to/demo_subset                  # output of make_demo_subset.sh
CPUS=8

## containers (same images the pipeline uses)
CONTAINER_DIR=/path/to/papio_containers
SIF_ALIGN=${CONTAINER_DIR}/bwa_samtools_fastp_fgbio_bcftools.sif   # samtools, bcftools, bwa-mem2
SIF_GATK=${CONTAINER_DIR}/gatk.sif
SIF_MSI=${CONTAINER_DIR}/msisensorpro.sif
SIF_STARFUSION=${CONTAINER_DIR}/starfusion.sif
PREP_GENOME_LIB=/usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl   # path inside SIF_STARFUSION

## full resources (your existing files)
REF_FASTA=/path/to/papio_resources/genome/GRCh38.d1.vd1.fa
ENSEMBL_GTF=/path/to/papio_resources/annotation/Homo_sapiens.GRCh38.115.gtf.gz
VEP_CACHE_DIR=/path/to/papio_resources/vep_cache            # contains homo_sapiens/115_GRCh38/
VEP_CACHE_VERSION=115
CTAT_LIB=/path/to/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir
CTAT_FUSION_ANNOT_LIB=/path/to/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.source/fusion_annot_lib.gz   # leave empty if absent
CTAT_ANNOT_FILTER_RULE=/path/to/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.source/AnnotFilterRule.pm    # leave empty if absent
REF_FUSIONS=/path/to/benchmark/reference_fusions.tsv        # GENE5<TAB>GENE3

## VCF resources to cut:  <bundle file name>  <full path>
VCFS="dbsnp.vcf.gz                /path/to/papio_resources/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz
mills.vcf.gz                      /path/to/papio_resources/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_indels.vcf.gz               /path/to/papio_resources/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz
hapmap.vcf.gz                     /path/to/papio_resources/gatk/hapmap_3.3.hg38.vcf.gz
omni.vcf.gz                       /path/to/papio_resources/gatk/1000G_omni2.5.hg38.vcf.gz
thousandg_snps.vcf.gz             /path/to/papio_resources/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz
germline_resource.vcf.gz          /path/to/papio_resources/gatk/af-only-gnomad.hg38.vcf.gz
panel_of_normals.vcf.gz           /path/to/papio_resources/gatk/1000g_pon.hg38.vcf.gz
pileup_variants.vcf.gz            /path/to/papio_resources/gatk/small_exac_common_3.hg38.vcf.gz
somalier_sites.vcf.gz             /path/to/papio_resources/somalier/sites.hg38.vcf.gz
clinvar.vcf.gz                    /path/to/papio_resources/clinvar/clinvar_YYYYMMDD.vcf.gz"

EXON_PAD=100        # call-region padding (bp), same as the full run
SCATTER_N=4
RNA_FLANK_BP=10000  # must match FLANK_BP used in make_demo_subset.sh
MAX_FILE_MB=95      # GitHub refuses files over 100 MB
MAX_BUNDLE_MB=150   # whole tarball

#------------------------------------------------------------------------------
##### Functions

## Function to run a tool inside a container (plain `singularity exec`)
sx() { local sif=$1; shift; singularity exec "$sif" "$@"; }

#------------------------------------------------------------------------------
##### Main Script

CHROM=$(awk -F'\t' '$1=="# chosen" {print $2}' "${DEMO_DIR}/demo_selection.tsv")
[ -n "$CHROM" ] || { echo "ERROR: no '# chosen' line in ${DEMO_DIR}/demo_selection.tsv" >&2; exit 1; }
echo "Test bundle v${BUNDLE_VERSION} for ${CHROM}"
mkdir -p "${BUNDLE_DIR}"/{reference,resources,targets/scatter,vep_cache,msi,ctat_mini,fastq,samplesheets,params_templates}
B=$BUNDLE_DIR
TMP="${B}/.tmp"; mkdir -p "$TMP"

## ---- 1. reference: demo chromosome only -------------------------------------
echo "== reference"
sx "$SIF_ALIGN" samtools faidx "$REF_FASTA" "$CHROM" > "${B}/reference/${CHROM}.fa"
sx "$SIF_ALIGN" samtools faidx "${B}/reference/${CHROM}.fa"
sx "$SIF_ALIGN" samtools dict "${B}/reference/${CHROM}.fa" -o "${B}/reference/${CHROM}.dict"
WINDOWS="${DEMO_DIR}/demo_windows.bed"
[ -s "$WINDOWS" ] && [ -s "${DEMO_DIR}/demo_targets.bed" ] || { echo "ERROR: run make_demo_subset.sh first (demo_windows.bed)" >&2; exit 1; }

## ---- 2. VCF resources cut to the demo windows, header contigs = cut reference ----
echo "== resources"
while read -r NAME SRC; do
    [ -z "$NAME" ] && continue
    ## reheader --fai keeps only the cut reference's contig(s); it cannot read stdin
    sx "$SIF_ALIGN" bcftools view -R "$WINDOWS" -Oz -o "${TMP}/cut.vcf.gz" "$SRC"
    ## no record in the windows: keep the chromosome's first 500 (an empty --custom VCF stops VEP)
    if [ "$(sx "$SIF_ALIGN" bcftools view -H "${TMP}/cut.vcf.gz" | awk 'END { print NR }')" -eq 0 ]; then
        { sx "$SIF_ALIGN" bcftools view -h "$SRC"
          sx "$SIF_ALIGN" bcftools view -H -r "$CHROM" "$SRC" | awk 'NR <= 500'; } | sx "$SIF_ALIGN" bgzip -c > "${TMP}/cut.vcf.gz"
    fi
    sx "$SIF_ALIGN" bcftools reheader --fai "${B}/reference/${CHROM}.fa.fai" -o "${B}/resources/${NAME}" "${TMP}/cut.vcf.gz"
    sx "$SIF_ALIGN" tabix -f -p vcf "${B}/resources/${NAME}"
    n=$(sx "$SIF_ALIGN" bcftools view -H "${B}/resources/${NAME}" | awk 'END { print NR }')
    echo "  ${NAME}: ${n} records"
    [ "$n" -gt 0 ] || [ "$NAME" != clinvar.vcf.gz ] || { echo "ERROR: ClinVar has no record on ${CHROM} (VEP --custom fails on an empty VCF)" >&2; exit 1; }
done <<< "$VCFS"

## ---- 3. targets: the demo targets (capture targets holding the demo variants) ----
echo "== targets"
NUM=${CHROM#chr}
cut -f1-3 "${DEMO_DIR}/demo_targets.bed" | LC_ALL=C sort -k1,1 -k2,2n \
| awk 'BEGIN{OFS="\t"} NR==1 {c=$1; s=$2; e=$3; next}
       $2 <= e { if ($3 > e) e = $3; next }
       { print c, s, e; c=$1; s=$2; e=$3 } END { print c, s, e }' > "${B}/targets/targets.bed"
awk -v p="$EXON_PAD" 'BEGIN{OFS="\t"} {s=$2-p; if (s<0) s=0; print $1, s, $3+p}' "${B}/targets/targets.bed" \
| awk 'BEGIN{OFS="\t"} NR==1 {c=$1; s=$2; e=$3; next}
       $2 <= e { if ($3 > e) e = $3; next }
       { print c, s, e; c=$1; s=$2; e=$3 } END { print c, s, e }' > "${B}/targets/targets_pad${EXON_PAD}.bed"
sx "$SIF_ALIGN" bgzip -c "${B}/targets/targets_pad${EXON_PAD}.bed" > "${B}/targets/targets_pad${EXON_PAD}.bed.gz"
sx "$SIF_ALIGN" tabix -f -p bed "${B}/targets/targets_pad${EXON_PAD}.bed.gz"
sx "$SIF_GATK" gatk BedToIntervalList -I "${B}/targets/targets.bed" \
    -SD "${B}/reference/${CHROM}.dict" -O "${B}/targets/targets.interval_list"
sx "$SIF_GATK" gatk BedToIntervalList -I "${B}/targets/targets_pad${EXON_PAD}.bed" \
    -SD "${B}/reference/${CHROM}.dict" -O "${B}/targets/targets_pad${EXON_PAD}.interval_list"
sx "$SIF_GATK" gatk SplitIntervals -R "${B}/reference/${CHROM}.fa" \
    -L "${B}/targets/targets_pad${EXON_PAD}.interval_list" \
    --scatter-count "$SCATTER_N" -O "${B}/targets/scatter"
echo "  $(wc -l < "${B}/targets/targets.bed") target intervals, $(awk '{s+=$3-$2} END{print s}' "${B}/targets/targets.bed") bp"

## gene BED for the chromosome (gene LOH, MSI gene map): chr, start, end, symbol, Ensembl ID
## (same layout as README Appendix A, step 6.2; genes without a name keep their ID as symbol)
zcat -f "$ENSEMBL_GTF" \
| awk -F'\t' -v n="$NUM" -v c="$CHROM" 'BEGIN{OFS="\t"}
    $1 == n && $3 == "gene" && match($9, /gene_id "[^"]+"/) {
        id = substr($9, RSTART + 9, RLENGTH - 10); sym = id
        if (match($9, /gene_name "[^"]+"/)) sym = substr($9, RSTART + 11, RLENGTH - 12)
        print c, $4 - 1, $5, sym, id }' \
| LC_ALL=C sort -k1,1 -k2,2n > "${B}/targets/genes.bed"

## ---- 4. VEP cache: this chromosome only -------------------------------------
echo "== VEP cache"
VSRC="${VEP_CACHE_DIR}/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38"
VDST="${B}/vep_cache/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38"
mkdir -p "$VDST"
find "$VSRC" -maxdepth 1 -type f -exec cp {} "$VDST/" \;          # info.txt, chr_synonyms.txt, ...
if   [ -d "${VSRC}/${NUM}" ];   then VC=$NUM
elif [ -d "${VSRC}/${CHROM}" ]; then VC=$CHROM
else echo "ERROR: no ${NUM}/ or ${CHROM}/ folder in ${VSRC}" >&2; exit 1; fi
mkdir -p "${VDST}/${VC}"
## transcript / regulation chunks (<start>-<end>[_reg].gz) that overlap a window
for f in "${VSRC}/${VC}"/*-*.gz; do
    [ -e "$f" ] || continue
    r=$(basename "$f" .gz); r=${r%%_*}; cs=${r%-*}; ce=${r#*-}
    awk -v cs="$cs" -v ce="$ce" '$2 < ce + 0 && $3 >= cs - 1 { f = 1 } END { exit !f }' "$WINDOWS" && cp "$f" "${VDST}/${VC}/"
done
## known variants (indexed cache): only the windows; the window names use the cache's chromosome name.
## The cut file is indexed like the original: index type and sequence / start / end columns read
## from its .csi or .tbi header (no assumption about the cache's column layout).
AV="${VSRC}/${VC}/all_vars.gz"
if [ -s "$AV" ]; then
    if   [ -e "${AV}.csi" ]; then IX=csi; read -r CS CB CE < <(gzip -dc "${AV}.csi" | od -An -t d4 -j 20 -N 12)
    elif [ -e "${AV}.tbi" ]; then IX=tbi; read -r CS CB CE < <(gzip -dc "${AV}.tbi" | od -An -t d4 -j 12 -N 12)
    else echo "ERROR: ${AV} has no .csi / .tbi index" >&2; exit 1; fi
    awk -v c="$VC" 'BEGIN { OFS = "\t" } { print c, $2 + 1, $3 }' "$WINDOWS" > "${TMP}/vep_regions.txt"
    ## windows are sorted and merged, so tabix -R output stays sorted; a record spanning two windows is printed twice
    sx "$SIF_ALIGN" tabix -R "${TMP}/vep_regions.txt" "$AV" | awk '!seen[$0]++' \
        | sx "$SIF_ALIGN" bgzip -c > "${VDST}/${VC}/all_vars.gz"
    if [ "$IX" = csi ]; then sx "$SIF_ALIGN" tabix -C -s "$CS" -b "$CB" -e "$CE" "${VDST}/${VC}/all_vars.gz"
    else sx "$SIF_ALIGN" tabix -s "$CS" -b "$CB" -e "$CE" "${VDST}/${VC}/all_vars.gz"; fi
    ## the cut file must answer a query exactly like the original
    q=$(awk -v c="$VC" 'NR == 1 { print c ":" $2 + 1 "-" $3; exit }' "$WINDOWS")
    [ "$(sx "$SIF_ALIGN" tabix "$AV" "$q" | md5sum)" = "$(sx "$SIF_ALIGN" tabix "${VDST}/${VC}/all_vars.gz" "$q" | md5sum)" ] \
        || { echo "ERROR: cut VEP all_vars.gz does not answer like the original for $q" >&2; exit 1; }
fi
echo "  $(du -sh "${B}/vep_cache" | cut -f1)"

## ---- 5. MSI site list -------------------------------------------------------
echo "== MSI sites"
sx "$SIF_MSI" msisensor-pro scan -d "${B}/reference/${CHROM}.fa" -o "${TMP}/msi_all.list"
## header + sites inside the windows (few windows: a plain scan is fast)
awk -F'\t' 'NR == FNR { n++; ws[n] = $2; we[n] = $3; next }
    FNR == 1 { print; next }
    { for (i = 1; i <= n; i++) if ($2 > ws[i] && $2 <= we[i]) { print; break } }' "$WINDOWS" "${TMP}/msi_all.list" > "${B}/msi/${CHROM}.msisensor.list"
echo "  $(( $(wc -l < "${B}/msi/${CHROM}.msisensor.list") - 1 )) microsatellites in the windows"

## ---- 6. mini CTAT library from the reference-fusion partner genes -----------
echo "== mini CTAT library"
CTAT_GTF="${CTAT_LIB}/ref_annot.gtf"
CTAT_FA="${CTAT_LIB}/ref_genome.fa"
grep -v '^#' "$REF_FUSIONS" | awk -F'\t' 'NF>=2 {print $1; print $2}' | LC_ALL=C sort -u > "${TMP}/genes.txt"

## partner-gene spans +/- flank, merged per chromosome -> mini contigs
awk -F'\t' -v F="$RNA_FLANK_BP" '
    NR == FNR { g[$1] = 1; next }
    $3 == "gene" && match($9, /gene_name "[^"]+"/) {
        n = substr($9, RSTART + 11, RLENGTH - 12)
        if (n in g) { s = $4 - F; if (s < 1) s = 1; print $1 "\t" s "\t" $5 + F }
    }' "${TMP}/genes.txt" "$CTAT_GTF" \
| LC_ALL=C sort -k1,1 -k2,2n \
| awk 'BEGIN{OFS="\t"} NR==1 {c=$1; s=$2; e=$3; next}
       $1 == c && $2 <= e { if ($3 > e) e = $3; next }
       { print c, s, e; c=$1; s=$2; e=$3 } END { print c, s, e }' > "${TMP}/mini_regions.tsv"     # 1-based closed

## FASTA: one contig per region, named <chrom>_<start>_<end> (no ':' for tool safety)
: > "${TMP}/mini_genome.fa"
while IFS=$'\t' read -r C S E; do
    sx "$SIF_ALIGN" samtools faidx "$CTAT_FA" "${C}:${S}-${E}" \
    | sed "1s/.*/>${C}_${S}_${E}/" >> "${TMP}/mini_genome.fa"
done < "${TMP}/mini_regions.tsv"

## GTF: keep genes lying fully inside a region, shift coordinates to the contig
awk -F'\t' 'BEGIN{OFS="\t"}
    NR == FNR { n++; rc[n]=$1; rs[n]=$2; re[n]=$3; next }
    /^#/ { next }
    {
        match($9, /gene_id "[^"]+"/); gid = substr($9, RSTART, RLENGTH)
        if ($3 == "gene") {
            for (i = 1; i <= n; i++)
                if ($1 == rc[i] && $4 >= rs[i] && $5 <= re[i]) { keep[gid] = i; break }
        }
        if (gid in keep) {
            i = keep[gid]
            $1 = rc[i] "_" rs[i] "_" re[i]; $4 = $4 - rs[i] + 1; $5 = $5 - rs[i] + 1
            print
        }
    }' "${TMP}/mini_regions.tsv" "$CTAT_GTF" > "${TMP}/mini_annot.gtf"
## (CTAT GTF lists the gene line before its transcripts/exons, so one pass suffices)
N_GENES=$(awk -F'\t' '$3 == "gene"' "${TMP}/mini_annot.gtf" | wc -l)
[ "$N_GENES" -gt 0 ] || { echo "ERROR: no genes kept for the mini CTAT GTF (does ${CTAT_GTF} have 'gene' lines?)" >&2; exit 1; }

FA_ARGS=""
if [ -n "$CTAT_FUSION_ANNOT_LIB" ] && [ -s "$CTAT_FUSION_ANNOT_LIB" ]; then
    ## only the entries of the mini library's genes (the full library is large)
    zcat -f "$CTAT_FUSION_ANNOT_LIB" | awk -F'\t' 'NR == FNR { g[$1] = 1; next }
        { split($1, p, "--"); if ((p[1] in g) || (p[2] in g)) print }' "${TMP}/genes.txt" - | gzip -c > "${TMP}/fusion_annot_lib.dat.gz"
    FA_ARGS="$FA_ARGS --fusion_annot_lib ${TMP}/fusion_annot_lib.dat.gz"
fi
[ -n "$CTAT_ANNOT_FILTER_RULE" ] && [ -s "$CTAT_ANNOT_FILTER_RULE" ] && FA_ARGS="$FA_ARGS --annot_filter_rule $CTAT_ANNOT_FILTER_RULE"
# shellcheck disable=SC2086
sx "$SIF_STARFUSION" "$PREP_GENOME_LIB" \
    --genome_fa "${TMP}/mini_genome.fa" \
    --gtf "${TMP}/mini_annot.gtf" \
    $FA_ARGS \
    --pfam_db current \
    --dfam_db human \
    --CPU "$CPUS" \
    --output_dir "${B}/ctat_mini/ctat_genome_lib_build_dir"
echo "  $(wc -l < "${TMP}/mini_regions.tsv") contigs, $(grep -c $'\tgene\t' "${TMP}/mini_annot.gtf") genes"
## not shipped: the STAR index (setup.sh --test rebuilds it) and the Pfam / Dfam databases and
## checkpoints of the library build (only the build uses them)
CB="${B}/ctat_mini/ctat_genome_lib_build_dir"
rm -rf "${CB}/ref_genome.fa.star.idx" "${CB}"/__* "${CB}"/*.ok
find "$CB" -maxdepth 1 \( -name 'Pfam-A.hmm*' -o -name '*dfam*.hmm*' -o -name '*Dfam*.hmm*' \) -exec rm -f {} +

## ---- 7. demo reads + samplesheets -------------------------------------------
echo "== demo reads"
cp "${DEMO_DIR}"/fastq/*.fastq.gz "${B}/fastq/"
cp "${DEMO_DIR}"/samplesheets/*.csv "${B}/samplesheets/"
cp "${DEMO_DIR}/demo_selection.tsv" "${DEMO_DIR}/demo_sizes.tsv" "${B}/"

## ---- 8. params templates (setup.sh replaces @BUNDLE@ @CONTAINERS@ @PIPELINE@ @RUN@) ----
echo "== params templates"
PAD=$EXON_PAD
cat > "${B}/params_templates/params_test_common.yaml" <<EOF
pipeline_root   : @PIPELINE@
container_dir   : @CONTAINERS@
bind_paths      : ['@BUNDLE@', '@RUN@']
fastq_root      : @BUNDLE@/fastq
dedup_mode         : markdup
umi_read_structure : null
bqsr               : false
EOF
cat > "${B}/params_templates/params_test_dna_refs.yaml" <<EOF
assay                : wes
ref_index_dir        : @BUNDLE@/reference/bwa_mem2
ref_fasta            : @BUNDLE@/reference/${CHROM}.fa
ref_fai              : @BUNDLE@/reference/${CHROM}.fa.fai
ref_dict             : @BUNDLE@/reference/${CHROM}.dict
interval_bed         : @BUNDLE@/targets/targets.bed
bait_interval_list   : @BUNDLE@/targets/targets_pad${PAD}.interval_list
target_interval_list : @BUNDLE@/targets/targets.interval_list
call_regions_gz      : @BUNDLE@/targets/targets_pad${PAD}.bed.gz
scatter_intervals    : @BUNDLE@/targets/scatter
gene_bed             : @BUNDLE@/targets/genes.bed
dbsnp                : @BUNDLE@/resources/dbsnp.vcf.gz
mills                : @BUNDLE@/resources/mills.vcf.gz
known_indels         : @BUNDLE@/resources/known_indels.vcf.gz
somalier_sites       : @BUNDLE@/resources/somalier_sites.vcf.gz
vep_cache_dir        : @BUNDLE@/vep_cache
vep_cache_version    : ${VEP_CACHE_VERSION}
clinvar_vcf          : @BUNDLE@/resources/clinvar.vcf.gz
annotsv_dir          : null      # AnnotSV annotations not bundled: SV/CNV annotation skipped
cosmic_vcf           : null      # COSMIC is licensed: not bundled
cosmic_cgc_tsv       : null
EOF
{
    cat "${B}/params_templates/params_test_common.yaml" "${B}/params_templates/params_test_dna_refs.yaml"
    cat <<EOF
work_root       : @RUN@/test_somatic
outdir          : @RUN@/test_somatic/Result
dna_samplesheet : @BUNDLE@/samplesheets/dna_somatic.csv
analysis        : somatic
start_at        : mapping
stop_at         : classification
germline_resource : @BUNDLE@/resources/germline_resource.vcf.gz
panel_of_normals  : @BUNDLE@/resources/panel_of_normals.vcf.gz
pileup_variants   : @BUNDLE@/resources/pileup_variants.vcf.gz
msisensor_list    : @BUNDLE@/msi/${CHROM}.msisensor.list
cnvkit_pon         : null    # one chromosome: CNV / purity / TMB / signatures are not meaningful
sequenza_gc_wiggle : null
compute_tmb        : false
mutational_signatures : false
EOF
} > "${B}/params_templates/params_test_somatic.yaml"
{
    cat "${B}/params_templates/params_test_common.yaml" "${B}/params_templates/params_test_dna_refs.yaml"
    cat <<EOF
work_root       : @RUN@/test_germline
outdir          : @RUN@/test_germline/Result
dna_samplesheet : @BUNDLE@/samplesheets/dna_germline.csv
analysis        : germline
start_at        : mapping
stop_at         : classification
germline_caller : haplotypecaller
germline_filter : nv
run_gcnv        : false
hapmap          : @BUNDLE@/resources/hapmap.vcf.gz
omni            : @BUNDLE@/resources/omni.vcf.gz
thousandg_snps  : @BUNDLE@/resources/thousandg_snps.vcf.gz
EOF
} > "${B}/params_templates/params_test_germline.yaml"
{
    cat "${B}/params_templates/params_test_common.yaml"
    cat <<EOF
work_root       : @RUN@/test_rna
outdir          : @RUN@/test_rna/Result
rna_samplesheet : @BUNDLE@/samplesheets/rna_fusion.csv
assay           : rna
analysis        : rna_fusion
start_at        : mapping
stop_at         : classification
ctat_lib        : @BUNDLE@/ctat_mini/ctat_genome_lib_build_dir
ref_fasta       : @BUNDLE@/ctat_mini/ctat_genome_lib_build_dir/ref_genome.fa
ref_fai         : @BUNDLE@/ctat_mini/ctat_genome_lib_build_dir/ref_genome.fa.fai
gtf             : @BUNDLE@/ctat_mini/ctat_genome_lib_build_dir/ref_annot.gtf
star_index      : null
use_arriba          : true     # runs without blacklist: its GRCh38 coordinates do not fit the mini genome
use_starfusion      : true
use_fusioncatcher   : false    # FusionCatcher DB (~19 GB) not bundled
use_fusioninspector : true
arriba_blacklist       : null
arriba_known_fusions   : null
arriba_protein_domains : null
fusioncatcher_db       : null
known_fusions_tsv      : null
fusion_min_callers     : 2
rna_shared_star        : false
EOF
} > "${B}/params_templates/params_test_rna.yaml"
rm -f "${B}/params_templates/params_test_common.yaml" "${B}/params_templates/params_test_dna_refs.yaml"

## ---- 9. manifest + tarball for Zenodo ---------------------------------------
echo "== package"
## reference shipped gzipped (setup.sh --test unpacks it and builds the bwa-mem2 index)
gzip -9 -f "${B}/reference/${CHROM}.fa"
## size limits: no file over MAX_FILE_MB
big=$(find "$B" -type f -size +"${MAX_FILE_MB}"M -printf '%s\t%p\n' | sort -rn)
[ -z "$big" ] || { echo "ERROR: files over ${MAX_FILE_MB} MB (GitHub limit 100 MB):" >&2; echo "$big" >&2; exit 1; }
echo "  largest files:"; find "$B" -type f -printf '%s\t%P\n' | sort -rn | head -n 8 | awk -F'\t' '{ printf "    %8.1f MB  %s\n", $1 / 1048576, $2 }'
rm -rf "$TMP"
( cd "$B" && find . -type f ! -name MANIFEST.md5 -print0 | LC_ALL=C sort -z | xargs -0 md5sum > MANIFEST.md5 )
printf 'bundle_version\t%s\nchromosome\t%s\nbuilt\t%s\n' "$BUNDLE_VERSION" "$CHROM" "$(date +%F)" > "${B}/BUNDLE_INFO.tsv"
TARBALL="$(dirname "$B")/$(basename "$B").tar.gz"
tar -C "$(dirname "$B")" -czf "$TARBALL" "$(basename "$B")"
md5sum "$TARBALL" | tee "${TARBALL}.md5"
echo "Test it: bash setup.sh --test <dir> --bundle ${TARBALL}"
echo "Then upload it to Zenodo and put its URL and MD5 into setup.sh (TEST_BUNDLE_URL, TEST_BUNDLE_MD5)."
du -sh "$TARBALL"
MB=$(( $(stat -c %s "$TARBALL") / 1048576 ))
[ "$MB" -le "$MAX_BUNDLE_MB" ] || { echo "ERROR: bundle is ${MB} MB (limit ${MAX_BUNDLE_MB}); lower N_MSI / MAX_RNA_PAIRS or MIN_SNV in make_demo_subset.sh" >&2; exit 1; }
