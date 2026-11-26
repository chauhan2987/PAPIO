#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: One-command setup for PAPIO after `git clone`.
#
#   bash setup.sh --test <dir> [--bundle <tar.gz file or URL>]
#       TEST RUN (installation and execution check): pulls the containers and
#       the small test bundle (demo windows of one chromosome, demo reads, mini
#       CTAT library), builds its bwa-mem2 and STAR indexes and writes
#       ready-to-run params files. One server: 16 CPUs, 64 GB RAM.
#
#   bash setup.sh --full <dir> [--kit ID --kit-regions BED --kit-probes BED] [--cosmic]
#                              [--no-fusioncatcher] [--no-annotsv] [--no-gnomad] [--no-salmon]
#       FULL RESOURCES (reproduces the paper): pulls the containers, downloads
#       every public GRCh38 resource (README Appendix A, sections 2-8), builds the
#       indexes and derived files, and writes params files with the paths.
#       --kit ...   capture-kit files from the design BEDs (SureDesign: <ID>_Regions.bed,
#                   probes <ID>_Covered.bed or MergedProbes): interval lists, calling
#                   shards, call regions and gCNV intervals (Appendix A, sections 4 and 5.2)
#       --cosmic    COSMIC v104 (licence): e-mail and password from COSMIC_EMAIL /
#                   COSMIC_PASSWORD, or asked for (Appendix A, section 7.4)
#       About 250 GB of downloads, 300 GB on disk at the end, 64 GB RAM for the
#       bwa-mem2 index; the steps resume where they stopped when re-run.
#
#   Options: --containers-only   only pull the containers
#            --no-containers     skip pulling (the .sif files are already in <dir>/containers)
#            --check-urls        only check that every download URL answers, then exit
#            --skip-url-check    do not check the URLs before downloading
#            --keep-downloads    keep the downloaded archives after unpacking
#            --no-salmon         skip the Salmon transcript FASTA (salmon_transcripts: null; Salmon
#                                quantifies STAR's transcriptome alignments, so it needs no index)
#            --threads N         threads for index building (default: all, at most 16)
#
#   Layout under <dir> (the same as README Appendix A):
#     containers/            *.sif
#     params/                params files with the paths filled in
#     references/downloads/  raw downloads (removed after use unless --keep-downloads)
#     references/genome/     GRCh38.d1.vd1.fa(.fai) .dict, bwa_mem2/
#     references/gatk/       GATK bundle VCFs restricted to the reference contigs (upstream names)
#     references/somalier/ clinvar/ gnomad/ vep_cache/ vep_plugins/ cosmic/ annotsv/
#     references/annotation/ gene, MANE CDS, arm and reportable-gene files, gCNV contig-ploidy priors
#     references/sequenza/ msi/
#     references/capture/<KIT>/  <KIT>_targets.bed, <KIT>_baits.bed, their .interval_list,
#                            <KIT>_targets_pad100.bed.gz, <KIT>_scatter_24_pad100/, <KIT>_scatter_150_pad100/,
#                            gcnv/ (<KIT>_preprocessed_intervals.interval_list, <KIT>_annotated_intervals.tsv,
#                            <KIT>_scatter_24_pad250/); built later: <KIT>_cnvkit_pon.cnn,
#                            <KIT>_inhouse_pon.vcf.gz
#     references/gencode/ star/ (StarIndex_gencode_v49/, gencode.v49.annotation.transcripts.fa) rna_qc/ ctat/ arriba/ fusion/ fusioncatcher/
#     references/trust4/ telescope/
#
#   Site overrides (mirrors, registry): export PAPIO_SETUP_CONF=<file>; the file is
#   sourced after the defaults below and may redefine any of the URL variables.
#
# Needs  : bash, curl, tar, gzip, md5sum, awk, singularity or apptainer;
#          nextflow to run the pipeline afterwards.
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
PAPIO_VERSION=1.0.0
PIPELINE_DIR=$(CDPATH= builtin cd "$(dirname "$0")" >/dev/null && builtin pwd -P)
NF_CONFIG="${PIPELINE_DIR}/nextflow.config"
## container registry: an image is pulled as ${REGISTRY}:<image name without .sif>
REGISTRY=oras://ghcr.io/chauhan2987/papio_containers
## test bundle on Zenodo (filled in after uploading the output of scripts/make_test_bundle.sh)
TEST_BUNDLE_URL=https://zenodo.org/records/XXXXXXX/files/papio_test_bundle_v${PAPIO_VERSION}.tar.gz
TEST_BUNDLE_MD5=XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## images the test run never uses (FusionCatcher off, HaplotypeCaller not DeepVariant, no signatures)
TEST_SKIP_IMAGES="fusioncatcher.sif deepvar_freebayes.sif sigprofiler.sif"

## FULL resources
URL_GENOME=https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834
URL_GATK_HG38=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
URL_GATK_SOM=https://storage.googleapis.com/gatk-best-practices/somatic-hg38
URL_SOMALIER=https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
URL_VEP_CACHE=https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
URL_NMD_PLUGIN=https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/115/NMD.pm
URL_CLINVAR=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
URL_GNOMAD_FAF=https://zenodo.org/records/23531816/files/variantsGnomADv4.1Joint.vcf.bgz
GNOMAD_FAF_MD5=92f6f7c9f7f375ad64b171de6fa7d2b0
URL_ENSEMBL_GTF=https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
URL_CYTOBAND=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
URL_GENCODE_GTF=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
URL_ARRIBA=https://github.com/suhrig/arriba/releases/download/v2.5.1/arriba_v2.5.1.tar.gz
URL_CTAT=https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
URL_FUSIONCATCHER=https://sourceforge.net/projects/fusioncatcher/files/data
URL_ANNOTSV=https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz
URL_SOMALIER_LABELS=https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv
URL_SOMALIER_1KG="https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1"
URL_TRUST4=https://raw.githubusercontent.com/liulab-dfci/TRUST4/master
URL_HERV=https://media.githubusercontent.com/media/mlbendall/telescope_annotation_db/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf
URL_COSMIC_API=https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted
COSMIC_VERSION=v104
## reportable genes for the germline LoF report and per-gene coverage (README Appendix A, step 6.2);
## written once to references/annotation/reportable_genes.txt - edit that file for your own list
REPORTABLE_GENES="VHL PBRM1 SETD2 BAP1 KDM5C SMARCB1 FH FLCN TSC1 TSC2 MTOR PTEN TP53 CDKN2A NF2 MET SDHA SDHB SDHC SDHD TFE3 TFEB"
## GATK bundle VCFs: <bundle> <file name>
GATK_VCFS="hg38 Homo_sapiens_assembly38.dbsnp138.vcf
hg38 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
hg38 Homo_sapiens_assembly38.known_indels.vcf.gz
hg38 hapmap_3.3.hg38.vcf.gz
hg38 1000G_omni2.5.hg38.vcf.gz
hg38 1000G_phase1.snps.high_confidence.hg38.vcf.gz
som  af-only-gnomad.hg38.vcf.gz
som  1000g_pon.hg38.vcf.gz
som  small_exac_common_3.hg38.vcf.gz"
## extra STAR genomeGenerate options (e.g. --genomeSAindexNbases for a small genome)
STAR_GENOME_EXTRA=""
FUSIONCATCHER_PARTS="human_v102.tar.gz.aa human_v102.tar.gz.ab human_v102.tar.gz.ac human_v102.tar.gz.ad"

# shellcheck disable=SC1090
[ -n "${PAPIO_SETUP_CONF:-}" ] && { [ -s "$PAPIO_SETUP_CONF" ] || { echo "ERROR: PAPIO_SETUP_CONF=$PAPIO_SETUP_CONF not found" >&2; exit 1; }; . "$PAPIO_SETUP_CONF"; }

#------------------------------------------------------------------------------
##### Functions

## Function to print usage and exit
usage() { sed -n '3,30p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

## Function to stop with a message
die() { echo "ERROR: $*" >&2; exit 1; }

## Function to print a step header
step() { echo "== $*"; }

## Function to pick the container runtime
runtime() {
    if command -v singularity >/dev/null 2>&1; then echo singularity
    elif command -v apptainer >/dev/null 2>&1; then echo apptainer
    else die "install singularity or apptainer first"; fi
}

## Function to read a container file name from nextflow.config (container_<key>)
sif_of() {
    local f
    f=$(sed -nE "s/^[[:space:]]*container_$1[[:space:]]*=[[:space:]]*'([^']+)'.*/\1/p" "$NF_CONFIG" | head -n 1)
    [ -n "$f" ] || die "container_$1 not found in $NF_CONFIG"
    echo "${CONTAINERS}/${f}"
}

## Function to run a command in a container: sx <container key> <command...>
sx() {
    local sif; sif=$(sif_of "$1"); shift
    [ -s "$sif" ] || die "container $sif missing (run without --no-containers, or pull it)"
    "$RT" exec --bind "$DIR" "$sif" "$@"
}

## Function to download a URL to a file; resumes a partial download, skips a finished one
fetch() {
    local url=$1 out=$2
    if [ -s "$out" ] && [ ! -e "${out}.part" ]; then echo "  present: $(basename "$out")"; return 0; fi
    echo "  downloading: $(basename "$out")"
    mkdir -p "$(dirname "$out")"; touch "${out}.part"
    ## -C - resumes; a server without range support restarts the file instead of failing
    curl -fL -# --retry 5 --retry-delay 10 -C - -o "$out" "$url" \
        || { rm -f "$out"; curl -fL -# --retry 5 --retry-delay 10 -o "$out" "$url"; } \
        || die "download failed: $url"
    rm -f "${out}.part"
}

## Function to check that a URL answers (HEAD; some servers only answer GET)
url_ok() {
    curl -sfIL --retry 2 -o /dev/null "$1" 2>/dev/null \
        || curl -sfL --retry 2 -r 0-0 -o /dev/null "$1" 2>/dev/null
}

## Function to list every URL the full setup downloads
full_urls() {
    local b name
    echo "$URL_GENOME"; echo "$URL_SOMALIER"; echo "$URL_VEP_CACHE"; echo "$URL_NMD_PLUGIN"; echo "$URL_CLINVAR"
    echo "$URL_ENSEMBL_GTF"; echo "$URL_CYTOBAND"; echo "$URL_GENCODE_GTF"; echo "$URL_ARRIBA"; echo "$URL_CTAT"
    while read -r b name; do
        [ "$b" = hg38 ] && echo "${URL_GATK_HG38}/${name}" || echo "${URL_GATK_SOM}/${name}"
    done <<< "$GATK_VCFS"
    echo "$URL_SOMALIER_LABELS"; echo "$URL_SOMALIER_1KG"
    echo "${URL_TRUST4}/hg38_bcrtcr.fa"; echo "${URL_TRUST4}/human_IMGT+C.fa"; echo "$URL_HERV"
    [ "$DO_GNOMAD" -eq 1 ] && { echo "$URL_GNOMAD_FAF"; echo "${URL_GNOMAD_FAF}.tbi"; }
    [ "$DO_ANNOTSV" -eq 1 ] && echo "$URL_ANNOTSV"
    if [ "$DO_FUSIONCATCHER" -eq 1 ]; then for p in $FUSIONCATCHER_PARTS; do echo "${URL_FUSIONCATCHER}/${p}/download"; done; fi
}

## Function to check every URL before anything is downloaded
check_urls() {
    local u bad=0
    step "checking download URLs"
    while read -r u; do
        if url_ok "$u"; then echo "  OK    $u"; else echo "  FAIL  $u"; bad=1; fi
    done < <(full_urls)
    [ "$bad" -eq 0 ] || die "some URLs do not answer (see FAIL above); fix them in setup.sh or PAPIO_SETUP_CONF"
}

## Function to pull the container images named in nextflow.config (test mode skips unused images)
pull_containers() {
    local img n=0
    mkdir -p "$CONTAINERS"
    for img in $(sed -nE "s/^[[:space:]]*container_[a-z0-9_]+[[:space:]]*=[[:space:]]*'([^']+\.sif)'.*/\1/p" "$NF_CONFIG" | sort -u); do
        n=$((n + 1))
        if [ "$MODE" = "test" ] && [[ " ${TEST_SKIP_IMAGES} " == *" ${img} "* ]]; then continue; fi
        if [ -s "${CONTAINERS}/${img}" ]; then echo "  present: ${img}"; continue; fi
        echo "  pulling: ${img}"
        "$RT" pull "${CONTAINERS}/${img}.part" "${REGISTRY}:${img%.sif}" \
            || die "could not pull ${REGISTRY}:${img%.sif}"
        mv "${CONTAINERS}/${img}.part" "${CONTAINERS}/${img}"
    done
    [ "$n" -gt 0 ] || die "no container_* entries found in $NF_CONFIG"
}

## Function to keep only the reference contigs in a VCF (bgzipped, indexed): reconcile <in> <out>
## (the GATK bundle and somalier VCFs carry ALT/HLA contigs that GRCh38.d1.vd1 does not have)
reconcile() {
    local in=$1 out=$2
    [ -s "$out" ] && [ -s "${out}.tbi" ] && { echo "  present: $(basename "$out")"; return 0; }
    echo "  reference contigs only: $(basename "$out")"
    sx align bcftools view -t "$CONTIGS" -Oz -o "${out}.tmp.vcf.gz" "$in"
    mv "${out}.tmp.vcf.gz" "$out"
    sx align tabix -f -p vcf "$out"
}

## Function to merge overlapping intervals of a sorted BED3 on stdin
merge3() {
    awk 'BEGIN { OFS = "\t" } { if ($1 == c && $2 <= e) { if ($3 > e) e = $3 } else { if (c != "") print c, s, e; c = $1; s = $2; e = $3 } }
         END { if (c != "") print c, s, e }'
}

## Function to write a Picard refFlat from a GENCODE GTF (exons in ascending order, also on the minus strand;
## non-coding transcripts get CDS start = CDS end = transcript end)
make_refflat() {
    local gtf=$1 out=$2
    awk -F'\t' 'BEGIN { OFS = "\t" }
        $3 != "exon" && $3 != "CDS" { next }
        { match($9, /transcript_id "[^"]+"/); t = substr($9, RSTART + 15, RLENGTH - 16)
          s = $4 - 1; e = $5 }
        $3 == "exon" { if (!(t in chr)) { ord[++n] = t; chr[t] = $1; str[t] = $7
                         g = t; if (match($9, /gene_name "[^"]+"/)) g = substr($9, RSTART + 11, RLENGTH - 12); gn[t] = g }
                       ne[t]++; es[t, ne[t]] = s; ee[t, ne[t]] = e; next }
        { if (!(t in c0) || s < c0[t]) c0[t] = s; if (!(t in c1) || e > c1[t]) c1[t] = e }
        END { for (i = 1; i <= n; i++) { t = ord[i]; k = ne[t]
                for (a = 1; a <= k; a++) { x[a] = es[t, a]; y[a] = ee[t, a] }
                for (a = 2; a <= k; a++) { vx = x[a]; vy = y[a]; b = a - 1
                  while (b >= 1 && x[b] > vx) { x[b + 1] = x[b]; y[b + 1] = y[b]; b-- }
                  x[b + 1] = vx; y[b + 1] = vy }
                sl = ""; el = ""; for (a = 1; a <= k; a++) { sl = sl x[a] ","; el = el y[a] "," }
                cs = (t in c0) ? c0[t] : y[k]; ce = (t in c1) ? c1[t] : y[k]
                print gn[t], t, chr[t], str[t], x[1], y[k], cs, ce, k, sl, el } }' "$gtf" > "${out}.tmp"
    mv "${out}.tmp" "$out"
}

## Function to write the rRNA interval list (dictionary header + rRNA / Mt_rRNA transcripts)
make_rrna_list() {
    local gtf=$1 dict=$2 out=$3
    { grep '^@' "$dict"
      awk -F'\t' '$3 == "transcript" && $9 ~ /gene_type "(Mt_)?rRNA(_pseudogene)?"/ {
          match($9, /transcript_id "[^"]+"/); print $1 "\t" $4 "\t" $5 "\t" $7 "\t" substr($9, RSTART + 15, RLENGTH - 16) }' "$gtf"
    } > "${out}.tmp"
    mv "${out}.tmp" "$out"
}

## Function to set one key of a params YAML: yset <file> <key> <value>
## (replaces the key's line, or appends it when the file does not have the key)
yset() {
    local f=$1 k=$2 v=$3
    if grep -qE "^${k}[[:space:]]*:" "$f"; then
        awk -v k="$k" -v v="$v" '$0 ~ "^" k "[[:space:]]*:" { sub(/:.*/, ": " v) } { print }' "$f" > "${f}.tmp" || die "cannot update $f"; mv "${f}.tmp" "$f"
    else
        printf '%s : %s\n' "$k" "$v" >> "$f"
    fi
}

## Function to read one key of a params YAML (value without the comment)
yget() { sed -nE "s/^$2[[:space:]]*:[[:space:]]*([^#]*[^[:space:]#]).*/\1/p" "$1" | head -n 1 | tr -d "'\""; }

## Function to write a clean BED3 (no track/browser/comment lines), sorted and merged: clean_bed <in> <out>
clean_bed() {
    awk 'BEGIN { OFS = "\t" } !/^(#|track|browser)/ && NF >= 3 { print $1, $2, $3 }' "$1" \
        | LC_ALL=C sort -k1,1 -k2,2n | merge3 > "${2}.tmp"
    mv "${2}.tmp" "$2"
}

## Function to build the capture-kit files (README Appendix A, sections 4.1-4.4 and 5.2)
make_kit() {
    local C="${R}/capture/${KIT_ID}" g d i n bad
    mkdir -p "$C"
    clean_bed "$KIT_REGIONS" "${C}/${KIT_ID}_targets.bed"
    clean_bed "$KIT_PROBES"  "${C}/${KIT_ID}_baits.bed"
    ## contig names must be those of the reference (chr1, not 1)
    bad=$(cut -f1 "${C}/${KIT_ID}_targets.bed" "${C}/${KIT_ID}_baits.bed" | sort -u \
          | awk -F'\t' 'NR == FNR { ok[$1] = 1; next } !($1 in ok)' <(sed -n 's/^@SQ.*SN:\([^[:space:]]*\).*/\1/p' "$DICT") - | head -n 3 || true)
    [ -z "$bad" ] || die "kit BED contigs not in the reference: $(echo $bad) (UCSC chr names expected)"
    [ -s "${C}/${KIT_ID}_baits.interval_list" ] || sx gatk gatk BedToIntervalList -I "${C}/${KIT_ID}_baits.bed" -SD "$DICT" -O "${C}/${KIT_ID}_baits.interval_list"
    [ -s "${C}/${KIT_ID}_targets.interval_list" ] || sx gatk gatk BedToIntervalList -I "${C}/${KIT_ID}_targets.bed" -SD "$DICT" -O "${C}/${KIT_ID}_targets.interval_list"
    ## padded calling shards: targets +-100 bp, overlaps merged, no target split between shards
    for n in 24 150; do
        d="${C}/${KIT_ID}_scatter_${n}_pad100"
        [ -n "$(ls "$d" 2>/dev/null)" ] && continue
        rm -rf "${d}.tmp"
        sx gatk gatk SplitIntervals -R "$REF" -L "${C}/${KIT_ID}_targets.interval_list" --interval-padding 100 \
            --interval-merging-rule ALL --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION --scatter-count "$n" -O "${d}.tmp"
        mv "${d}.tmp" "$d"
    done
    ## call regions for Strelka2 / Manta, built from the somatic shards
    if [ ! -s "${C}/${KIT_ID}_targets_pad100.bed.gz.tbi" ]; then
        grep -hv '^@' "${C}/${KIT_ID}_scatter_24_pad100"/*.interval_list | awk 'BEGIN { OFS = "\t" } { print $1, $2 - 1, $3 }' \
            | LC_ALL=C sort -k1,1 -k2,2n | merge3 > "${C}/${KIT_ID}_targets_pad100.bed"
        sx align bgzip -f "${C}/${KIT_ID}_targets_pad100.bed"
        sx align tabix -f -p bed "${C}/${KIT_ID}_targets_pad100.bed.gz"
    fi
    ## germline copy-number (gCNV) intervals; run_gcnv itself stays off (needs >= 30 normals)
    g="${C}/gcnv"; mkdir -p "$g"
    [ -s "${g}/${KIT_ID}_preprocessed_intervals.interval_list" ] || sx gatk gatk PreprocessIntervals -R "$REF" -L "${C}/${KIT_ID}_targets.interval_list" \
        --bin-length 0 --padding 250 --interval-merging-rule OVERLAPPING_ONLY -O "${g}/${KIT_ID}_preprocessed_intervals.interval_list"
    [ -s "${g}/${KIT_ID}_annotated_intervals.tsv" ] || sx gatk gatk AnnotateIntervals -R "$REF" -L "${g}/${KIT_ID}_preprocessed_intervals.interval_list" \
        --interval-merging-rule OVERLAPPING_ONLY -O "${g}/${KIT_ID}_annotated_intervals.tsv"
    if [ -z "$(ls "${g}/${KIT_ID}_scatter_24_pad250" 2>/dev/null)" ]; then
        rm -rf "${DL}/${KIT_ID}_gcnv_scatter_tmp" "${g}/${KIT_ID}_scatter_24_pad250"; mkdir -p "${g}/${KIT_ID}_scatter_24_pad250"
        sx gatk gatk IntervalListTools --INPUT "${g}/${KIT_ID}_preprocessed_intervals.interval_list" \
            --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_COUNT 24 --OUTPUT "${DL}/${KIT_ID}_gcnv_scatter_tmp"
        i=0; for d in "${DL}/${KIT_ID}_gcnv_scatter_tmp"/temp_*; do i=$((i + 1)); cp "${d}/scattered.interval_list" "$(printf '%s/%s_scatter_24_pad250/shard_%04d.interval_list' "$g" "$KIT_ID" "$i")"; done
        rm -rf "${DL}/${KIT_ID}_gcnv_scatter_tmp"
    fi
}

## Function to download one COSMIC file: cosmic_get <path in the COSMIC bucket> <output file>
## (the scripted download returns a temporary link for the authenticated user)
cosmic_get() {
    local link
    [ -s "$2" ] && { echo "  present: $(basename "$2")"; return 0; }
    link=$(curl -sf -H "Authorization: Basic ${COSMIC_AUTH}" "${URL_COSMIC_API}?path=$1&bucket=downloads" \
           | sed -n 's/.*"url" *: *"\([^"]*\)".*/\1/p') || true
    [ -n "$link" ] || die "COSMIC refused the download of $1 (check e-mail, password and licence)"
    fetch "$link" "$2"
}

## Function to fill the params templates for the full resources
write_full_params() {
    local f out k v
    mkdir -p "${DIR}/params"
    for f in "${PIPELINE_DIR}"/params/params_{somatic,germline,rna_fusion}.yaml; do
        out="${DIR}/params/$(basename "$f")"
        ## the templates use the README layout (${projectDir}/../references etc.); write this install's folders
        sed -e "s#\${projectDir}/\.\./references#${R}#g" -e "s#\${projectDir}/\.\./containers#${CONTAINERS}#g" \
            -e "s#\${projectDir}/\.\./runs#${DIR}/runs#g" -e "s#\${projectDir}/\.\./demo#${DIR}/demo#g" \
            -e "s#clinvar_YYYYMMDD#clinvar_${CV_DATE}#g" "$f" > "$out"
        ## optional resources that are not on disk: switched off instead of failing at run time
        if [ -n "$KIT_ID" ]; then
            case "$f" in *rna_fusion*) ;;
                *)  C="${R}/capture/${KIT_ID}"
                    yset "$out" interval_bed         "${C}/${KIT_ID}_targets.bed"
                    yset "$out" bait_interval_list   "${C}/${KIT_ID}_baits.interval_list"
                    yset "$out" target_interval_list "${C}/${KIT_ID}_targets.interval_list"
                    yset "$out" call_regions_gz      "${C}/${KIT_ID}_targets_pad100.bed.gz"
                    yset "$out" scatter_intervals    "${C}/${KIT_ID}_scatter_24_pad100"
                    yset "$out" germline_scatter_intervals "${C}/${KIT_ID}_scatter_150_pad100"
                    case "$f" in *germline*)
                        yset "$out" preprocessed_intervals   "${C}/gcnv/${KIT_ID}_preprocessed_intervals.interval_list"
                        yset "$out" gcnv_annotated_intervals "${C}/gcnv/${KIT_ID}_annotated_intervals.tsv"
                        yset "$out" gcnv_scatter_intervals   "${C}/gcnv/${KIT_ID}_scatter_24_pad250"
                        yset "$out" contig_ploidy_priors     "${R}/annotation/GRCh38_contig_ploidy_priors.tsv" ;;
                    esac ;;
            esac
        fi
        if [ "$DO_COSMIC" -eq 1 ]; then
            yset "$out" cosmic_cgc_tsv "${R}/cosmic/Cosmic_CancerGeneCensus_${COSMIC_VERSION}_GRCh38.tsv.gz"
            case "$f" in *rna_fusion*) ;; *) yset "$out" cosmic_vcf "${R}/cosmic/Cosmic_GenomeScreensMutant_Normal_${COSMIC_VERSION}_GRCh38.chr_dedup.vcf.gz" ;; esac
        fi
        for k in gnomad_faf_vcf annotsv_dir fusioncatcher_db verifybamid_db cnvkit_pon cosmic_vcf cosmic_cgc_tsv \
                 salmon_transcripts trust4_ref_dir telescope_gtf somalier_ancestry_labels somalier_1kg_dir coverage_report_bed germline_report_genes; do
            v=$(yget "$out" "$k")
            if [ -n "$v" ] && [ "$v" != null ] && [ ! -e "$v" ]; then yset "$out" "$k" "null    # not found by setup.sh: $v"; fi
        done
        case "$f" in *rna_fusion*) [ "$(yget "$out" fusioncatcher_db)" != null ] || yset "$out" use_fusioncatcher false ;; esac
    done
}

#------------------------------------------------------------------------------
##### Main Script

MODE=""; DIR=""; PULL=1; ONLY_CONTAINERS=0; CHECK_ONLY=0; CHECK_URLS=1; KEEP=0; BUNDLE_SRC=""
KIT_ID=""; KIT_REGIONS=""; KIT_PROBES=""; DO_COSMIC=0; DO_SALMON=1
DO_FUSIONCATCHER=1; DO_ANNOTSV=1; DO_GNOMAD=1
THREADS=$(nproc 2>/dev/null || echo 4); [ "$THREADS" -gt 16 ] && THREADS=16
while [ $# -gt 0 ]; do
    case "$1" in
        --test|--full)      MODE=${1#--}; [ $# -ge 2 ] || usage; DIR=$2; shift 2 ;;
        --bundle)           [ $# -ge 2 ] || usage; BUNDLE_SRC=$2; shift 2 ;;
        --threads)          [ $# -ge 2 ] || usage; THREADS=$2; shift 2 ;;
        --no-containers)    PULL=0; shift ;;
        --containers-only)  ONLY_CONTAINERS=1; shift ;;
        --check-urls)       CHECK_ONLY=1; shift ;;
        --skip-url-check)   CHECK_URLS=0; shift ;;
        --keep-downloads)   KEEP=1; shift ;;
        --no-fusioncatcher) DO_FUSIONCATCHER=0; shift ;;
        --no-annotsv)       DO_ANNOTSV=0; shift ;;
        --no-gnomad)        DO_GNOMAD=0; shift ;;
        --no-salmon)        DO_SALMON=0; shift ;;
        --cosmic)           DO_COSMIC=1; shift ;;
        --kit)              [ $# -ge 2 ] || usage; KIT_ID=$2; shift 2 ;;
        --kit-regions)      [ $# -ge 2 ] || usage; KIT_REGIONS=$2; shift 2 ;;
        --kit-probes)       [ $# -ge 2 ] || usage; KIT_PROBES=$2; shift 2 ;;
        *) usage ;;
    esac
done
[ -n "$MODE" ] && [ -n "$DIR" ] || usage
if [ -n "$KIT_ID$KIT_REGIONS$KIT_PROBES" ]; then
    [ -n "$KIT_ID" ] && [ -s "$KIT_REGIONS" ] && [ -s "$KIT_PROBES" ] || die "--kit needs --kit-regions and --kit-probes (existing BED files)"
    KIT_REGIONS=$(CDPATH= builtin cd "$(dirname "$KIT_REGIONS")" >/dev/null && builtin pwd -P)/$(basename "$KIT_REGIONS")
    KIT_PROBES=$(CDPATH= builtin cd "$(dirname "$KIT_PROBES")" >/dev/null && builtin pwd -P)/$(basename "$KIT_PROBES")
fi
for t in curl tar gzip md5sum awk; do command -v "$t" >/dev/null 2>&1 || die "'$t' not found on PATH"; done
[ -s "$NF_CONFIG" ] || die "nextflow.config not found next to setup.sh ($PIPELINE_DIR)"
mkdir -p "$DIR"; DIR=$(CDPATH= builtin cd "$DIR" >/dev/null && builtin pwd -P)
CONTAINERS="${DIR}/containers"

echo "== PAPIO ${PAPIO_VERSION} setup (${MODE}) in ${DIR}"
command -v nextflow >/dev/null 2>&1 || echo "  NOTE: nextflow not on PATH (needed to run; https://www.nextflow.io)"
if [ "$MODE" = full ] && [ "$CHECK_ONLY" -eq 1 ]; then check_urls; exit 0; fi
RT=$(runtime)

## ---- containers -------------------------------------------------------------
if [ "$PULL" -eq 1 ]; then step "containers -> ${CONTAINERS}"; pull_containers; fi
[ "$ONLY_CONTAINERS" -eq 1 ] && exit 0

## ---- TEST -------------------------------------------------------------------
if [ "$MODE" = "test" ]; then
    step "test bundle"
    SRC=${BUNDLE_SRC:-$TEST_BUNDLE_URL}
    case "$SRC" in
        *XXXXXXX*) die "the test bundle is not published yet (TEST_BUNDLE_URL in setup.sh is a placeholder).
       Build it with scripts/make_test_bundle.sh and pass it: setup.sh --test $DIR --bundle <papio_test_bundle_*.tar.gz>" ;;
        http://*|https://*|ftp://*|file://*) TB="${DIR}/$(basename "$SRC")"; fetch "$SRC" "$TB" ;;
        *) [ -s "$SRC" ] || die "test bundle not found: $SRC"; TB=$(CDPATH= builtin cd "$(dirname "$SRC")" >/dev/null && builtin pwd -P)/$(basename "$SRC") ;;
    esac
    ## checksum: the published MD5, or the .md5 file written next to a locally built bundle
    if [ -z "$BUNDLE_SRC" ]; then
        [ "$(md5sum "$TB" | cut -d' ' -f1)" = "$TEST_BUNDLE_MD5" ] || die "MD5 mismatch for $TB - delete it and re-run"
    elif [ -s "${TB}.md5" ]; then
        [ "$(md5sum "$TB" | cut -d' ' -f1)" = "$(cut -d' ' -f1 "${TB}.md5")" ] || die "MD5 mismatch for $TB (${TB}.md5)"
    fi
    BNAME=$(tar -tzf "$TB" | head -n 1 | cut -d/ -f1) || true
    [ -n "$BNAME" ] || die "cannot read $TB"
    BUNDLE="${DIR}/${BNAME}"
    [ -d "$BUNDLE" ] || tar -C "$DIR" -xzf "$TB"
    ( cd "$BUNDLE" && md5sum -c --quiet MANIFEST.md5 ) || die "bundle contents corrupted - delete $BUNDLE and re-run"
    ## indexes are not shipped (bundle size): built here once, a few minutes
    step "test indexes (bwa-mem2, STAR)"
    for z in "${BUNDLE}"/reference/*.fa.gz; do
        [ -e "$z" ] || continue
        FA=${z%.gz}; CH=$(basename "$FA" .fa)
        [ -s "$FA" ] || { gzip -dc "$z" > "${FA}.tmp" || die "cannot unpack $z"; mv "${FA}.tmp" "$FA"; }
        mkdir -p "${BUNDLE}/reference/bwa_mem2"
        [ -s "${BUNDLE}/reference/bwa_mem2/${CH}.fa.bwt.2bit.64" ] || sx align bwa-mem2 index -p "${BUNDLE}/reference/bwa_mem2/${CH}.fa" "$FA"
    done
    CB="${BUNDLE}/ctat_mini/ctat_genome_lib_build_dir"
    if [ -s "${CB}/ref_genome.fa" ] && [ ! -s "${CB}/ref_genome.fa.star.idx/SA" ]; then
        ## --genomeSAindexNbases scaled to the mini genome (STAR manual: min(14, log2(length)/2 - 1))
        NB=$(awk '!/^>/ { n += length($0) } END { b = int(log(n) / log(2) / 2 - 1); print (b > 14 ? 14 : (b < 4 ? 4 : b)) }' "${CB}/ref_genome.fa")
        rm -rf "${CB}/ref_genome.fa.star.idx"; mkdir -p "${CB}/ref_genome.fa.star.idx"
        sx star STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "${CB}/ref_genome.fa.star.idx" \
            --genomeFastaFiles "${CB}/ref_genome.fa" --sjdbGTFfile "${CB}/ref_annot.gtf" --sjdbOverhang 100 \
            --genomeSAindexNbases "$NB" --outTmpDir "${CB}/_star_tmp"
        rm -rf "${CB}/_star_tmp"
    fi
    mkdir -p "${DIR}/runs"
    for m in somatic germline rna; do
        T="${BUNDLE}/params_templates/params_test_${m}.yaml"
        [ -s "$T" ] || die "missing $T in the bundle"
        sed -e "s#@BUNDLE@#${BUNDLE}#g" -e "s#@CONTAINERS@#${CONTAINERS}#g" \
            -e "s#@PIPELINE@#$(dirname "$PIPELINE_DIR")#g" -e "s#@RUN@#${DIR}/runs#g" "$T" > "${DIR}/params_test_${m}.yaml"
    done
    [ "$PULL" -eq 1 ] || for s in $(sed -nE "s/^[[:space:]]*container_[a-z0-9_]+[[:space:]]*=[[:space:]]*'([^']+\.sif)'.*/\1/p" "$NF_CONFIG" | sort -u); do
        [[ " ${TEST_SKIP_IMAGES} " == *" ${s} "* ]] || [ -s "${CONTAINERS}/${s}" ] || echo "  WARN: ${CONTAINERS}/${s} missing"
    done
    echo
    echo "Setup complete. Check, then run each mode (add -profile slurm on SLURM; no profile = LSF):"
    echo
    for m in somatic germline rna; do
        echo "  bash ${PIPELINE_DIR}/scripts/check_installation.sh ${DIR}/params_test_${m}.yaml"
        echo "  nextflow run ${PIPELINE_DIR}/main.nf -params-file ${DIR}/params_test_${m}.yaml -profile local"
    done
    EXP="${BUNDLE}/expected"; [ -d "$EXP" ] || EXP="${PIPELINE_DIR}/tests/expected"
    if [ -d "$EXP" ]; then
        echo; echo "Compare with the expected outputs:"
        for m in somatic germline rna; do
            echo "  bash ${PIPELINE_DIR}/scripts/compare_to_expected.sh check ${DIR}/runs/test_${m}/Result ${EXP}/${m}"
        done
    fi
    exit 0
fi

## ---- FULL -------------------------------------------------------------------
[ "$PULL" -eq 1 ] || [ -d "$CONTAINERS" ] || die "--no-containers but ${CONTAINERS} does not exist"
[ "$CHECK_URLS" -eq 1 ] && check_urls
if [ "$DO_COSMIC" -eq 1 ]; then
    CE=${COSMIC_EMAIL:-}; CP=${COSMIC_PASSWORD:-}
    [ -n "$CE" ] || read -r -p "COSMIC e-mail: " CE
    [ -n "$CP" ] || { read -r -s -p "COSMIC password: " CP; echo; }
    COSMIC_AUTH=$(printf '%s:%s' "$CE" "$CP" | base64 | tr -d '\n')
fi
R="${DIR}/references"; DL="${R}/downloads"
mkdir -p "$R"/{genome,gatk,somalier,clinvar,gnomad,annotation,vep_cache,vep_plugins,msi,sequenza,gencode,star,rna_qc,ctat,arriba,fusion} "$DL"
df -Pk "$R" | awk 'NR == 2 && $4 < 300 * 1024 * 1024 { printf "  WARN: %.0f GB free under %s; the full setup needs about 300 GB\n", $4 / 1048576, $6 }'
G="${R}/genome"; REF="${G}/GRCh38.d1.vd1.fa"; DICT="${G}/GRCh38.d1.vd1.dict"

step "genome GRCh38.d1.vd1 (bwa-mem2 index: ~64 GB RAM, 1-2 h)"
if [ ! -s "$REF" ]; then
    fetch "$URL_GENOME" "${DL}/GRCh38.d1.vd1.fa.tar.gz"
    tar -C "$G" -xzf "${DL}/GRCh38.d1.vd1.fa.tar.gz"
    [ -s "$REF" ] || die "GRCh38.d1.vd1.fa not found after unpacking ${DL}/GRCh38.d1.vd1.fa.tar.gz"
    [ "$KEEP" -eq 1 ] || rm -f "${DL}/GRCh38.d1.vd1.fa.tar.gz"
fi
[ -s "${REF}.fai" ] || sx align samtools faidx "$REF"
[ -s "$DICT" ]      || sx gatk gatk CreateSequenceDictionary -R "$REF" -O "$DICT"
mkdir -p "${G}/bwa_mem2"
[ -s "${G}/bwa_mem2/GRCh38.d1.vd1.fa.bwt.2bit.64" ] || sx align bwa-mem2 index -p "${G}/bwa_mem2/GRCh38.d1.vd1.fa" "$REF"
CONTIGS=$(sed -n 's/^@SQ.*SN:\([^[:space:]]*\).*/\1/p' "$DICT" | paste -sd, -)

step "GATK resource bundle and somalier sites (restricted to the reference contigs)"
while read -r B NAME; do
    [ "$B" = hg38 ] && U="${URL_GATK_HG38}/${NAME}" || U="${URL_GATK_SOM}/${NAME}"
    OUT="${R}/gatk/${NAME%.gz}.gz"                     # dbSNP is shipped uncompressed
    [ -s "$OUT" ] && [ -s "${OUT}.tbi" ] && { echo "  present: $(basename "$OUT")"; continue; }
    fetch "$U" "${DL}/${NAME}"
    reconcile "${DL}/${NAME}" "$OUT"
    [ "$KEEP" -eq 1 ] || rm -f "${DL}/${NAME}"
done <<< "$GATK_VCFS"
if [ ! -s "${R}/somalier/sites.hg38.vcf.gz.tbi" ]; then
    fetch "$URL_SOMALIER" "${DL}/sites.hg38.vcf.gz"
    reconcile "${DL}/sites.hg38.vcf.gz" "${R}/somalier/sites.hg38.vcf.gz"
fi

step "somalier ancestry labels and 1000 Genomes profiles"
fetch "$URL_SOMALIER_LABELS" "${R}/somalier/ancestry-labels-1kg.tsv"
if [ ! -d "${R}/somalier/1kg-somalier" ]; then
    fetch "$URL_SOMALIER_1KG" "${DL}/1kg.somalier.tar.gz"
    rm -rf "${DL}/1kg_tmp"; mkdir -p "${DL}/1kg_tmp"; tar -C "${DL}/1kg_tmp" -xzf "${DL}/1kg.somalier.tar.gz"
    d=$(find "${DL}/1kg_tmp" -name '*.somalier' -print -quit); [ -n "$d" ] || die "no .somalier files in 1kg.somalier.tar.gz"
    mv "$(dirname "$d")" "${R}/somalier/1kg-somalier"; rm -rf "${DL}/1kg_tmp"
    [ "$KEEP" -eq 1 ] || rm -f "${DL}/1kg.somalier.tar.gz"
fi

step "ClinVar (chromosome names chr*, release date in the file name)"
fetch "$URL_CLINVAR" "${DL}/clinvar.vcf.gz"
## header read in full: 'head' closing the pipe would make zcat fail under pipefail
CV_DATE=$(zcat "${DL}/clinvar.vcf.gz" | awk -F= '/^##fileDate=/ { gsub(/-/, "", $2); print $2; exit } !/^#/ { exit }' || true)
[ -n "$CV_DATE" ] || die "no ##fileDate in ${DL}/clinvar.vcf.gz"
CV_OUT="${R}/clinvar/clinvar_${CV_DATE}.vcf.gz"
if [ ! -s "${CV_OUT}.tbi" ]; then
    { for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' "$c" "$c"; done; printf 'MT\tchrM\n'; } > "${R}/clinvar/chr_rename.txt"
    sx align bcftools annotate --rename-chrs "${R}/clinvar/chr_rename.txt" -Oz -o "${CV_OUT}.tmp.vcf.gz" "${DL}/clinvar.vcf.gz"
    mv "${CV_OUT}.tmp.vcf.gz" "$CV_OUT"
    sx align tabix -f -p vcf "$CV_OUT"
fi

if [ "$DO_GNOMAD" -eq 1 ]; then
    step "gnomAD v4.1 filtering allele frequency (33.6 GB)"
    GF="${R}/gnomad/variantsGnomADv4.1Joint.vcf.bgz"
    if [ ! -s "${GF}.tbi" ] || [ -e "${GF}.part" ]; then
        fetch "$URL_GNOMAD_FAF" "$GF"
        [ "$(md5sum "$GF" | cut -d' ' -f1)" = "$GNOMAD_FAF_MD5" ] || die "MD5 mismatch for $GF - delete it and re-run"
        fetch "${URL_GNOMAD_FAF}.tbi" "${GF}.tbi"
    fi
    sx align bcftools view -h "$GF" | grep '^##INFO=<ID=fafmax_faf95_max_j,' >/dev/null \
        || die "$GF has no INFO/fafmax_faf95_max_j"
fi

step "VEP 115 cache and NMD plugin"
if [ ! -d "${R}/vep_cache/homo_sapiens/115_GRCh38" ]; then
    fetch "$URL_VEP_CACHE" "${DL}/homo_sapiens_vep_115_GRCh38.tar.gz"
    tar -C "${R}/vep_cache" -xzf "${DL}/homo_sapiens_vep_115_GRCh38.tar.gz"
    [ "$KEEP" -eq 1 ] || rm -f "${DL}/homo_sapiens_vep_115_GRCh38.tar.gz"
fi
fetch "$URL_NMD_PLUGIN" "${R}/vep_plugins/NMD.pm"

step "Ensembl 115 gene BED and MANE CDS BED (TMB territory); chromosome arms"
fetch "$URL_ENSEMBL_GTF" "${DL}/Homo_sapiens.GRCh38.115.gtf.gz"
if [ ! -s "${R}/annotation/ensembl115_genes.bed" ]; then
    zcat "${DL}/Homo_sapiens.GRCh38.115.gtf.gz" | awk -F'\t' '$3 == "gene" && $1 ~ /^([0-9]+|X|Y|MT)$/ {
        id = ""; sym = ""; chr = ($1 == "MT") ? "chrM" : "chr" $1
        if (match($9, /gene_id "[^"]+"/))   id  = substr($9, RSTART + 9, RLENGTH - 10)
        if (match($9, /gene_name "[^"]+"/)) sym = substr($9, RSTART + 11, RLENGTH - 12)
        if (sym == "") sym = id
        printf "%s\t%d\t%d\t%s\t%s\n", chr, $4 - 1, $5, sym, id }' \
        | LC_ALL=C sort -k1,1 -k2,2n > "${R}/annotation/ensembl115_genes.bed.tmp"
    mv "${R}/annotation/ensembl115_genes.bed.tmp" "${R}/annotation/ensembl115_genes.bed"
fi
if [ ! -s "${R}/annotation/ensembl115_MANE_CDS.bed" ]; then
    zcat "${DL}/Homo_sapiens.GRCh38.115.gtf.gz" | awk -F'\t' 'BEGIN { OFS = "\t" } $3 == "CDS" && $9 ~ /tag "MANE_Select"/ {
        c = $1; if (c == "MT") c = "M"; if (c !~ /^([0-9]+|X|Y|M)$/) next
        print "chr" c, $4 - 1, $5 }' | LC_ALL=C sort -k1,1 -k2,2n | merge3 > "${R}/annotation/ensembl115_MANE_CDS.bed.tmp"
    mv "${R}/annotation/ensembl115_MANE_CDS.bed.tmp" "${R}/annotation/ensembl115_MANE_CDS.bed"
fi
[ -s "${R}/annotation/ensembl115_genes.bed" ] && [ -s "${R}/annotation/ensembl115_MANE_CDS.bed" ] \
    || die "gene or MANE CDS BED is empty (unexpected Ensembl GTF layout)"
## reportable genes (written once; edit for your own list) and their MANE coding exons (column 4 = gene)
RG="${R}/annotation/reportable_genes.txt"
[ -s "$RG" ] || printf '%s\n' $REPORTABLE_GENES > "$RG"
if [ ! -s "${R}/annotation/ensembl115_reportable_genes_MANE_CDS.bed" ] || [ "$RG" -nt "${R}/annotation/ensembl115_reportable_genes_MANE_CDS.bed" ]; then
    zcat "${DL}/Homo_sapiens.GRCh38.115.gtf.gz" | awk -F'\t' -v list="$RG" 'BEGIN { OFS = "\t"; while ((getline g < list) > 0) G[toupper(g)] = 1 }
        $3 == "CDS" && $9 ~ /tag "MANE_Select"/ && match($9, /gene_name "[^"]+"/) {
          g = substr($9, RSTART + 11, RLENGTH - 12); c = $1; if (c == "MT") c = "M"
          if (c ~ /^([0-9]+|X|Y|M)$/ && (toupper(g) in G)) print "chr" c, $4 - 1, $5, g }' \
        | LC_ALL=C sort -k1,1 -k2,2n -u > "${R}/annotation/ensembl115_reportable_genes_MANE_CDS.bed.tmp"
    mv "${R}/annotation/ensembl115_reportable_genes_MANE_CDS.bed.tmp" "${R}/annotation/ensembl115_reportable_genes_MANE_CDS.bed"
fi
if [ ! -s "${R}/annotation/GRCh38_arms.bed" ]; then
    fetch "$URL_CYTOBAND" "${DL}/cytoBand.hg38.txt.gz"
    zcat "${DL}/cytoBand.hg38.txt.gz" | awk -F'\t' 'BEGIN { OFS = "\t" }
        $1 ~ /^chr([0-9]+|X|Y)$/ { k = $1 SUBSEP substr($4, 1, 1)
          if (!(k in s) || $2 < s[k]) s[k] = $2
          if (!(k in e) || $3 > e[k]) e[k] = $3 }
        END { for (k in s) { split(k, a, SUBSEP); arm = substr(a[1], 4) a[2]
                if (arm ~ /^(13|14|15|21|22)p$/) continue      # acrocentric p arms: no assembled sequence
                print a[1], s[k], e[k], arm } }' | LC_ALL=C sort -k1,1 -k2,2n > "${R}/annotation/GRCh38_arms.bed"
fi
## gCNV contig-ploidy priors: genome-level (expected copies of each chromosome), the same for every kit
[ -s "${R}/annotation/GRCh38_contig_ploidy_priors.tsv" ] || awk 'BEGIN { OFS = "\t"
    print "CONTIG_NAME", "PLOIDY_PRIOR_0", "PLOIDY_PRIOR_1", "PLOIDY_PRIOR_2", "PLOIDY_PRIOR_3"
    for (i = 1; i <= 22; i++) print "chr" i, "0.01", "0.01", "0.97", "0.01"
    print "chrX", "0.01", "0.49", "0.49", "0.01"; print "chrY", "0.495", "0.495", "0.01", "0.01" }' > "${R}/annotation/GRCh38_contig_ploidy_priors.tsv"

step "Sequenza GC wiggle and MSIsensor-pro site list"
SQ="${R}/sequenza/GRCh38.d1.vd1.gc50Base.wig.gz"
[ -s "$SQ" ] || { sx sequenza sequenza-utils gc_wiggle --fasta "$REF" -w 50 -o "${SQ}.tmp.gz" || die "Sequenza gc_wiggle failed"; mv "${SQ}.tmp.gz" "$SQ"; }
MS="${R}/msi/GRCh38.d1.vd1.msisensor.list"
[ -s "$MS" ] || { sx msisensor msisensor-pro scan -d "$REF" -o "${MS}.tmp" || die "MSIsensor-pro scan failed"; mv "${MS}.tmp" "$MS"; }

if [ -n "$KIT_ID" ]; then
    step "capture kit ${KIT_ID}: interval lists, calling shards, call regions, gCNV intervals"
    make_kit
fi

if [ "$DO_COSMIC" -eq 1 ]; then
    step "COSMIC ${COSMIC_VERSION} (licence)"
    CO="${R}/cosmic"; mkdir -p "$CO"
    CGC="${CO}/Cosmic_CancerGeneCensus_${COSMIC_VERSION}_GRCh38.tsv.gz"
    CV="${CO}/Cosmic_GenomeScreensMutant_Normal_${COSMIC_VERSION}_GRCh38.chr_dedup.vcf.gz"
    if [ ! -s "$CGC" ]; then
        cosmic_get "grch38/cosmic/${COSMIC_VERSION}/Cosmic_CancerGeneCensus_Tsv_${COSMIC_VERSION}_GRCh38.tar" "${DL}/cosmic_cgc.tar"
        tar -C "$CO" -xf "${DL}/cosmic_cgc.tar" "$(basename "$CGC")" || die "$(basename "$CGC") not in the COSMIC CGC archive"
    fi
    if [ ! -s "${CV}.tbi" ]; then
        cosmic_get "grch38/cosmic/${COSMIC_VERSION}/Cosmic_GenomeScreensMutant_Normal_Vcf_${COSMIC_VERSION}_GRCh38.tar" "${DL}/cosmic_gs.tar"
        tar -C "$DL" -xf "${DL}/cosmic_gs.tar" "Cosmic_GenomeScreensMutant_Normal_${COSMIC_VERSION}_GRCh38.vcf.gz" \
            || die "genome-screens VCF not in the COSMIC archive"
        ## chr names, sorted, transcript-level duplicate rows collapsed to one record per variant
        sx align bash -c "set -euo pipefail
            bcftools annotate --rename-chrs '${R}/clinvar/chr_rename.txt' '${DL}/Cosmic_GenomeScreensMutant_Normal_${COSMIC_VERSION}_GRCh38.vcf.gz' -Ou \
              | bcftools sort -T '${DL}' -Ou | bcftools norm -d exact -Oz -o '${CV}.tmp.vcf.gz'"
        mv "${CV}.tmp.vcf.gz" "$CV"; sx align tabix -f -p vcf "$CV"
        [ "$KEEP" -eq 1 ] || rm -f "${DL}/cosmic_gs.tar" "${DL}/Cosmic_GenomeScreensMutant_Normal_${COSMIC_VERSION}_GRCh38.vcf.gz"
    fi
fi

if [ "$DO_ANNOTSV" -eq 1 ]; then
    step "AnnotSV 3.5 annotations"
    if [ ! -d "${R}/annotsv/Annotations_Human" ]; then
        fetch "$URL_ANNOTSV" "${DL}/Annotations_Human_3.5.tar.gz"
        mkdir -p "${R}/annotsv"; tar -C "${R}/annotsv" -xzf "${DL}/Annotations_Human_3.5.tar.gz"
        [ -d "${R}/annotsv/Annotations_Human" ] || die "Annotations_Human/ not found after unpacking AnnotSV"
        [ "$KEEP" -eq 1 ] || rm -f "${DL}/Annotations_Human_3.5.tar.gz"
    fi
fi

step "RNA: GENCODE v49 STAR index (~40 GB RAM), refFlat, rRNA list; CTAT v44; Arriba database"
GTF49="${R}/gencode/gencode.v49.annotation.gtf"
if [ ! -s "$GTF49" ]; then
    fetch "$URL_GENCODE_GTF" "${DL}/gencode.v49.annotation.gtf.gz"
    zcat "${DL}/gencode.v49.annotation.gtf.gz" > "${GTF49}.tmp" || die "cannot unpack the GENCODE GTF"; mv "${GTF49}.tmp" "$GTF49"
fi
SI="${R}/star/StarIndex_gencode_v49"
if [ ! -s "${SI}/SA" ] || [ ! -s "${SI}/genomeParameters.txt" ]; then
    rm -rf "$SI" "${DL}/star_index_tmp"; mkdir -p "$SI"      # STAR needs a new --outTmpDir
    sx star STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$SI" --genomeFastaFiles "$REF" \
        --sjdbGTFfile "$GTF49" --sjdbOverhang 100 --limitGenomeGenerateRAM 60000000000 --outTmpDir "${DL}/star_index_tmp" \
        ${STAR_GENOME_EXTRA}
    rm -rf "${DL}/star_index_tmp"
fi
## Salmon transcript FASTA: same GTF + genome as the STAR index (Salmon quantifies STAR's
## transcriptome alignments; every transcript of that BAM must be in the FASTA)
TXFA="${R}/star/gencode.v49.annotation.transcripts.fa"
if [ "$DO_SALMON" -eq 1 ] && [ ! -s "$TXFA" ]; then
    cp "${PIPELINE_DIR}/scripts/gtf_to_transcripts.sh" "${DL}/gtf_to_transcripts.sh"     # inside the bound folder
    TMPDIR="$DL" sx align bash "${DL}/gtf_to_transcripts.sh" --gtf "$GTF49" --fasta "$REF" --out "$TXFA" \
        || die "transcript FASTA for Salmon not built"
    rm -f "${DL}/gtf_to_transcripts.sh"
fi
[ -s "${R}/rna_qc/gencode.v49.refFlat" ]           || make_refflat "$GTF49" "${R}/rna_qc/gencode.v49.refFlat"
[ -s "${R}/rna_qc/gencode.v49.rrna.interval_list" ] || make_rrna_list "$GTF49" "$DICT" "${R}/rna_qc/gencode.v49.rrna.interval_list"
CT="${R}/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play"
if [ ! -s "${CT}/ctat_genome_lib_build_dir/ref_annot.gtf" ]; then
    fetch "$URL_CTAT" "${DL}/$(basename "$URL_CTAT")"
    tar -C "${R}/ctat" -xzf "${DL}/$(basename "$URL_CTAT")"
    [ -s "${CT}/ctat_genome_lib_build_dir/ref_annot.gtf" ] || die "CTAT library layout not as expected under ${CT}"
    [ "$KEEP" -eq 1 ] || rm -f "${DL}/$(basename "$URL_CTAT")"
fi
if [ ! -s "${R}/fusion/known_fusions_hg38_GRCh38_v2.5.1.tsv" ]; then
    fetch "$URL_ARRIBA" "${DL}/arriba_v2.5.1.tar.gz"
    tar -C "$DL" -xzf "${DL}/arriba_v2.5.1.tar.gz" --wildcards '*/database/*hg38*'
    find "$DL" -path '*/database/*hg38*' -type f -exec mv {} "${R}/arriba/" \;
    for f in blacklist_hg38_GRCh38_v2.5.1.tsv.gz known_fusions_hg38_GRCh38_v2.5.1.tsv.gz protein_domains_hg38_GRCh38_v2.5.1.gff3; do
        [ -s "${R}/arriba/${f}" ] || die "Arriba database file ${f} missing"
    done
    ## known_fusions_tsv for the fusion annotation step: uncompressed copy of Arriba's list
    zcat "${R}/arriba/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz" > "${R}/fusion/known_fusions_hg38_GRCh38_v2.5.1.tsv"
    [ "$KEEP" -eq 1 ] || rm -rf "${DL}/arriba_v2.5.1.tar.gz" "${DL}"/arriba_v2.5.1
fi
step "RNA immune resources: TRUST4, HERV (Telescope)"
mkdir -p "${R}/trust4" "${R}/telescope"
fetch "${URL_TRUST4}/hg38_bcrtcr.fa"    "${R}/trust4/hg38_bcrtcr.fa"
fetch "${URL_TRUST4}/human_IMGT+C.fa"   "${R}/trust4/human_IMGT+C.fa"
fetch "$URL_HERV"                       "${R}/telescope/HERV_rmsk.hg38.v2.gtf"
grep -q '^>' "${R}/trust4/hg38_bcrtcr.fa" || die "${R}/trust4/hg38_bcrtcr.fa is not a FASTA"
awk -F'\t' '!/^#/ { exit !($3 != "") }' "${R}/telescope/HERV_rmsk.hg38.v2.gtf" || die "HERV GTF looks wrong (Git LFS pointer?)"

if [ "$DO_FUSIONCATCHER" -eq 1 ] && [ ! -d "${R}/fusioncatcher/human_v102" ]; then
    step "FusionCatcher human_v102 (~30 GB)"
    mkdir -p "${R}/fusioncatcher"
    for p in $FUSIONCATCHER_PARTS; do fetch "${URL_FUSIONCATCHER}/${p}/download" "${DL}/${p}"; done
    ## shell glob order (aa, ab, ...) is the part order
    cat "${DL}"/human_v102.tar.gz.a? | tar -C "${R}/fusioncatcher" -xz
    [ -d "${R}/fusioncatcher/human_v102" ] || die "human_v102/ not found after unpacking FusionCatcher"
    [ "$KEEP" -eq 1 ] || rm -f "${DL}"/human_v102.tar.gz.a?
fi

step "params files -> ${DIR}/params"
write_full_params
if [ -s "${PIPELINE_DIR}/scripts/make_reference_manifest.sh" ]; then
    bash "${PIPELINE_DIR}/scripts/make_reference_manifest.sh" "${DIR}/params/params_somatic.yaml" "${DIR}/references_manifest.tsv" --no-md5 >/dev/null \
        || echo "  WARN: reference manifest not written"
fi

cat <<EOF

Public references are in ${R}; params files in ${DIR}/params.
Still needed (cannot be downloaded automatically); the params files point at these paths:
  1. $( [ -n "$KIT_ID" ] && echo "capture kit ${KIT_ID}: built under ${R}/capture/${KIT_ID} and set in the DNA params files" \
        || echo "capture-kit files: re-run with --kit <ID> --kit-regions <Regions.bed> --kit-probes <Covered.bed>
     (design BEDs from Agilent SureDesign; SEQC2 demonstration kit S07604624)")
  2. demonstration FASTQs: bash ${PIPELINE_DIR}/scripts/download_demo_fastq.sh ${DIR}/demo
     (or set dna_samplesheet / rna_samplesheet and fastq_root to your own data)
  3. optional, set to null until present: CNVkit panel of normals for the kit, VerifyBamID2 resource,
     in-house Mutect2 panel of normals (main_pon.nf, from your normals; README Appendix A, step 4.5)$( [ "$DO_COSMIC" -eq 1 ] || echo ", COSMIC (--cosmic)")$( [ "$DO_ANNOTSV" -eq 1 ] || echo ", AnnotSV annotations")$( [ "$DO_GNOMAD" -eq 1 ] || echo ", gnomAD filtering AF")$( [ "$DO_SALMON" -eq 1 ] || echo ", Salmon transcript FASTA (--no-salmon)")
  4. reportable genes: ${R}/annotation/reportable_genes.txt holds the README example list (RCC);
     edit it and re-run setup to rebuild the per-gene coverage BED
Then check and run:
  bash ${PIPELINE_DIR}/scripts/check_installation.sh ${DIR}/params/params_somatic.yaml
EOF
