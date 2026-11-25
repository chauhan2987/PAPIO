#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Date   : 2026-11-12
# Purpose: Pre-flight check for a PAPIO run. Reads the same params YAML you pass
#          to `nextflow run -params-file`, and checks - for the assay/analysis
#          and stages that YAML selects - that software, containers, references,
#          indexes and input FASTQs are present and consistent.
#          Nothing is submitted; nothing is modified.
#
# Usage  : bash check_installation.sh <params.yaml> [path/to/nextflow.config]
# Exit   : 0 = no FAIL (WARNs allowed), 1 = at least one FAIL
#
# Note   : reads flat "key: value" YAML (the format of the PAPIO params files).
#------------------------------------------------------------------------------
set -uo pipefail

###### GLOBAL VARIABLES
PARAMS=${1:?usage: check_installation.sh <params.yaml> [nextflow.config]}
NF_CONFIG=${2:-$(dirname "$0")/../nextflow.config}
NF_TESTED="26.04.1"
N_PASS=0; N_WARN=0; N_FAIL=0

[ -s "$PARAMS" ] || { echo "ERROR: params file not found: $PARAMS" >&2; exit 1; }

#------------------------------------------------------------------------------
##### Functions

## Function to read a top-level scalar from the params YAML (default if absent/null)
yget() {
    local v
    v=$(awk -v k="$1" '
        $0 ~ "^"k"[[:space:]]*:" {
            sub(/^[^:]*:[[:space:]]*/, ""); sub(/[[:space:]]+#.*$/, "")
            gsub(/^["\047]|["\047]$/, ""); print; exit }' "$PARAMS")
    if [ -z "$v" ] || [ "$v" = "null" ] || [ "$v" = "~" ]; then v="${2:-}"; fi
    echo "$v"
}

## Function to read a container filename (YAML override, else nextflow.config default)
cget() {
    local v; v=$(yget "$1")
    if [ -z "$v" ] && [ -s "$NF_CONFIG" ]; then
        v=$(grep -E "^[[:space:]]*$1[[:space:]]*=" "$NF_CONFIG" | head -1 | sed -E "s/.*=[[:space:]]*'([^']+)'.*/\1/")
    fi
    echo "$v"
}

## Functions to report
pass() { echo "  PASS  $*"; N_PASS=$((N_PASS+1)); }
warn() { echo "  WARN  $*"; N_WARN=$((N_WARN+1)); }
fail() { echo "  FAIL  $*"; N_FAIL=$((N_FAIL+1)); }
section() { echo; echo "== $*"; }

## Function to report a missing item as FAIL (level=req) or WARN (level=opt)
miss() { if [ "$1" = "req" ]; then fail "$2"; else warn "$2"; fi; }

## Function to check a file param:  chk_file <param> <req|opt> [why]
chk_file() {
    local key=$1 level=$2 why=${3:-} p; p=$(yget "$key")
    if [ -z "$p" ]; then miss "$level" "$key not set${why:+ ($why)}"; return 1; fi
    if [ -s "$p" ]; then pass "$key"; return 0; fi
    miss "$level" "$key: file missing or empty: $p"; return 1
}

## Function to check a directory param: chk_dir <param> <req|opt> [why]
chk_dir() {
    local key=$1 level=$2 why=${3:-} p; p=$(yget "$key")
    if [ -z "$p" ]; then miss "$level" "$key not set${why:+ ($why)}"; return 1; fi
    if [ -d "$p" ] && [ -n "$(ls -A "$p" 2>/dev/null)" ]; then pass "$key"; return 0; fi
    miss "$level" "$key: directory missing or empty: $p"; return 1
}

## Function to check a bgzipped VCF/BED param and its tabix/csi index
chk_indexed() {
    local key=$1 level=$2 why=${3:-} p
    chk_file "$key" "$level" "$why" || return 1
    p=$(yget "$key")
    if [ -s "${p}.tbi" ] || [ -s "${p}.csi" ]; then pass "$key index"
    else miss "$level" "$key: no .tbi/.csi index next to $p"; fi
}

## Function: does the run touch a stage? (same ranks as utils.stageRankMap)
rank() {
    case "$1" in
        mapping) echo 1;; markdup) echo 2;; metrics) echo 3;; bqsr) echo 4;;
        variant_calling) echo 5;; annotation) echo 6;; classification) echo 7;;
        *) echo 0;;
    esac
}
runs() { local r; r=$(rank "$1"); [ "$START_R" -le "$r" ] && [ "$r" -le "$STOP_R" ]; }

## Function to compare dotted versions: ver_ge A B  -> true if A >= B
ver_ge() { [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -1)" = "$2" ]; }

#------------------------------------------------------------------------------
##### Main Script

ASSAY=$(yget assay wes)
ANALYSIS=$(yget analysis germline)
START_AT=$(yget start_at mapping)
STOP_AT=$(yget stop_at annotation)
START_R=$(rank "$START_AT"); STOP_R=$(rank "$STOP_AT")

echo "PAPIO pre-flight check"
echo "  params   : $PARAMS"
echo "  assay    : $ASSAY   analysis: $ANALYSIS   stages: $START_AT -> $STOP_AT"
[ "$START_R" -gt 0 ] && [ "$STOP_R" -gt 0 ] && [ "$START_R" -le "$STOP_R" ] \
    || { echo "ERROR: invalid start_at/stop_at" >&2; exit 1; }

## ---- Software ---------------------------------------------------------------
section "Software"
if command -v nextflow >/dev/null 2>&1; then
    NF_VER=$(nextflow -version 2>/dev/null | grep -oE 'version [0-9]+\.[0-9]+\.[0-9]+' | head -1 | awk '{print $2}')
    if [ -n "$NF_VER" ] && ver_ge "$NF_VER" "$NF_TESTED"; then pass "nextflow $NF_VER"
    else warn "nextflow ${NF_VER:-unknown}: PAPIO is tested with ${NF_TESTED} (strict syntax)"; fi
else
    fail "nextflow not on PATH"
fi
if command -v singularity >/dev/null 2>&1; then pass "singularity $(singularity --version 2>/dev/null | head -1)"
elif command -v apptainer >/dev/null 2>&1; then pass "apptainer $(apptainer --version 2>/dev/null | head -1)"
else fail "neither singularity nor apptainer on PATH"; fi
BIN_DIR="$(dirname "$NF_CONFIG")/bin"
if [ -d "$BIN_DIR" ]; then
    n_bad=0
    for f in "$BIN_DIR"/*; do [ -x "$f" ] || { warn "bin script not executable: $f (chmod +x)"; n_bad=$((n_bad + 1)); }; done
    [ "$n_bad" -eq 0 ] && pass "bin/ helper scripts present and executable ($(ls "$BIN_DIR" | wc -l))"
else
    fail "bin/ folder missing next to nextflow.config ($BIN_DIR)"
fi

## ---- Run paths --------------------------------------------------------------
section "Run paths"
PR=$(yget pipeline_root)
if [ -n "$PR" ] && [ -s "${PR}/main.nf" ]; then pass "pipeline_root (main.nf found)"
else fail "pipeline_root: no main.nf in ${PR:-not set}"; fi
for key in work_root outdir; do
    p=$(yget "$key")
    if [ -z "$p" ]; then fail "$key not set"; continue; fi
    ## nearest existing ancestor must be writable (Nextflow creates the rest)
    a=$p; while [ ! -e "$a" ] && [ "$a" != "/" ]; do a=$(dirname "$a"); done
    if [ -d "$a" ] && [ -w "$a" ]; then
        if [ "$a" = "$p" ]; then pass "$key (writable)"; else pass "$key (will be created under $a)"; fi
    else
        fail "$key not writable: nearest existing parent $a"
    fi
done
grep -qE '^bind_paths[[:space:]]*:' "$PARAMS" && pass "bind_paths set" \
    || fail "bind_paths not set (required: list of directories containers must see)"

## ---- Containers -------------------------------------------------------------
section "Containers"
CDIR=$(yget container_dir)
case "${ASSAY}:${ANALYSIS}" in
    rna:*)          REQ="container_star container_arriba container_starfusion container_fusioncatcher container_gatk container_annotation container_qc container_align" ;;
    *:somatic)      REQ="container_align container_gatk container_qc container_manta container_strelka container_strelka_af container_muse container_varscan container_msisensor container_annotation container_annotsv container_sequenza container_vcf2maf container_somalier" ;;
    *:germline)     REQ="container_align container_gatk container_qc container_manta container_annotation container_annotsv container_somalier" ;;
    *)              REQ="container_align container_gatk container_qc" ;;
esac
## drop images whose step is switched off in this params file
drop() { REQ=$(echo " $REQ " | sed "s/ $1 / /"); }
[ "$(yget use_fusioncatcher true)" = "false" ] && drop container_fusioncatcher
[ "$(yget use_arriba true)" = "false" ]        && drop container_arriba
[ "$(yget use_starfusion true)" = "false" ] && [ "$(yget use_fusioninspector true)" = "false" ] && drop container_starfusion
[ "$(yget use_muse_in_consensus true)" = "false" ]     && drop container_muse
[ "$(yget use_varscan_in_consensus true)" = "false" ]  && drop container_varscan
[ "$(yget use_strelka2_in_consensus true)" = "false" ] && drop container_strelka
[ "$(yget generate_maf true)" = "false" ]              && drop container_vcf2maf
[ "$ANALYSIS" = "germline" ] && [ "$(yget germline_caller haplotypecaller)" = "deepvariant" ] && REQ="$REQ container_deepvariant"
[ "$(yget mutational_signatures false)" = "true" ] && REQ="$REQ container_signatures"
if [ -z "$CDIR" ] || [ ! -d "$CDIR" ]; then
    fail "container_dir missing: ${CDIR:-not set}"
else
    for c in $REQ; do
        f=$(cget "$c")
        if [ -z "$f" ]; then warn "$c: no filename in params or $NF_CONFIG"
        elif [ -s "${CDIR}/${f}" ]; then pass "$c -> $f"
        else fail "$c: ${CDIR}/${f} not found"; fi
    done
fi

## ---- Core reference ---------------------------------------------------------
section "Reference genome"
chk_file ref_fasta req
chk_file ref_fai   req
if [ "$ASSAY" != "rna" ]; then
    chk_file ref_dict req
    ## the .dict and .fai must list the same contigs in the same order
    FAI=$(yget ref_fai); DICT=$(yget ref_dict)
    if [ -s "$FAI" ] && [ -s "$DICT" ]; then
        if diff -q <(cut -f1 "$FAI") <(grep '^@SQ' "$DICT" | tr '\t' '\n' | grep '^SN:' | cut -c4-) >/dev/null; then
            pass "ref_dict and ref_fai contigs identical"
        else
            fail "ref_dict and ref_fai contig lists differ (dict built from a different FASTA?)"
        fi
    fi
    if runs mapping; then
        IDX=$(yget ref_index_dir)
        if [ -n "$IDX" ] && ls "$IDX"/*.bwt.2bit.64 >/dev/null 2>&1 && ls "$IDX"/*.0123 >/dev/null 2>&1; then
            pass "ref_index_dir has bwa-mem2 index"
        else
            fail "ref_index_dir: no bwa-mem2 index (*.0123, *.bwt.2bit.64) in ${IDX:-not set}"
        fi
    fi
fi

## ---- DNA arm ----------------------------------------------------------------
if [ "$ASSAY" != "rna" ]; then
    section "DNA resources ($ASSAY, $ANALYSIS)"
    if [ "$ASSAY" = "wes" ]; then
        chk_file interval_bed req "WES capture targets"
        chk_file bait_interval_list   opt "CollectHsMetrics skipped"
        chk_file target_interval_list opt "CollectHsMetrics skipped"
        chk_indexed call_regions_gz   opt "Manta/Strelka2 run without call regions"
    fi
    runs metrics && chk_dir verifybamid_db opt "contamination QC skipped"
    runs metrics && chk_file coverage_report_bed opt "no per-gene coverage of reportable genes"
    if runs bqsr && [ "$(yget bqsr false)" = "true" ]; then
        chk_indexed dbsnp req "BQSR"; chk_indexed mills req "BQSR"; chk_indexed known_indels req "BQSR"
    fi
    if runs variant_calling; then
        chk_dir scatter_intervals req "interval scatter"
        chk_indexed somalier_sites opt "sample-swap check skipped"
        if [ -n "$(yget somalier_ancestry_labels)" ]; then chk_file somalier_ancestry_labels opt; chk_dir somalier_1kg_dir opt "ancestry skipped"; fi
        if [ "$ANALYSIS" = "somatic" ]; then
            chk_indexed germline_resource req "Mutect2"
            chk_indexed panel_of_normals  opt "more Mutect2 false positives"
            chk_indexed pileup_variants   opt "contamination filtering skipped"
            chk_file cnvkit_pon     opt "CNVkit uses the matched normal as its reference"
            chk_file msisensor_list opt "MSI skipped"
            chk_file sequenza_gc_wiggle opt "Sequenza purity/LOH skipped"
        fi
        if [ "$ANALYSIS" = "germline" ]; then
            [ "$(yget germline_caller haplotypecaller)" = "haplotypecaller" ] && chk_indexed dbsnp req "HaplotypeCaller"
            for k in hapmap omni thousandg_snps mills; do chk_indexed "$k" opt "germline filtering needs it"; done
            if [ "$(yget run_gcnv true)" = "true" ]; then
                chk_file preprocessed_intervals req "gCNV"; chk_file contig_ploidy_priors req "gCNV"
                chk_dir  gcnv_scatter_intervals req "gCNV"
            fi
            [ "$ASSAY" = "wgs" ] && [ "$(yget run_roh true)" = "true" ] && chk_indexed roh_af_file opt "ROH skipped"
        fi
    fi
    if runs annotation; then
        chk_dir vep_cache_dir req "VEP"
        VC=$(yget vep_cache_dir); VV=$(yget vep_cache_version 115)
        if [ -n "$VC" ] && [ -d "$VC" ]; then
            [ -d "${VC}/homo_sapiens/${VV}_GRCh38" ] && pass "VEP cache release ${VV}_GRCh38 present" \
                || fail "VEP cache: ${VC}/homo_sapiens/${VV}_GRCh38 not found (vep_cache_version=${VV})"
        fi
        chk_indexed clinvar_vcf opt "ClinVar annotation skipped"
        chk_indexed cosmic_vcf  opt "COSMIC recurrence skipped"
        chk_dir annotsv_dir opt "SV/CNV annotation skipped"
        chk_indexed gnomad_faf_vcf opt "population AF = VEP MAX_AF instead of the gnomAD filtering AF"
        VP=$(yget vep_plugins_dir)
        if [ -n "$VP" ] && [ -s "${VP}/NMD.pm" ]; then pass "vep_plugins_dir NMD.pm"
        else miss opt "vep_plugins_dir/NMD.pm not found (NMD escape judged from the last exon only)"; fi
        if [ "$ANALYSIS" = "somatic" ]; then
            chk_file tmb_coding_bed opt "TMB territory = capture BED only"
            chk_file gene_bed opt "no gene-level LOH / MSI tables"
            chk_file arm_bed  opt "no chromosome-arm table"
        fi
    fi
    runs classification && chk_file cosmic_cgc_tsv opt "no CGC gene roles in classification"
fi

## ---- RNA arm ----------------------------------------------------------------
if [ "$ASSAY" = "rna" ]; then
    section "RNA resources"
    chk_file gtf req
    SI=$(yget star_index); CL=$(yget ctat_lib)
    ## CTAT library: needed by STAR-Fusion, FusionInspector, rna_use_ctat_reference, or when star_index is unset
    if [ "$(yget use_starfusion true)" = "true" ] || [ "$(yget use_fusioninspector true)" = "true" ] \
       || [ "$(yget rna_use_ctat_reference true)" = "true" ] || [ -z "$SI" ]; then
        chk_dir ctat_lib req
    fi
    [ -z "$SI" ] && SI="${CL}/ref_genome.fa.star.idx"
    [ -s "${SI}/SA" ] && pass "STAR index: $SI" || fail "STAR index missing (no SA file): $SI"
    ## GTF and FASTA must use the same contig names (STAR-Fusion kickstart yields 0 fusions otherwise)
    G=$(yget gtf); FAI=$(yget ref_fai)
    if [ -s "$G" ] && [ -s "$FAI" ]; then
        C=$( (zcat -f "$G" 2>/dev/null || true) | grep -v '^#' | head -1 | cut -f1)
        cut -f1 "$FAI" | grep -qx "$C" && pass "GTF contig naming matches FASTA ($C)" \
            || fail "GTF contig '$C' not in ref_fai - chromosome naming mismatch"
    fi
    if runs variant_calling; then
        [ "$(yget use_arriba true)" = "true" ] && { chk_file arriba_blacklist opt "Arriba precision drops sharply"
            chk_file arriba_known_fusions opt; chk_file arriba_protein_domains opt; }
        [ "$(yget use_fusioncatcher true)" = "true" ] && chk_dir fusioncatcher_db req "or set use_fusioncatcher: false"
        chk_file known_fusions_tsv opt "known-fusion rescue limited"
    fi
    runs metrics && { chk_file refflat opt "RNA-seq metrics skipped"; chk_file rrna_interval_list opt "%rRNA skipped"; }
    if runs mapping; then
        [ -z "$(yget salmon_index)" ] || fail "salmon_index is no longer used: remove it and set salmon_transcripts (README Appendix A, step 8.7)"
        TX=$(yget salmon_transcripts); GTFX=$(yget gtf)
        [ "$(yget rna_use_ctat_reference true)" = true ] && GTFX="$(yget ctat_lib)/ref_annot.gtf"
        if [ -z "$TX" ]; then miss opt "salmon_transcripts not set (Salmon quantification skipped)"
        elif [ ! -s "$TX" ]; then fail "salmon_transcripts: $TX not found"
        elif [ -s "$GTFX" ]; then
            ## every transcript of the STAR index's GTF must be in the FASTA (else Salmon stops)
            NMISS=$(awk -F'\t' 'NR == FNR { if (sub(/^>/, "")) { split($0, a, /[ \t]/); h[a[1]] = 1 }; next }
                $3 == "exon" && match($9, /transcript_id "[^"]+"/) { t = substr($9, RSTART + 15, RLENGTH - 16); if (!(t in h) && !(t in seen)) { seen[t] = 1; n++ } }
                END { print n + 0 }' "$TX" <(zcat -f "$GTFX"))
            [ "$NMISS" -eq 0 ] && pass "salmon_transcripts matches the STAR index's GTF" \
                || fail "salmon_transcripts lacks ${NMISS} transcripts of ${GTFX}: rebuild it from that GTF (scripts/gtf_to_transcripts.sh)"
        else pass "salmon_transcripts"; fi
        T4=$(yget trust4_ref_dir)
        if [ -n "$T4" ]; then [ -s "${T4}/hg38_bcrtcr.fa" ] && [ -s "${T4}/human_IMGT+C.fa" ] && pass "trust4_ref_dir" \
            || fail "trust4_ref_dir: hg38_bcrtcr.fa / human_IMGT+C.fa missing in $T4"
        else miss opt "trust4_ref_dir not set (TRUST4 skipped)"; fi
        chk_file telescope_gtf opt "Telescope HERV expression skipped"
    fi
    runs annotation && chk_file cosmic_cgc_tsv opt "fusion tiering limited"
fi

## ---- Samplesheet + FASTQs ---------------------------------------------------
section "Samplesheet and inputs"
if [ "$ASSAY" = "rna" ]; then SHEET=$(yget samplesheet "$(yget rna_samplesheet)")
else SHEET=$(yget samplesheet "$(yget dna_samplesheet)"); fi
if [ -z "$SHEET" ] || [ ! -s "$SHEET" ]; then
    fail "samplesheet not found: ${SHEET:-not set}"
else
    pass "samplesheet: $SHEET"
    HDR=$(head -1 "$SHEET" | tr -d '\r')
    REQ_COLS="patient_id sample_type sex batch flowcell run_num lane_num chunk read1"
    [ "$ASSAY" != "rna" ] && REQ_COLS="$REQ_COLS read2"
    for col in $REQ_COLS; do
        echo ",$HDR," | grep -q ",$col," || fail "samplesheet missing column: $col"
    done
    if runs mapping; then
        FQ=$(yget fastq_root)
        [ -d "$FQ" ] || fail "fastq_root missing: ${FQ:-not set}"
        N_MISS=0; N_ROWS=0
        while IFS=, read -r -a F; do
            N_ROWS=$((N_ROWS+1))
            for col in read1 read2; do
                i=$(echo "$HDR" | tr ',' '\n' | grep -nx "$col" | cut -d: -f1)
                [ -z "$i" ] && continue
                v=$(echo "${F[$((i-1))]:-}" | tr -d '\r')
                [ -z "$v" ] && continue
                [ -s "${FQ}/${v}" ] || { N_MISS=$((N_MISS+1)); [ "$N_MISS" -le 5 ] && fail "FASTQ missing: ${FQ}/${v}"; }
            done
        done < <(tail -n +2 "$SHEET")
        [ "$N_MISS" -eq 0 ] && pass "all FASTQs present (${N_ROWS} rows)"
        [ "$N_MISS" -gt 5 ] && fail "... ${N_MISS} FASTQs missing in total"
    else
        IB=$(yget input_bam_dir); IV=$(yget input_vcf_dir)
        if [ "$START_AT" = "annotation" ]; then
            [ -d "$IV" ] && pass "input_vcf_dir" || fail "start_at=annotation needs input_vcf_dir"
        else
            [ -d "$IB" ] && pass "input_bam_dir" || fail "start_at=${START_AT} needs input_bam_dir"
        fi
    fi
fi

## ---- Summary ----------------------------------------------------------------
echo
echo "Summary: ${N_PASS} PASS, ${N_WARN} WARN, ${N_FAIL} FAIL"
if [ "$N_FAIL" -gt 0 ]; then
    echo "Fix the FAIL items before launching."; exit 1
fi
echo "Ready to launch."
