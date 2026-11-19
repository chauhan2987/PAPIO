#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Date   : 2026-11-12
# Purpose: Write a reference manifest (TSV) for every resource path set in a
#          PAPIO params YAML: path, size, MD5 (files) or a content fingerprint
#          (directories), plus the version/source columns to fill in once.
#          Re-run later and `diff` the two manifests to prove references are
#          unchanged between runs.
#
# Usage  : bash make_reference_manifest.sh <params.yaml> [out.tsv] [--no-md5]
#          --no-md5 skips checksums (fast; sizes only).
#------------------------------------------------------------------------------
set -uo pipefail

###### GLOBAL VARIABLES
PARAMS=${1:?usage: make_reference_manifest.sh <params.yaml> [out.tsv] [--no-md5]}
OUT=${2:-references_manifest.tsv}
DO_MD5=1
for a in "$@"; do [ "$a" = "--no-md5" ] && DO_MD5=0; done
[ "$OUT" = "--no-md5" ] && OUT=references_manifest.tsv

## every params key that points at a reference resource (order = manifest order)
KEYS="ref_fasta ref_fai ref_dict ref_index_dir
interval_bed bait_interval_list target_interval_list call_regions_gz scatter_intervals germline_scatter_intervals
dbsnp mills known_indels germline_resource panel_of_normals pileup_variants
hapmap omni thousandg_snps roh_af_file somalier_sites verifybamid_db
cnvkit_pon msisensor_list sequenza_gc_wiggle gene_bed
preprocessed_intervals contig_ploidy_priors gcnv_scatter_intervals gcnv_annotated_intervals
vep_cache_dir vep_plugins_dir clinvar_vcf cosmic_vcf cosmic_cgc_tsv annotsv_dir
gtf ctat_lib star_index refflat rrna_interval_list
arriba_blacklist arriba_known_fusions arriba_protein_domains fusioncatcher_db known_fusions_tsv"

#------------------------------------------------------------------------------
##### Functions

## Function to read a top-level scalar from the params YAML
yget() {
    awk -v k="$1" '
        $0 ~ "^"k"[[:space:]]*:" {
            sub(/^[^:]*:[[:space:]]*/, ""); sub(/[[:space:]]+#.*$/, "")
            gsub(/^["\047]|["\047]$/, ""); if ($0 != "null" && $0 != "~") print; exit }' "$PARAMS"
}

## Function to pull a version hint from the file itself where one is recorded
version_hint() {
    local p=$1 key=$2
    case "$key" in
        vep_cache_dir) find "$p/homo_sapiens" -mindepth 1 -maxdepth 1 -type d -name '*_GRCh38' -printf '%f\n' 2>/dev/null | sort | paste -sd, - ;;
        ## CTAT: .../GRCh38_gencode_v44_CTAT_lib_<date>.plug-n-play/ctat_genome_lib_build_dir
        ctat_lib) basename "$(cd "$p/.." 2>/dev/null && pwd)" ;;
        ## FusionCatcher: .../human_v102 ; AnnotSV: last two path components
        fusioncatcher_db|annotsv_dir) echo "$(basename "$(dirname "$p")")/$(basename "$p")" ;;
        *) case "$p" in
               *.vcf.gz|*.vcf)
                   (zcat -f "$p" 2>/dev/null | head -200 || true) \
                   | grep -m1 -E '^##(fileDate|source|reference)=' | cut -c3- ;;
               *.gtf|*.gtf.gz)
                   (zcat -f "$p" 2>/dev/null | head -5 || true) \
                   | grep -m1 -iE '^#!?(description|genome-version|provider)' | cut -c2- ;;
           esac ;;
    esac
}

#------------------------------------------------------------------------------
##### Main Script

printf 'param\tpath\ttype\tsize_bytes\tmd5_or_fingerprint\tversion_hint\tversion\tsource_url\tdownload_date\n' > "$OUT"

for key in $KEYS; do
    p=$(yget "$key")
    [ -z "$p" ] && continue
    if [ -f "$p" ]; then
        size=$(stat -c %s "$p")
        if [ "$DO_MD5" -eq 1 ]; then sum=$(md5sum "$p" | cut -d' ' -f1); else sum="NA"; fi
        printf '%s\t%s\tfile\t%s\t%s\t%s\t\t\t\n' "$key" "$p" "$size" "$sum" "$(version_hint "$p" "$key")" >> "$OUT"
    elif [ -d "$p" ]; then
        ## directory fingerprint = MD5 of the sorted (relative path, size) listing:
        ## changes if any file is added, removed or resized; cheap for 30+ GB dirs
        size=$(du -sb "$p" 2>/dev/null | cut -f1)
        sum=$(cd "$p" && find -L . -type f -printf '%P\t%s\n' 2>/dev/null | LC_ALL=C sort | md5sum | cut -d' ' -f1)
        printf '%s\t%s\tdir\t%s\tlisting:%s\t%s\t\t\t\n' "$key" "$p" "$size" "$sum" "$(version_hint "$p" "$key")" >> "$OUT"
    else
        printf '%s\t%s\tMISSING\tNA\tNA\t\t\t\t\n' "$key" "$p" >> "$OUT"
    fi
done

echo "Wrote $OUT ($(($(wc -l < "$OUT") - 1)) resources)."
echo "Fill the version / source_url / download_date columns once; they are not auto-detectable for every resource."
