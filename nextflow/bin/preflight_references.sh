#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Check that resource VCFs and BEDs match the reference (contig names and
#          primary-chromosome lengths) before any job is submitted.
#          Called by the PAPIO process PREFLIGHT_REFERENCES.
# Usage  : preflight_references.sh <ref.dict> <vcf_dir> <bed_dir>
# Output : preflight_references.txt (one OK / WARN / FAIL line per file); exit 1 on any FAIL
#------------------------------------------------------------------------------
set -uo pipefail

###### GLOBAL VARIABLES
DICT=$1; VCF_DIR=$2; BED_DIR=$3
OUT="preflight_references.txt"
fail=0

##### Functions
## Function to record a result line
note() { echo "$1" | tee -a "$OUT" >&2; }

## Function to test whether a contig is in the reference (reads the whole table: pipefail-safe)
in_dict() { awk -F'\t' -v c="$1" '$1 == c { f = 1 } END { exit !f }' dict.tsv; }

## Function to check one VCF: primary-chromosome lengths and the first record's contig
check_vcf() {
    local v=$1 name hdr bad first
    name=$(basename "$v")
    hdr=$(bcftools view -h "$v" 2>/dev/null) || { note "FAIL  $name: not readable by bcftools"; fail=1; return; }
    bad=$(printf '%s\n' "$hdr" | sed -n 's/^##contig=<ID=\([^,>]*\),length=\([0-9]*\).*/\1\t\2/p' \
          | awk -F'\t' 'NR == FNR { L[$1] = $2; next }
              $1 ~ /^(chr)?([0-9]+|X|Y)$/ && ($1 in L) && L[$1] != $2 { printf "%s(%s vs ref %s) ", $1, $2, L[$1] }' dict.tsv -)
    first=$(set +o pipefail; bcftools view -H "$v" 2>/dev/null | head -1 | cut -f1)
    if [ -n "$bad" ]; then
        note "FAIL  $name: contig lengths differ from the reference (different build?): $bad"; fail=1
    elif [ -n "$first" ] && ! in_dict "$first"; then
        case "$name" in
            *clinvar*|*ClinVar*|*cosmic*|*COSMIC*) note "WARN  $name: contig '$first' not in reference (OK only if VEP --synonyms maps it)" ;;
            *) note "FAIL  $name: contig '$first' not in the reference .dict (chr naming / build mismatch)"; fail=1 ;;
        esac
    else
        note "OK    $name"
    fi
}

## Function to check one BED: the first region's contig
check_bed() {
    local b=$1 name cat first
    name=$(basename "$b")
    case "$name" in *.gz) cat="zcat" ;; *) cat="cat" ;; esac
    first=$(set +o pipefail; $cat "$b" | grep -v -E '^(#|track|browser)' | head -1 | cut -f1)
    if [ -z "$first" ]; then note "FAIL  $name: empty"; fail=1
    elif ! in_dict "$first"; then note "FAIL  $name: contig '$first' not in the reference .dict"; fail=1
    else note "OK    $name"; fi
}

##### Main Script
awk -F'\t' '$1 == "@SQ" { n = ""; l = ""
    for (i = 2; i <= NF; i++) { if ($i ~ /^SN:/) n = substr($i, 4); if ($i ~ /^LN:/) l = substr($i, 4) }
    print n "\t" l }' "$DICT" > dict.tsv
: > "$OUT"
for v in "$VCF_DIR"/*; do [ -e "$v" ] && check_vcf "$v"; done
for b in "$BED_DIR"/*; do [ -e "$b" ] && check_bed "$b"; done
[ "$fail" -eq 0 ] || { echo "PREFLIGHT_REFERENCES failed - see messages above" >&2; exit 1; }
exit 0
