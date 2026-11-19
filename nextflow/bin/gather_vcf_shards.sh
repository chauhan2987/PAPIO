#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Concatenate per-shard VCFs (plain text) into one bgzipped, indexed VCF in
#          reference (.fai) contig order. Used by the PAPIO process VARSCAN_MERGE.
# Usage  : gather_vcf_shards.sh <out.vcf.gz> <ref.fa.fai> <shard.vcf> [...]
#          The header is taken from the first shard that has one; shards without calls
#          may be header-less.
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
OUT=$1; FAI=$2; shift 2

##### Main Script
hdr=""
for f in "$@"; do if grep -q '^#CHROM' "$f"; then hdr="$f"; break; fi; done
{
    if [ -n "$hdr" ]; then
        grep '^#' "$hdr"
    else
        printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    fi
    ## body of every shard, ordered by reference contig rank then position
    cat "$@" | { grep -v '^#' || true; } \
        | awk -F'\t' 'BEGIN { OFS = "\t" } NR == FNR { rank[$1] = NR; next }
                       { print (($1 in rank) ? rank[$1] : 999999), $0 }' "$FAI" - \
        | LC_ALL=C sort -t$'\t' -k1,1n -k3,3n | cut -f2-
} | bgzip > "$OUT"
bcftools index -t "$OUT"
