#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Pad a BED, clamp it to the contig lengths of a FASTA index and merge overlaps.
#          Used by the PAPIO process DEEPVARIANT (DeepVariant has no padding option).
# Usage  : pad_merge_bed.sh <in.bed[.gz]> <ref.fa.fai> <padding_bp> <out.bed>
#          Regions on contigs absent from the .fai are dropped; an empty result is an error.
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
IN=$1; FAI=$2; PAD=$3; OUT=$4

##### Main Script
case "$IN" in *.gz) cat="zcat" ;; *) cat="cat" ;; esac
$cat "$IN" \
    | awk -F'\t' -v p="$PAD" -v fai="$FAI" 'BEGIN { OFS = "\t"; while ((getline l < fai) > 0) { split(l, a, "\t"); len[a[1]] = a[2] } }
        /^(#|track|browser)/ { next }
        ($1 in len) { s = $2 - p; if (s < 0) s = 0; e = $3 + p; if (e > len[$1]) e = len[$1]; print $1, s, e }' \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | awk -F'\t' 'BEGIN { OFS = "\t" }
        { if ($1 == c && $2 <= e) { if ($3 > e) e = $3 } else { if (c != "") print c, s, e; c = $1; s = $2; e = $3 } }
        END { if (c != "") print c, s, e }' > "$OUT"
[ -s "$OUT" ] || { echo "ERROR: padded regions BED is empty (contig names of BED vs FASTA?)" >&2; exit 1; }
