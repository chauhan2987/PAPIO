#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Transcript (cDNA) FASTA from a GTF and its genome, for Salmon
#          alignment mode on STAR's Aligned.toTranscriptome.out.bam. Built from
#          the SAME GTF + genome as the STAR index, so every transcript in the
#          BAM header is in the FASTA (Salmon stops on a missing one).
#          Exons of each transcript are joined in transcript order; minus-strand
#          transcripts are reverse-complemented. Names = GTF transcript_id.
# Usage  : gtf_to_transcripts.sh --gtf <annotation.gtf[.gz]> --fasta <genome.fa> \
#              --out <transcripts.fa>
# Needs  : samtools (faidx), awk, sort; genome.fa.fai next to the FASTA
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
GTF=""
FASTA=""
OUT=""
TMP=""

##### Functions
## Function to print usage and exit
usage() { sed -n '/^# Usage/,/^# Needs/p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 2; }

## Function to parse the options (a missing value is an error)
parse_args() {
    while [ $# -gt 0 ]; do
        case "$1" in
            --gtf|--fasta|--out)
                if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "ERROR: $1 needs a value" >&2; exit 2; fi
                case "$1" in --gtf) GTF=$2 ;; --fasta) FASTA=$2 ;; --out) OUT=$2 ;; esac
                shift 2 ;;
            -h|--help) usage ;;
            *) echo "ERROR: unknown option $1" >&2; usage ;;
        esac
    done
    [ -n "$GTF" ] && [ -n "$FASTA" ] && [ -n "$OUT" ] || usage
    [ -s "$GTF" ] || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
    [ -s "$FASTA" ] || { echo "ERROR: FASTA not found: $FASTA" >&2; exit 1; }
    [ -s "${FASTA}.fai" ] || { echo "ERROR: ${FASTA}.fai not found" >&2; exit 1; }
}

## Function to list exons (transcript order index, transcript_id, contig, start, end, strand)
## of the contigs present in the genome, sorted by transcript then start
list_exons() {
    zcat -f "$GTF" | awk -F'\t' 'NR == FNR { ok[$1] = 1; next }
        $3 == "exon" && ($1 in ok) {
            if (!match($9, /transcript_id "[^"]+"/)) next
            t = substr($9, RSTART + 15, RLENGTH - 16)
            if (!(t in idx)) idx[t] = ++n
            print idx[t] "\t" t "\t" $1 "\t" $4 "\t" $5 "\t" $7
        }' "${FASTA}.fai" - \
    | LC_ALL=C sort -t$'\t' -k1,1n -k4,4n -S 1G -T "$TMP" > "${TMP}/exons.tsv"
    [ -s "${TMP}/exons.tsv" ] || { echo "ERROR: no exon of the GTF is on a contig of $FASTA (contig names differ?)" >&2; exit 1; }
}

## Function to fetch each distinct exon once per strand (minus strand reverse-complemented)
fetch_exons() {
    awk -F'\t' '$6 != "-" { print $3 ":" $4 "-" $5 }' "${TMP}/exons.tsv" | LC_ALL=C sort -u -S 1G -T "$TMP" > "${TMP}/plus.txt"
    awk -F'\t' '$6 == "-" { print $3 ":" $4 "-" $5 }' "${TMP}/exons.tsv" | LC_ALL=C sort -u -S 1G -T "$TMP" > "${TMP}/minus.txt"
    : > "${TMP}/plus.fa"; : > "${TMP}/minus.fa"
    [ ! -s "${TMP}/plus.txt" ]  || samtools faidx "$FASTA" -r "${TMP}/plus.txt" > "${TMP}/plus.fa"
    [ ! -s "${TMP}/minus.txt" ] || samtools faidx -i "$FASTA" -r "${TMP}/minus.txt" > "${TMP}/minus.fa"
}

## Function to join the exons of each transcript (minus strand: last exon first) and write the FASTA
assemble() {
    awk -F'\t' '
        FILENAME == ARGV[1] || FILENAME == ARGV[2] {
            if (/^>/) { k = substr($0, 2); sub(/\/rc$/, "", k); k = k (FILENAME == ARGV[2] ? "-" : "+"); next }
            s[k] = s[k] $0; next
        }
        function flush() { if (cur != "") { print ">" cur; for (i = 1; i <= length(seq); i += 60) print substr(seq, i, 60) } }
        {
            if ($2 != cur) { flush(); cur = $2; seq = "" }
            k = $3 ":" $4 "-" $5 ($6 == "-" ? "-" : "+")
            if (!(k in s)) { print "ERROR: no sequence for exon " k " of " $2 > "/dev/stderr"; bad = 1; exit 1 }
            seq = ($6 == "-") ? s[k] seq : seq s[k]
        }
        END { if (!bad) flush() }' "${TMP}/plus.fa" "${TMP}/minus.fa" "${TMP}/exons.tsv" > "${OUT}.tmp"
    mv "${OUT}.tmp" "$OUT"
}

##### Main Script
parse_args "$@"
TMP=$(mktemp -d "${TMPDIR:-.}/gtf2tx.XXXXXX"); trap 'rm -rf "$TMP"' EXIT
list_exons
fetch_exons
assemble
echo "transcripts: $(grep -c '^>' "$OUT") -> $OUT" >&2
