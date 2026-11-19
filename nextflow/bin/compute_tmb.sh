#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Tumor mutational burden: PASS nonsynonymous coding variants per megabase of the
#          territory = capture BED, intersected with a coding BED when given, and with the bases
#          callable in both tumor and normal when their mosdepth depth classes are given.
#          WGS reports NA. Called by the PAPIO process COMPUTE_TMB.
# Usage  : compute_tmb.sh --sample ID --vcf VEP.vcf.gz --assay wes|wgs --capture BED [--coding CDS.bed]
#                   [--tumor-depth T.quantized.bed.gz --normal-depth N.quantized.bed.gz
#                    --min-tumor-depth N --min-normal-depth N] --cutoff N --include-splice true|false
#          Depth thresholds act on the mosdepth classes (0:1:10:20:30:50:100): a threshold that
#          is not a class boundary is rounded up to the next one.
# Output : <sample>_tmb.txt
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; VCF=""; ASSAY="wes"; CAPTURE=""; CODING=""; CUTOFF="10"; INCL_SPLICE="false"
T_DEPTH=""; N_DEPTH=""; MIN_T_DP="20"; MIN_N_DP="10"

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "compute_tmb.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)            SAMPLE=$2 ;;
        --vcf)               VCF=$2 ;;
        --assay)             ASSAY=$2 ;;
        --capture)           CAPTURE=$2 ;;
        --coding)            CODING=$2 ;;
        --cutoff)            CUTOFF=$2 ;;
        --tumor-depth)      T_DEPTH=$2 ;;
        --normal-depth)     N_DEPTH=$2 ;;
        --min-tumor-depth)  MIN_T_DP=$2 ;;
        --min-normal-depth) MIN_N_DP=$2 ;;
        --include-splice)    INCL_SPLICE=$2 ;;
        *) echo "compute_tmb.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

## nonsynonymous consequences counted (splice donor/acceptor optional)
TERM_RE="missense_variant|stop_gained|stop_lost|start_lost|frameshift_variant|inframe_insertion|inframe_deletion|protein_altering_variant"
[ "$INCL_SPLICE" = "true" ] && TERM_RE="${TERM_RE}|splice_acceptor_variant|splice_donor_variant"

## WGS: the interval BED is the CNVkit access BED, not a coding territory - a coding
## numerator over it is not a TMB. NA in the same layout until a WGS callable-territory
## denominator is implemented.
if [ "$ASSAY" = "wgs" ]; then
    echo "COMPUTE_TMB: assay=wgs - TMB not computed (no valid WGS denominator yet)" >&2
    {   echo -e "sample_id\t${SAMPLE}"
        echo -e "nonsynonymous_coding_variants\tNA"
        echo -e "covered_megabases\tNA"
        echo -e "tmb_mut_per_mb\tNA"
        echo -e "tmb_cutoff\t${CUTOFF}"
        echo -e "tmb_status\tNA"
        echo -e "include_splice\t${INCL_SPLICE}"
        echo -e "territory\tNA"
        echo -e "consequence_source\tNA"
    } > "${SAMPLE}_tmb.txt"
    exit 0
fi
[ -n "$CAPTURE" ] && [ -s "$CAPTURE" ] || { echo "ERROR: COMPUTE_TMB on WES needs the capture BED (interval_bed)" >&2; exit 1; }

# Merge a BED (plain or .gz) into sorted, non-overlapping BED3 on stdout.
merge_bed() {
    case "$1" in *.gz) CAT="zcat" ;; *) CAT="cat" ;; esac
    $CAT "$1" | { grep -v -E '^(#|track|browser)' || true; } \
        | awk -F'\t' 'NF >= 3 { print $1 "\t" $2 "\t" $3 }' \
        | LC_ALL=C sort -k1,1 -k2,2n \
        | awk -F'\t' 'BEGIN { OFS = "\t" }
            { if ($1 == c && $2 <= e) { if ($3 > e) e = $3 }
              else { if (c != "") print c, s, e; c = $1; s = $2; e = $3 } }
            END { if (c != "") print c, s, e }'
}

# ---- TERRITORY (denominator AND numerator region): capture, intersected with
#      the coding BED when provided (sweep over +1/-1 boundary events). ----
# Contig naming (chr1 vs 1, chrM vs MT): bring every BED to the VCF's naming. A mismatch
# would otherwise give an empty territory and a silent NA.
bcftools view -h ${VCF} | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' > vcf.contigs
harmonize() {
    if [ ! -s vcf.contigs ]; then cat "$1"; return; fi
    awk -F'\t' 'BEGIN { OFS = "\t" } NR == FNR { ok[$1] = 1; next }
        { c = $1
          if (!(c in ok)) { d = (c ~ /^chr/) ? substr(c, 4) : "chr" c
                            if (d == "chrMT") d = "chrM"; else if (d == "M") d = "MT"
                            if (d in ok) c = d; else next }
          $1 = c; print }' vcf.contigs "$1" | LC_ALL=C sort -k1,1 -k2,2n
}
# Intersection of two sorted, merged BED3 files (sweep over +1/-1 boundary events).
intersect() {
    {   awk -F'\t' '{ print $1 "\t" $2 "\t1\tA"; print $1 "\t" $3 "\t-1\tA" }' "$1"
        awk -F'\t' '{ print $1 "\t" $2 "\t1\tB"; print $1 "\t" $3 "\t-1\tB" }' "$2"
    } | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
      | awk -F'\t' 'BEGIN { OFS = "\t" }
          $1 != c { c = $1; a = 0; b = 0; open = 0 }
          { if ($4 == "A") a += $3; else b += $3
            if (!open && a > 0 && b > 0) { st = $2; open = 1 }
            else if (open && !(a > 0 && b > 0)) { if ($2 > st) print c, st, $2; open = 0 } }'
}
# Bases of a mosdepth quantized BED whose depth class starts at >= the threshold.
callable() {
    zcat -f "$1" | awk -F'\t' -v m="$2" '{ split($4, q, ":"); if (q[1] + 0 >= m + 0) print $1 "\t" $2 "\t" $3 }' > callable.tmp
    merge_bed callable.tmp
}
merge_bed "${CAPTURE}" > capture.raw.bed
harmonize capture.raw.bed > capture.merged.bed
cp capture.merged.bed territory.bed
if [ -n "${CODING}" ]; then
    merge_bed "${CODING}" > coding.raw.bed
    harmonize coding.raw.bed > coding.merged.bed
    intersect territory.bed coding.merged.bed > t.tmp && mv t.tmp territory.bed
    territory="capture_x_coding_bed"
else
    territory="capture_bed_only(includes_non-coding_target;set_params.tmb_coding_bed)"
fi
if [ -n "${T_DEPTH}" ] && [ -n "${N_DEPTH}" ]; then
    callable "${T_DEPTH}" "${MIN_T_DP}" > t_call.raw.bed; harmonize t_call.raw.bed > t_call.bed
    callable "${N_DEPTH}" "${MIN_N_DP}" > n_call.raw.bed; harmonize n_call.raw.bed > n_call.bed
    bed_mb=$(awk -F'\t' '{ t += $3 - $2 } END { printf "%.6f", t / 1e6 }' territory.bed)
    intersect territory.bed t_call.bed > t.tmp && intersect t.tmp n_call.bed > territory.bed
    territory="${territory}_x_callable(tumor_dp>=${MIN_T_DP},normal_dp>=${MIN_N_DP};bed_mb=${bed_mb})"
    if [ ! -s territory.bed ]; then
        ## the BED territory exists but no base reaches the depth thresholds: a failed sample, not a TMB of 0
        echo "WARN: no territory base has tumor depth >= ${MIN_T_DP} and normal depth >= ${MIN_N_DP} - TMB NA" >&2
        {   echo -e "sample_id\t${SAMPLE}"
            echo -e "nonsynonymous_coding_variants\tNA"
            echo -e "covered_megabases\t0"
            echo -e "tmb_mut_per_mb\tNA"
            echo -e "tmb_cutoff\t${CUTOFF}"
            echo -e "tmb_status\tNA"
            echo -e "include_splice\t${INCL_SPLICE}"
            echo -e "territory\t${territory}"
            echo -e "consequence_source\tNA"
        } > "${SAMPLE}_tmb.txt"
        exit 0
    fi
fi
if [ ! -s territory.bed ]; then
    echo "ERROR: TMB territory is empty - interval_bed / tmb_coding_bed share no contig (or no overlap) with the VCF." >&2
    echo "       VCF contigs: $(head -3 vcf.contigs | tr '\n' ' ')  capture: $(cut -f1 capture.raw.bed | sort -u | head -3 | tr '\n' ' ')" >&2
    exit 1
fi
covered_mb=$(awk -F'\t' '{ t += $3 - $2 } END { printf "%.6f", t / 1e6 }' territory.bed)

# ---- NUMERATOR: PASS nonsynonymous variants inside the territory. One consequence
#      per variant: VEP PICK transcript (MANE-first pick_order) when annotated,
#      else worst-across-transcripts (VCFs annotated before PICK was enabled). ----
csq_fields=$(bcftools +split-vep -l ${VCF} | awk '{ print $2 }')
if grep -qx PICK <<< "$csq_fields"; then
    csq_src="VEP_PICK"
    # PICK is per allele+gene: count each VARIANT once if any picked consequence qualifies.
    nonsyn=$(bcftools view -f PASS,. -T territory.bed -Ou ${VCF} \
        | bcftools +split-vep -d -f '%PICK\t%CHROM\t%POS\t%REF\t%ALT\t%Consequence\n' - \
        | awk -F'\t' -v re="(${TERM_RE})" '$1 == "1" && $6 ~ re { print $2 "\t" $3 "\t" $4 "\t" $5 }' \
        | sort -u | wc -l)
else
    csq_src="worst_across_transcripts"
    nonsyn=$(bcftools view -f PASS,. -T territory.bed -Ou ${VCF} \
        | bcftools +split-vep -s worst -f '%Consequence\n' - \
        | grep -E -c "(${TERM_RE})" || true)
fi

# ---- TMB = numerator / denominator ----
tmb=$(awk -v n="$nonsyn" -v mb="$covered_mb" 'BEGIN{ if(mb>0) printf "%.4f", n/mb; else print "NA" }')
status=$(awk -v t="$tmb" -v c="${CUTOFF}" 'BEGIN{ if(t=="NA"){print "NA"} else if(t+0>=c+0){print "TMB-High"} else {print "TMB-Low"} }')

{   echo -e "sample_id\t${SAMPLE}"
    echo -e "nonsynonymous_coding_variants\t$nonsyn"
    echo -e "covered_megabases\t$covered_mb"
    echo -e "tmb_mut_per_mb\t$tmb"
    echo -e "tmb_cutoff\t${CUTOFF}"
    echo -e "tmb_status\t$status"
    echo -e "include_splice\t${INCL_SPLICE}"
    echo -e "territory\t$territory"
    echo -e "consequence_source\t$csq_src"
} > ${SAMPLE}_tmb.txt
