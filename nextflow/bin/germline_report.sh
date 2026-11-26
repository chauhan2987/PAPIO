#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Consolidated germline report for one sample: ClinVar-flagged SNV/indels, LoF in
#          reportable genes (CGC GERMLINE + user list), AnnotSV SV and CNV rows (ACMG 4/5 or
#          reportable gene) and runs of homozygosity. Called by the PAPIO process
#          CONSOLIDATE_GERMLINE_REPORT.
# Usage  : germline_report.sh --sample ID --clin-min pathogenic --lof-max-af F --min-stars N
#                   [--cgc CGC.tsv[.gz]] [--gene-list genes.txt]
#                   Inputs are read from the working directory: snv.vcf.gz, sv_annotsv.tsv, cnv_annotsv.tsv, roh.txt
# Output : <sample>_germline_consolidated.tsv
#------------------------------------------------------------------------------
set -uo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; CLIN_MIN="pathogenic"; LOF_MAX_AF="0.01"; MIN_STARS="1"; CGC=""; GENE_LIST=""

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "germline_report.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)        SAMPLE=$2 ;;
        --clin-min)      CLIN_MIN=$2 ;;
        --lof-max-af)    LOF_MAX_AF=$2 ;;
        --min-stars)     MIN_STARS=$2 ;;
        --cgc)           CGC=$2 ;;
        --gene-list)     GENE_LIST=$2 ;;
        *) echo "germline_report.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

REPORT=${SAMPLE}_germline_consolidated.tsv

{
  echo "## ============================================================"
  echo "## CONSOLIDATED GERMLINE REPORT - SAMPLE: ${SAMPLE}"
  echo "## Generated: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "## NOTE: computational summary for review; not a final ACMG clinical classification."
  echo "## ============================================================"
} > $REPORT

# -- [1] SNV/indel - pull ClinVar-flagged calls from the VEP VCF ----------
# Same sources as CLASSIFY_VARIANTS: ClinVar from the dated ClinVar VCF (ClinVar_CLNSIG) when
# present, else the VEP cache CLIN_SIG; population AF = gnomAD grpmax filtering AF when
# annotated, else MAX_AF; one consequence per variant = the VEP PICK transcript when present.
csq_names=$(bcftools +split-vep -l snv.vcf.gz 2>/dev/null | awk '{ print $2 }' || true)
has() { grep -qx "$1" <<< "$csq_names"; }
if   has ClinVar_CLNSIG; then clin_tag='%ClinVar_CLNSIG'; clin_src="ClinVar VCF (ClinVar_CLNSIG)"
elif has CLIN_SIG;       then clin_tag='%CLIN_SIG';       clin_src="VEP cache CLIN_SIG"
else clin_tag=''; clin_src="none"; fi
## a missing field would make split-vep fail and hide behind "(query failed)": write '.'
if   has gnomAD_fafmax_faf95_max_j; then af_tag='%gnomAD_fafmax_faf95_max_j'; af_name='gnomAD_FAF95_GRPMAX'
elif has MAX_AF;                        then af_tag='%MAX_AF';                        af_name='MAX_AF'
else                                         af_tag='.';                              af_name='MAX_AF'; fi
if has ClinVar_CLNREVSTAT; then rev_tag='%ClinVar_CLNREVSTAT'; else rev_tag='.'; fi
if has EXON;         then tag_EXON='%EXON';                 else tag_EXON='.'; fi
if has CDS_position; then tag_CDS_position='%CDS_position'; else tag_CDS_position='.'; fi
if has NMD;          then tag_NMD='%NMD';                   else tag_NMD='.'; fi

## one consequence per variant: PICK transcript, else the worst across transcripts
csq_rows() {
    if has PICK; then
        bcftools view -f PASS,. -Ou snv.vcf.gz 2>/dev/null \
        | bcftools +split-vep -d -f "%PICK\t$1" - 2>/dev/null | awk -F'\t' '$1 == "1"' | cut -f2-
    else
        bcftools view -f PASS,. -Ou snv.vcf.gz 2>/dev/null \
        | bcftools +split-vep -s worst -d -f "$1" - 2>/dev/null
    fi
}

printf '\n## [1] GERMLINE SNV/INDEL  (VEP+ClinVar; significance matching "%s"; source: %s)\n' "${CLIN_MIN}" "$clin_src" >> $REPORT
if [ -s snv.vcf.gz ]; then
  printf 'CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tIMPACT\t%s\tClinVar\tGT\n' "$af_name" >> $REPORT
  if [ -n "$clin_tag" ]; then
    # PASS-only; 'conflicting...pathogenicity' contains the substring 'pathogenic' and
    # must not be reported as pathogenic.
    csq_rows "%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t${af_tag}\t${clin_tag}\t[%GT]\t${rev_tag}\n" \
      | awk -F'\t' -v term="${CLIN_MIN}" -v minst="${MIN_STARS}" -v hasrev="$rev_tag" 'BEGIN { OFS = "\t" }
          function stars(s) { s = tolower(s)
              if (s == "" || s == "." || s ~ /no_assertion|no_classification/) return 0
              if (s ~ /practice_guideline/) return 4; if (s ~ /expert_panel/) return 3
              if (s ~ /multiple_submitters/) return 2; return 1 }
          tolower($9) ~ tolower(term) && tolower($9) !~ /conflicting/ &&
          (hasrev == "." || stars($11) >= minst + 0) { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }' \
      | sort -u >> $REPORT \
      || echo "(SNV query failed)" >> $REPORT
  else
    echo "(VEP VCF has no ClinVar_CLNSIG / CLIN_SIG field; run ANNOTATE_VCF with ClinVar)" >> $REPORT
  fi
else echo "(no annotated SNV/indel VCF)" >> $REPORT; fi

# -- [1b] LoF in reportable genes (tumour-agnostic, parameter-driven gene set) -----
#   gene set = CGC genes flagged GERMLINE (params.cosmic_cgc_tsv) + params.germline_report_genes.
#   NMD_ESCAPE as in CLASSIFY_VARIANTS (ClinGen SVI PVS1): a nonsense/frameshift in the last
#   exon, or in the last 50 nt of the penultimate exon (NMD plugin), of a multi-exon transcript.
#   Not LOFTEE-filtered - every call needs review.
: > report_genes.txt
if [ -n "${CGC}" ]; then
    case "${CGC}" in *.gz) CATC="zcat" ;; *) CATC="cat" ;; esac
    "$CATC" "${CGC}" | awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) { h = toupper($i); gsub(/[ ]+/, "_", h)
                      if (h == "GENE_SYMBOL") g = i; if (h == "GERMLINE") m = i }
                  next }
        g && m && tolower($m) ~ /^y/ && $g != "" { print toupper($g) }' >> report_genes.txt || true
fi
if [ -n "${GENE_LIST}" ]; then
    grep -v -E '^[[:space:]]*(#|$)' "${GENE_LIST}" | awk '{ print toupper($1) }' >> report_genes.txt || true
fi
sort -u -o report_genes.txt report_genes.txt
n_genes=$(grep -c . report_genes.txt || true)
printf '\n## [1b] GERMLINE LoF (HIGH impact) IN REPORTABLE GENES  (%s genes; %s < %s; not LOFTEE-filtered)\n' "$n_genes" "$af_name" "${LOF_MAX_AF}" >> $REPORT
if [ "$n_genes" -eq 0 ]; then
    echo "(no reportable gene set: set params.germline_report_genes and/or a CGC file with a GERMLINE column)" >> $REPORT
elif [ -s snv.vcf.gz ]; then
    printf 'CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tIMPACT\t%s\tClinVar\tGT\tNMD_ESCAPE\n' "$af_name" >> $REPORT
    csq_rows "%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t${af_tag}\t${clin_tag:-.}\t[%GT]\t${tag_EXON}\t${tag_CDS_position}\t${tag_NMD}\n" \
      | awk -F'\t' -v maxaf="${LOF_MAX_AF}" 'BEGIN { OFS = "\t"; while ((getline g < "report_genes.txt") > 0) G[g] = 1 }
          $7 == "HIGH" && $6 ~ /(stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|transcript_ablation)/ &&
          (toupper($5) in G) && ($8 == "" || $8 == "." || $8 + 0 < maxaf + 0) {
              ne = split($11, ex, "/"); split($12, cp, "-")
              if (ne == 2) { nr = split(ex[1], er, "-"); ex[1] = er[nr] }
              if ($6 !~ /(stop_gained|frameshift_variant)/) esc = "NA"
              else if (ne != 2) esc = "NA"
              else if (ex[2] + 0 <= 1) esc = "no"
              else if (ex[1] + 0 == ex[2] + 0) esc = "yes"
              else if (ex[1] + 0 == ex[2] - 1 && $13 == "NMD_escaping_variant" && cp[1] ~ /^[0-9]+$/ && cp[1] + 0 > 101) esc = "yes"
              else esc = "no"
              print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, esc }' \
      | sort -u >> $REPORT \
      || echo "(LoF query failed)" >> $REPORT
else echo "(no annotated SNV/indel VCF)" >> $REPORT; fi

# AnnotSV rows worth reading: one row per SV (Annotation_mode=full) that passed the caller's
# filters (FILTER PASS or '.', as in the somatic report), kept when ACMG class 4/5 OR
# overlapping a reportable gene. Unknown layout -> all rows.
filter_annotsv() {
    awk -F'\t' 'BEGIN { while ((getline g < "report_genes.txt") > 0) G[g] = 1 }
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; print; next }
        { if (("Annotation_mode" in h) && $(h["Annotation_mode"]) != "full") next
          if (("FILTER" in h) && $(h["FILTER"]) != "PASS" && $(h["FILTER"]) != ".") next
          keep = !(("ACMG_class" in h) || ("Gene_name" in h))
          if (("ACMG_class" in h) && $(h["ACMG_class"]) ~ /(^|=)(4|5)$/) keep = 1
          if (!keep && ("Gene_name" in h)) {
              n = split(toupper($(h["Gene_name"])), gs, /[;,]/)
              for (j = 1; j <= n; j++) if (gs[j] in G) { keep = 1; break } }
          if (keep) print }' "$1"
}

# -- [2] SV - AnnotSV (Manta diploidSV), ACMG class 4/5 highlighted ------
printf '\n## [2] GERMLINE SV  (AnnotSV, Manta diploidSV; PASS, full rows, ACMG 4/5 or reportable gene)\n' >> $REPORT
if [ -s sv_annotsv.tsv ] && [ "$(wc -l < sv_annotsv.tsv)" -gt 1 ]; then
  filter_annotsv sv_annotsv.tsv >> $REPORT
else echo "(no annotated SV output)" >> $REPORT; fi

# -- [3] CNV - AnnotSV (gCNV segments) -----------------------------------
printf '\n## [3] GERMLINE CNV  (AnnotSV, GATK gCNV calls (non-reference segments); full rows, ACMG 4/5 or reportable gene)\n' >> $REPORT
if [ -s cnv_annotsv.tsv ] && [ "$(wc -l < cnv_annotsv.tsv)" -gt 1 ]; then
  filter_annotsv cnv_annotsv.tsv >> $REPORT
else echo "(no annotated CNV output)" >> $REPORT; fi

# -- [4] ROH - bcftools roh region records -------------------------------
printf '\n## [4] RUNS OF HOMOZYGOSITY  (bcftools roh, region-level)\n' >> $REPORT
if [ -s roh.txt ]; then cat roh.txt >> $REPORT; else echo "(no ROH output)" >> $REPORT; fi
