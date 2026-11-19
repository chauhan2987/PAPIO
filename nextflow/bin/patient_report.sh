#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Consolidated somatic report for one tumor: classified SNV/indels, copy number,
#          LOH, structural variants, MSI, TMB and contamination. Called by the PAPIO process
#          CONSOLIDATE_PATIENT_REPORT.
# Usage  : patient_report.sh --patient ID --sample ID --sex 0|1|NA --sex-label female|male|unknown
#                   --max-contam F [--amp-min CN]
#                   Inputs are read from the working directory: class.tsv, cnv.cns, seq_segments.txt,
#                   sv_annot.tsv, msi_summary, tmb.txt, contamination.table (absent = not run)
# Output : <patient>_consolidated.tsv
#------------------------------------------------------------------------------
set -uo pipefail

###### GLOBAL VARIABLES
PATIENT=""; SAMPLE=""; SEX="NA"; SEX_LAB="unknown"; MAX_CONTAM="0.03"; AMP_MIN=""

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "patient_report.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --patient)       PATIENT=$2 ;;
        --sample)        SAMPLE=$2 ;;
        --sex)           SEX=$2 ;;
        --sex-label)     SEX_LAB=$2 ;;
        --max-contam)    MAX_CONTAM=$2 ;;
        --amp-min)       AMP_MIN=$2 ;;
        *) echo "patient_report.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

REPORT=${PATIENT}_consolidated.tsv

{
  echo "## ============================================================"
  echo "## CONSOLIDATED VARIANT REPORT - PATIENT: ${PATIENT}   (tumor sample: ${SAMPLE})"
  echo "## Generated: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "## NOTE: computational summary for review; not a final clinical tier."
  echo "## ============================================================"
} > $REPORT

# true when a tabular input has at least one data row after its header line
# (a header-only Sequenza / CNVkit / contamination file means "no result").
has_rows() { [ -s "$1" ] && awk 'NR > 1 { f = 1; exit } END { exit !f }' "$1"; }

printf '\n## [1] SOMATIC SNV/INDEL  (07_classification; PASS calls only)\n' >> $REPORT
if [ -s class.tsv ]; then grep -v '^#' class.tsv >> $REPORT; else echo "(no classification output)" >> $REPORT; fi

# CN calls relative to the EXPECTED copy number: autosomes = tumour ploidy baseline,
# chrX/chrY sex-aware (male X/Y = baseline/2; female Y skipped; unknown sex -> X/Y skipped).
# Source: Sequenza (purity/ploidy-fitted) when present, else CNVkit (NOT purity-adjusted).
if has_rows seq_segments.txt; then
  printf '\n## [2] COPY NUMBER  (Sequenza segments; purity/ploidy-fitted; calls vs sex-aware expected CN)\n' >> $REPORT
  awk -v src=sequenza -v sex="${SEX}" -v sexlab="${SEX_LAB}" -v ampmin="${AMP_MIN}" -f "$(command -v cn_calls.awk)" seq_segments.txt seq_segments.txt >> $REPORT
elif has_rows cnv.cns; then
  printf '\n## [2] COPY NUMBER  (CNVkit call.cns; NOT purity-adjusted; calls vs sex-aware expected CN)\n' >> $REPORT
  awk -v src=cnvkit -v sex="${SEX}" -v sexlab="${SEX_LAB}" -v ampmin="${AMP_MIN}" -f "$(command -v cn_calls.awk)" cnv.cns >> $REPORT
else
  printf '\n## [2] COPY NUMBER\n(no Sequenza or CNVkit output)\n' >> $REPORT
fi

printf '\n## [3] LOH  (Sequenza; minor allele B==0; male chrX/chrY = HEMIZYGOUS, expected; unknown sex = SEX_UNKNOWN)\n' >> $REPORT
if has_rows seq_segments.txt; then
  awk -F'\t' -v sex="${SEX}" 'BEGIN{OFS="\t"; print "CHROM","START","END","CNt","A","B","STATUS"}
    NR==1{ for(i=1;i<=NF;i++) h[$i]=i; next }
    { b=$(h["B"]); if(b=="" || b=="NA")next; c=$(h["chromosome"])
      xy = (c ~ /^(chr)?[XY]$/)
      st = (xy && sex == "1") ? "HEMIZYGOUS" : (xy && sex != "0") ? "SEX_UNKNOWN" : "LOH"
      if(b+0==0 && $(h["CNt"])!="") print c,$(h["start.pos"]),$(h["end.pos"]),$(h["CNt"]),$(h["A"]),b,st }' seq_segments.txt >> $REPORT
else echo "(no Sequenza output)" >> $REPORT; fi

printf '\n## [4] STRUCTURAL VARIANTS  (Manta + AnnotSV; one row/SV, PASS. ACMG = germline framework, advisory)\n' >> $REPORT
if [ -s sv_annot.tsv ] && [ "$(grep -vc '^AnnotSV_ID' sv_annot.tsv)" -gt 0 ]; then
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;print"CHROM","START","END","SVTYPE","GENE(S)","ACMG";next}
    (("Annotation_mode" in h)?$(h["Annotation_mode"]):"full")!="full"{next}
    (("FILTER" in h)&&$(h["FILTER"])!="PASS"&&$(h["FILTER"])!="."){next}
    {g=(("Gene_name" in h)?$(h["Gene_name"]):".");g=(g==""?".":g)
     a=(("ACMG_class" in h)?$(h["ACMG_class"]):".");a=(a==""||a=="NA"?".":a)
     print $(h["SV_chrom"]),$(h["SV_start"]),$(h["SV_end"]),$(h["SV_type"]),g,a}' OFS='\t' sv_annot.tsv >> $REPORT
else echo "(no annotated SV output)" >> $REPORT; fi

printf '\n## [5] MSI  (MSIsensor-pro)\n' >> $REPORT
if [ -s msi_summary ]; then cat msi_summary >> $REPORT; else echo "(no MSI output)" >> $REPORT; fi

printf '\n## [6] TMB\n' >> $REPORT
if [ -s tmb.txt ]; then cat tmb.txt >> $REPORT; else echo "(no TMB output)" >> $REPORT; fi

printf '\n## [7] TUMOUR CROSS-SAMPLE CONTAMINATION  (GATK CalculateContamination; threshold %s)\n' "${MAX_CONTAM}" >> $REPORT
if has_rows contamination.table; then
    awk -F'\t' -v max="${MAX_CONTAM}" 'NR == 1 { next }
        { printf "sample\t%s\ncontamination\t%s\nerror\t%s\nstatus\t%s\n", $1, $2, $3, ($2 + 0 > max + 0 ? "WARN" : "PASS") }' \
        contamination.table >> $REPORT
else echo "(no contamination estimate - pileup_variants not set, or no result rows)" >> $REPORT; fi
