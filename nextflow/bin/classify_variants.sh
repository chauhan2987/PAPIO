#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Computational pre-classification of PASS variants for review (somatic: CGC role-aware
#          oncogenicity points; germline: ClinVar P/LP, BA1/BS1, LoF in germline TSGs).
#          Extracts one consequence per variant with bcftools +split-vep and scores it with
#          classify_variants.awk. Called by the PAPIO process CLASSIFY_VARIANTS.
# Usage  : classify_variants.sh --sample ID --vcf VEP.vcf.gz --mode somatic|germline [--cgc CGC.tsv[.gz]]
#                   [--af-thresholds GENE_BA1_BS1.tsv] --af-field MAX_AF --af-cutoff F
#                   --cosmic-min N --cosmic-strong N --min-stars N --ba1 F --germ-low-vaf F --mosaic-min-dp N
# Output : <sample>_preclassification.tsv
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; VCF=""; MODE="somatic"; CGC=""; AF_THR=""; AF_FIELD="MAX_AF"; AF_CUTOFF="0.01"; COSMIC_MIN="10"; COSMIC_STRONG="50"; MIN_STARS="1"; BA1="0.05"; GERM_LOW_VAF="0.30"; MOSAIC_MIN_DP="20"

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "classify_variants.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)           SAMPLE=$2 ;;
        --vcf)              VCF=$2 ;;
        --mode)             MODE=$2 ;;
        --cgc)              CGC=$2 ;;
        --af-thresholds)    AF_THR=$2 ;;
        --af-field)         AF_FIELD=$2 ;;
        --af-cutoff)        AF_CUTOFF=$2 ;;
        --cosmic-min)       COSMIC_MIN=$2 ;;
        --cosmic-strong)    COSMIC_STRONG=$2 ;;
        --min-stars)        MIN_STARS=$2 ;;
        --ba1)              BA1=$2 ;;
        --germ-low-vaf)     GERM_LOW_VAF=$2 ;;
        --mosaic-min-dp)    MOSAIC_MIN_DP=$2 ;;
        *) echo "classify_variants.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

# CGC gene -> role/tier, columns located by HEADER NAME (layout differs between
# CGC releases; fixed indices fail silently). Missing CGC -> empty lookup.
if [ -n "${CGC}" ]; then
    case "${CGC}" in *.gz) CATC="zcat" ;; *) CATC="cat" ;; esac
    "$CATC" "${CGC}" | awk -F'\t' '
        NR==1 { for (i = 1; i <= NF; i++) { h = toupper($i); gsub(/[ ]+/, "_", h)
                    if (h == "GENE_SYMBOL")    g = i
                    if (h == "ROLE_IN_CANCER") r = i
                    if (h == "TIER")           t = i
                    if (h == "GERMLINE")       m = i }
                if (!g || !r || !t) { print "ERROR: CGC header lacks GENE_SYMBOL/ROLE_IN_CANCER/TIER" > "/dev/stderr"; bad = 1; exit 1 }
                next }
        $g != "" { gl = !m ? "NA" : (tolower($m) ~ /^y/ ? "Y" : "N")
                   print $g "\t" $r "\t" $t "\t" gl }
        END { if (bad) exit 1 }' > cgc_lookup.tsv
else
    : > cgc_lookup.tsv
fi

# Gene-specific germline BA1/BS1 (ClinGen VCEP specifications): GENE<TAB>BA1<TAB>BS1.
if [ -n "${AF_THR}" ]; then
    awk -F'\t' '/^#/ || NF == 0 || toupper($1) == "GENE" { next }
        NF < 3 || $2 !~ /^[0-9.eE-]+$/ || $3 !~ /^[0-9.eE-]+$/ || $2 + 0 < $3 + 0 {
            print "ERROR: germline_af_thresholds line " NR " is not GENE<TAB>BA1<TAB>BS1 with BA1 >= BS1" > "/dev/stderr"; bad = 1; exit 1 }
        { print toupper($1) "\t" $2 "\t" $3 }
        END { if (bad) exit 1 }' "${AF_THR}" > af_thr.tsv
else
    : > af_thr.tsv
fi

fields="$(bcftools +split-vep -l "${VCF}" 2>/dev/null | awk '{print $2}')"

# Every optional field gets a column ('.' when absent) so column positions never shift.
if grep -qx "COSMIC_GENOME_SCREEN_SAMPLE_COUNT" <<< "$fields"; then
    COSMIC_FMT='\t%COSMIC_GENOME_SCREEN_SAMPLE_COUNT'
else
    COSMIC_FMT='\t.'
fi

# ClinVar oncogenicity classification (somatic; ClinVar releases from 2024).
if grep -qx "ClinVar_ONC" <<< "$fields"; then
    ONC_FMT='\t%ClinVar_ONC'
else
    echo "WARN: ClinVar_ONC not in the VEP annotation (older ClinVar release); no ClinVar oncogenicity evidence" >&2
    ONC_FMT='\t.'
fi
if grep -qx "ClinVar_ONCREVSTAT" <<< "$fields"; then
    ONC_FMT="${ONC_FMT}\t%ClinVar_ONCREVSTAT"
else
    ONC_FMT="${ONC_FMT}\t."
fi

if grep -qx "ClinVar_CLNSIG" <<< "$fields"; then
    CLIN_FMT='%ClinVar_CLNSIG'
elif grep -qx "CLIN_SIG" <<< "$fields"; then
    echo "WARN: ClinVar_CLNSIG not in VCF; using VEP cache CLIN_SIG" >&2
    CLIN_FMT='%CLIN_SIG'
else
    echo "WARN: no ClinVar field in the VEP annotation; clinical significance = '.'" >&2
    CLIN_FMT='.'
fi

if grep -qx "ClinVar_CLNREVSTAT" <<< "$fields"; then
    REVSTAT_FMT='%ClinVar_CLNREVSTAT'
else
    REVSTAT_FMT='.'
fi

# Population AF: gnomAD grpmax FAF95 when annotated (ClinGen SVI field for BA1/BS1), else
# the VEP cache field (MAX_AF by default, which includes founder and small groups). A VEP
# run without either lacks the field: split-vep would abort, so write '.' instead.
if grep -qx "gnomAD_fafmax_faf95_max_joint" <<< "$fields"; then
    AF_NAME='gnomAD_FAF95_GRPMAX'; AF_FMT='%gnomAD_fafmax_faf95_max_joint'
elif grep -qx "${AF_FIELD}" <<< "$fields"; then
    AF_NAME="${AF_FIELD}"; AF_FMT="%${AF_FIELD}"
else
    echo "WARN: no gnomAD FAF and no ${AF_FIELD} in the VEP annotation; population AF = '.' (no AF-based benign calls)" >&2
    AF_NAME="${AF_FIELD}"; AF_FMT='.'
fi

# Exon number (n/N) for the last-exon (NMD-escape) check on truncating variants.
if grep -qx "EXON" <<< "$fields"; then EXON_FMT='\t%EXON'; else EXON_FMT='\t.'; fi
if grep -qx "CDS_position" <<< "$fields"; then EXON_FMT="${EXON_FMT}\t%CDS_position"; else EXON_FMT="${EXON_FMT}\t."; fi
if grep -qx "NMD" <<< "$fields"; then EXON_FMT="${EXON_FMT}\t%NMD"; else EXON_FMT="${EXON_FMT}\t."; fi

# Sample allele depths (VAF; origin flags). The sample is this run's sample, by exact name.
# A multi-sample VCF without that name (e.g. TUMOR/NORMAL re-entry VCFs) gives no VAF rather
# than a guess that could be the normal; a single-sample VCF is used whatever its name.
bcftools query -l "${VCF}" > samples.txt
SMP="${SAMPLE}"
if ! grep -qxF "$SMP" samples.txt; then
    if [ "$(awk 'END { print NR }' samples.txt)" -eq 1 ]; then
        SMP="$(cat samples.txt)"; echo "WARN: sample ${SAMPLE} not in VCF; VAF from its only sample $SMP" >&2
    else
        SMP=''
    fi
fi
# header read in full first: 'grep -q' closing the pipe early on a large header would
# make the test fail under pipefail and silently drop VAF
vcf_hdr=$(bcftools view -h "${VCF}")
if [ -n "$SMP" ] && grep -q '^##FORMAT=<ID=AD,' <<< "$vcf_hdr"; then
    SEL=(-s "$SMP"); AD_FMT='\t[%AD]'
else
    echo "WARN: no sample or no FORMAT/AD; VAF = '.'" >&2
    SEL=(); AD_FMT='\t.'
fi

# PASS-only: consensus/filter run in 'tag' mode, so failed calls are still in
# the VCF; they must not be pre-classified or reported.
# One consequence per variant: the VEP PICK transcript (MANE-first pick_order) when
# present, else worst-across-transcripts (VCFs annotated before PICK was enabled).
CSQ_FMT="%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t${AF_FMT}\t${CLIN_FMT}\t${REVSTAT_FMT}${COSMIC_FMT}${ONC_FMT}${EXON_FMT}${AD_FMT}\n"
if grep -qx "PICK" <<< "$fields"; then
    bcftools view ${SEL[@]+"${SEL[@]}"} -f PASS,. -Ou "${VCF}" \
    | bcftools +split-vep -d -f "%PICK\t${CSQ_FMT}" - \
    | awk -F'\t' '$1 == "1"' | cut -f2- | sort -u -T . > variants_raw.tsv
else
    echo "WARN: no VEP PICK field; using worst consequence across transcripts" >&2
    bcftools view ${SEL[@]+"${SEL[@]}"} -f PASS,. -Ou "${VCF}" \
    | bcftools +split-vep -s worst -d -f "${CSQ_FMT}" - \
    | sort -u -T . > variants_raw.tsv
fi


awk -F'\t' -v AF_CUTOFF="${AF_CUTOFF}" -v COSMIC_MIN="${COSMIC_MIN}" -v COSMIC_STRONG="${COSMIC_STRONG}" -v MIN_STARS="${MIN_STARS}" \
    -v MODE="${MODE}" -v BA1="${BA1}" \
    -v GERM_LOW_VAF="${GERM_LOW_VAF}" -v MOSAIC_MIN_DP="${MOSAIC_MIN_DP}" \
    -f "$(command -v classify_variants.awk)" variants_raw.tsv > body.tsv

{
    echo "# COMPUTATIONAL PRE-CLASSIFICATION FOR PATHOLOGIST REVIEW - NOT A FINAL AMP/ASCO/CAP CLINICAL TIER"
    echo "# Final tiering requires a curated knowledge base (OncoKB/CIViC) and expert review."
    echo "# ClinVar source: ClinVar_CLNSIG when present, else VEP cache CLIN_SIG (germline; not used for somatic oncogenicity)."
    echo "# CLINVAR_ONC: ClinVar oncogenicity classification (somatic; star-gated); used for somatic oncogenicity."
    echo "# Rules (${MODE}): somatic = CGC role-aware, ClinGen/CGC/VICC-style points (COSMIC >= ${COSMIC_MIN} recurrent, >= ${COSMIC_STRONG} hotspot); germline = ClinVar P/LP (star-gated) + BA1 AF > ${BA1}; LoF review in germline TSGs."
    echo "# Population AF column: ${AF_NAME} (gnomAD_FAF95_GRPMAX = gnomAD v4.1 joint grpmax filtering AF)."
    echo "# LOW_IMPACT = no coding/splice consequence (VEP LOW/MODIFIER) and no other evidence: not assessed, not a benign call."
    echo "# VAF from sample ${SMP:-none} FORMAT/AD. ORIGIN (${MODE}): flag only, does not change PRE_CLASSIFICATION."
    if [ -s af_thr.tsv ]; then echo "# Germline BA1/BS1: gene-specific for $(awk 'END { print NR }' af_thr.tsv) genes (germline_af_thresholds); else ${BA1}/${AF_CUTOFF}."; fi
    printf 'CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tIMPACT\t%s\tClinVar\tCLINVAR_STARS\tCOSMIC_RECURRENCE\tCGC_ROLE\tCGC_TIER\tPRE_CLASSIFICATION\tCLINVAR_ONC\tCLINVAR_ONC_STARS\tVAF\tORIGIN\n' "${AF_NAME}"
    cat body.tsv
} > ${SAMPLE}_preclassification.tsv
