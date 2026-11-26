#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Somatic SNV/indel call set from Mutect2, Strelka2, MuSE and VarScan2 (Called by the
#          PAPIO process ENSEMBLE_CONSENSUS). Each caller is reduced to its somatic PASS calls.
#          Backbone mutect2: Mutect2 PASS records carry an N-of-M vote. Backbone union: true N-of-M
#          vote with Mutect2 as one vote; Mutect2-filtered records outvoted by the others are kept
#          (FILTER in M2_FILTER) and sites Mutect2 never emitted are added in the Mutect2 layout.
#          Without --mutect2 (union only) the VCF has its own header (normal, tumour; GT:AD:AF:DP).
#          One caller with threshold 1 (somatic_consensus: false) = that caller's calls, no vote.
#          A disabled caller (no file) lowers the threshold; one that ran and called nothing votes no.
#          Sub-threshold sites are LOWCONSENSUS (tag) or dropped (remove).
# Usage  : ensemble_consensus.sh --sample ID --normal ID --fasta REF [--mutect2 VCF]
#                   [--strelka-snv VCF --strelka-indel VCF] [--muse VCF] [--varscan-snv VCF --varscan-indel VCF]
#                   [--prefix NAME] --min-snv N --min-indel N --muse-exclude Tier5 --filter-mode tag|remove --backbone mutect2|union
# Output : <prefix>.vcf.gz (+ .tbi), <prefix>_stats.txt; prefix default <sample>_consensus (--prefix)
#------------------------------------------------------------------------------
set -euo pipefail

###### GLOBAL VARIABLES
SAMPLE=""; NORMAL=""; FASTA=""; MUTECT2=""; STRELKA_SNV=""; STRELKA_INDEL=""; MUSE=""; VARSCAN_SNV=""; VARSCAN_INDEL=""; PREFIX=""; MIN_SNV="2"; MIN_INDEL="2"; MUSE_EXCLUDE="Tier5"; FILTER_MODE="tag"; BACKBONE="mutect2"

##### Main Script
while [ $# -gt 0 ]; do
    ## every option takes a value; a missing one would shift the rest of the line
    if [ $# -lt 2 ] || [ "${2#--}" != "$2" ]; then echo "ensemble_consensus.sh: option $1 needs a value" >&2; exit 2; fi
    case "$1" in
        --sample)           SAMPLE=$2 ;;
        --normal)           NORMAL=$2 ;;
        --fasta)            FASTA=$2 ;;
        --mutect2)          MUTECT2=$2 ;;
        --strelka-snv)      STRELKA_SNV=$2 ;;
        --strelka-indel)    STRELKA_INDEL=$2 ;;
        --muse)             MUSE=$2 ;;
        --varscan-snv)      VARSCAN_SNV=$2 ;;
        --varscan-indel)    VARSCAN_INDEL=$2 ;;
        --prefix)           PREFIX=$2 ;;
        --min-snv)          MIN_SNV=$2 ;;
        --min-indel)        MIN_INDEL=$2 ;;
        --muse-exclude)     MUSE_EXCLUDE=$2 ;;
        --filter-mode)      FILTER_MODE=$2 ;;
        --backbone)         BACKBONE=$2 ;;
        *) echo "ensemble_consensus.sh: unknown option $1" >&2; exit 2 ;;
    esac
    shift 2
done

    PREFIX=${PREFIX:-${SAMPLE}_consensus}

    case "${BACKBONE}" in
        mutect2|union) ;;
        *) echo "ERROR: consensus_backbone must be 'mutect2' or 'union' (got '${BACKBONE}')" >&2; exit 1 ;;
    esac
    NO_M2=0; [ -s "${MUTECT2}" ] || NO_M2=1
    if [ "$NO_M2" -eq 1 ] && [ "${BACKBONE}" != "union" ]; then
        echo "ERROR: without a Mutect2 VCF only the union backbone can vote" >&2; exit 1
    fi

    # Select one caller VCF (PASS / somatic-filtered), split multiallelics and
    # left-align, so alleles match across callers. Writes BCF to $3.
    select_norm() {
        local invcf="$1" mode="$2" out="$3" sel vs_hdr
        case "$mode" in
          muse)    sel=(bcftools view -e "FILTER==\"${MUSE_EXCLUDE}\"" -O u "$invcf") ;;
          varscan)
              # capture the header first: 'grep -q' closing the pipe early under
              # pipefail (large contig header) would make this test spuriously false.
              vs_hdr=$(bcftools view -h "$invcf")
              # VarScan declares SS as Type=String; compare as a string then (a numeric
              # compare on a String tag is a bcftools error), as a number otherwise.
              if grep -q 'ID=SS,.*Type=String' <<< "$vs_hdr"; then
                  sel=(bcftools view -i 'INFO/SS=="2"' -f PASS,. -O u "$invcf")
              elif grep -q 'ID=SS,' <<< "$vs_hdr"; then
                  sel=(bcftools view -i 'INFO/SS==2' -f PASS,. -O u "$invcf")
              else
                  sel=(bcftools view -f PASS,. -O u "$invcf")
              fi ;;
          *)       sel=(bcftools view -f PASS,. -O u "$invcf") ;;
        esac
        "${sel[@]}" | bcftools norm -m -any -f ${FASTA} -O b -o "$out" -
    }

    # Per-variant caller support (CHROM POS REF ALT CALLER). Selected BCFs are
    # kept (sel.list) so the union backbone can take records from them.
    : > support.raw; : > sel.list
    add_caller() {
        local invcf="$1" caller="$2" mode="$3" out
        [ -s "$invcf" ] || return 0
        out="sel_$(wc -l < sel.list)_${caller}.bcf"
        select_norm "$invcf" "$mode" "$out"
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t${caller}\n" "$out" >> support.raw
        printf '%s\t%s\n' "$out" "$caller" >> sel.list
    }
    add_caller "${MUTECT2}"   mutect2  pass
    add_caller "${STRELKA_SNV}"   strelka2 pass
    add_caller "${STRELKA_INDEL}" strelka2 pass
    add_caller "${MUSE}"      muse     muse
    add_caller "${VARSCAN_SNV}"   varscan  varscan
    add_caller "${VARSCAN_INDEL}" varscan  varscan

    # aggregate - CHROM POS REF ALT CALLERS(csv) NCALLERS(distinct)
    sort -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 support.raw | awk -F'\t' 'BEGIN{OFS="\t"}
      { k=$1 FS $2 FS $3 FS $4
        if(k!=prev && prev!=""){ print pc,pp,pr,pa,cs,nc; cs="";nc=0 }
        if($5!="" && index(","cs",", ","$5",")==0){ cs=(cs==""?$5:cs","$5); nc++ }
        prev=k;pc=$1;pp=$2;pr=$3;pa=$4 }
      END{ if(prev!="") print pc,pp,pr,pa,cs,nc }' > agg.tsv

    # Mutect2 records used as backbone: PASS only ('mutect2'), or every record ('union',
    # so a site Mutect2 filtered can still be kept when the other callers outvote it).
    if [ "$NO_M2" -eq 1 ]; then
        :
    elif [ "${BACKBONE}" = "union" ]; then
        bcftools norm -m -any -f ${FASTA} -O z -o m2.norm.vcf.gz ${MUTECT2}
    else
        bcftools view -f PASS,. -O u ${MUTECT2} \
          | bcftools norm -m -any -f ${FASTA} -O z -o m2.norm.vcf.gz -
    fi
    [ "$NO_M2" -eq 1 ] || bcftools index -f -t m2.norm.vcf.gz

    # Mutect2 MNPs (equal-length REF/ALT > 1 bp) are kept as ONE record (correct codon-level
    # annotation). Other callers report the constituent SNVs, so an MNP's vote = callers that
    # support EVERY constituent SNV (intersection), + mutect2 only if Mutect2 PASSed it.
    : > m2_mnp.tsv
    [ "$NO_M2" -eq 1 ] || bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n' m2.norm.vcf.gz \
      | awk -F'\t' 'BEGIN{OFS="\t"} length($3) > 1 && length($3) == length($4) { print $1, $2, $3, $4, (($5 == "PASS" || $5 == ".") ? 1 : 0) }' > m2_mnp.tsv
    awk -F'\t' 'BEGIN { OFS = "\t"; no = split("mutect2,strelka2,muse,varscan", ord, ",") }
      FILENAME == ARGV[1] { mnp[$1 FS $2 FS $3 FS $4] = $5; next }
      { k = $1 FS $2 FS $3 FS $4; cs[k] = $5; if (!(k in mnp)) print }
      END {
        for (k in mnp) {
          split(k, f, FS)
          for (o in inter) delete inter[o]; first = 1
          for (i = 1; i <= length(f[3]); i++) {
            rb = substr(f[3], i, 1); ab = substr(f[4], i, 1); if (rb == ab) continue
            for (o in cur) delete cur[o]
            m = split(cs[f[1] FS (f[2] + i - 1) FS rb FS ab], c, ",")
            for (j = 1; j <= m; j++) if (c[j] != "" && c[j] != "mutect2") cur[c[j]] = 1
            if (first) { for (o in cur) inter[o] = 1; first = 0 }
            else       { for (o in inter) if (!(o in cur)) delete inter[o] }
          }
          if (mnp[k] == 1) inter["mutect2"] = 1
          out = ""; nc = 0
          for (x = 1; x <= no; x++) if (ord[x] in inter) { out = (out == "" ? ord[x] : out "," ord[x]); nc++ }
          if (nc > 0) print f[1], f[2], f[3], f[4], out, nc
        }
      }' m2_mnp.tsv agg.tsv | sort -k1,1 -k2,2n -k3,3 -k4,4 > agg.mnp.tsv
    mv agg.mnp.tsv agg.tsv

    # Clamp thresholds to the callers that RAN (input VCF supplied), per variant type.
    # A caller that ran and called nothing still votes "no"; only a disabled caller
    # (input []) is excluded. MuSE calls SNVs only, so it cannot vote on indels.
    n_snv_callers=$((1 - NO_M2)); n_indel_callers=$((1 - NO_M2))   # Mutect2, when it ran
    if [ -s "${STRELKA_SNV}" ] || [ -s "${STRELKA_INDEL}" ]; then
        n_snv_callers=$((n_snv_callers + 1)); n_indel_callers=$((n_indel_callers + 1)); fi
    if [ -s "${MUSE}" ]; then n_snv_callers=$((n_snv_callers + 1)); fi
    if [ -s "${VARSCAN_SNV}" ] || [ -s "${VARSCAN_INDEL}" ]; then
        n_snv_callers=$((n_snv_callers + 1)); n_indel_callers=$((n_indel_callers + 1)); fi
    eff_snv=${MIN_SNV};     [ "$n_snv_callers"   -lt "$eff_snv" ]   && eff_snv=$n_snv_callers
    eff_indel=${MIN_INDEL}; [ "$n_indel_callers" -lt "$eff_indel" ] && eff_indel=$n_indel_callers
    [ "$eff_snv"   -lt 1 ] && eff_snv=1
    [ "$eff_indel" -lt 1 ] && eff_indel=1

    # union only: sites with no Mutect2 PASS vote that still reach the threshold (without
    # Mutect2 and in tag mode, every site, sub-threshold ones as LOWCONSENSUS).
    # Source record = first caller in priority strelka2 > varscan > muse.
    alltag=0; [ "$NO_M2" -eq 1 ] && [ "${FILTER_MODE}" = "tag" ] && alltag=1
    awk -F'\t' -v es="$eff_snv" -v ei="$eff_indel" -v alltag="$alltag" 'BEGIN { OFS = "\t"; np = split("strelka2,varscan,muse", pr, ",") }
      index("," $5 ",", ",mutect2,") { next }
      { t = (length($3) == length($4)) ? es : ei; ok = ($6 + 0 >= t + 0); if (!ok && !alltag) next
        for (i = 1; i <= np; i++) if (index("," $5 ",", "," pr[i] ",")) { print $0, pr[i], (ok ? "PASS" : "LOWCONSENSUS"); break } }' agg.tsv > fill.tsv
    bgzip -f fill.tsv; tabix -f -s1 -b2 -e2 fill.tsv.gz
    bgzip -f agg.tsv;  tabix -f -s1 -b2 -e2 agg.tsv.gz

    cat > cons_hdr.txt <<'HDR'
##INFO=<ID=CALLERS,Number=.,Type=String,Description="Somatic callers supporting this variant">
##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description="Number of somatic callers supporting this variant">
##INFO=<ID=CONS_SOURCE,Number=1,Type=String,Description="Caller whose record was used when Mutect2 has no record at the site (union backbone); absent = Mutect2 record">
##INFO=<ID=M2_FILTER,Number=1,Type=String,Description="Original Mutect2 FILTER ('|'-separated) of a Mutect2-filtered record kept by the caller vote (union backbone)">
HDR
    { cat cons_hdr.txt; echo '##INFO=<ID=CONS_FILTER,Number=1,Type=String,Description="FILTER of a consensus site without a Mutect2 record">'; } > fill_hdr.txt

    parts=()
    if [ "$NO_M2" -eq 0 ]; then
        bcftools annotate -a agg.tsv.gz -h cons_hdr.txt \
            -c CHROM,POS,REF,ALT,INFO/CALLERS,INFO/NCALLERS \
            m2.norm.vcf.gz -O z -o m2.annot.vcf.gz
        bcftools index -f -t m2.annot.vcf.gz

        # union: keep Mutect2-filtered records only when the vote reaches the threshold;
        # set them PASS and keep the Mutect2 FILTER in INFO/M2_FILTER for review.
        if [ "${BACKBONE}" = "union" ]; then
            keep="FILTER=\"PASS\" || FILTER=\".\" || ((TYPE=\"snp\" || TYPE=\"mnp\") && INFO/NCALLERS>=$eff_snv) || (TYPE!=\"snp\" && TYPE!=\"mnp\" && INFO/NCALLERS>=$eff_indel)"
            bcftools view -i "$keep" m2.annot.vcf.gz \
              | awk -F'\t' 'BEGIN { OFS = "\t" } /^#/ { print; next }
                  $7 != "PASS" && $7 != "." { f = $7; gsub(/;/, "|", f); $8 = ($8 == "." ? "" : $8 ";") "M2_FILTER=" f; $7 = "PASS" }
                  { print }' \
              | bcftools view -O z -o m2.bb.vcf.gz -
            bcftools index -f -t m2.bb.vcf.gz
        else
            mv m2.annot.vcf.gz m2.bb.vcf.gz; mv m2.annot.vcf.gz.tbi m2.bb.vcf.gz.tbi
        fi

        # split by type, apply per-type threshold, recombine (all Mutect2 - concat-safe).
        # 'other' (complex / non-SNV-non-indel records) goes with the indel threshold.
        if [ "${FILTER_MODE}" = "tag" ]; then
            bcftools view -v snps,mnps -O u m2.bb.vcf.gz | bcftools filter -m + -s LOWCONSENSUS -e "INFO/NCALLERS < $eff_snv"   -O z -o part.snv.vcf.gz
            bcftools view -v indels,other -O u m2.bb.vcf.gz | bcftools filter -m + -s LOWCONSENSUS -e "INFO/NCALLERS < $eff_indel" -O z -o part.indel.vcf.gz
        else
            bcftools view -v snps,mnps -i "INFO/NCALLERS >= $eff_snv"   -O z -o part.snv.vcf.gz   m2.bb.vcf.gz
            bcftools view -v indels,other -i "INFO/NCALLERS >= $eff_indel" -O z -o part.indel.vcf.gz m2.bb.vcf.gz
        fi
        bcftools index -f -t part.snv.vcf.gz
        bcftools index -f -t part.indel.vcf.gz
        parts=(part.snv.vcf.gz part.indel.vcf.gz)
    fi   # Mutect2 records

    # union: sites Mutect2 never emitted, written in the Mutect2 layout (same header and
    # sample order; FORMAT GT:AD:AF:DP from the source caller) so FILTER_SOMATIC_VCF,
    # VEP and vcf2maf treat them the same. Sites overlapping a KEPT Mutect2 record are
    # excluded (that record already represents the locus, e.g. an MNP or a deletion);
    # a dropped Mutect2 record does not mask other callers' calls.
    n_fill=0
    if [ "$NO_M2" -eq 1 ]; then
        ## the output VCF of its own: normal, tumour; FORMAT as the filled records
        {   echo '##fileformat=VCFv4.2'
            echo '##FILTER=<ID=PASS,Description="All filters passed">'
            echo '##FILTER=<ID=LOWCONSENSUS,Description="Called by fewer somatic callers than the consensus threshold">'
            cat cons_hdr.txt
            echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
            echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref, alt)">'
            echo '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">'
            echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">'
            awk -F'\t' '{ print "##contig=<ID=" $1 ",length=" $2 ">" }' "${FASTA}.fai"
            printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s\n' "${NORMAL}" "${SAMPLE}"
        } > nom2_header.vcf
    fi
    if [ "${BACKBONE}" = "union" ] && { [ "$NO_M2" -eq 1 ] || [ "$(zcat fill.tsv.gz | wc -l)" -gt 0 ]; }; then
        mask=()
        if [ "$NO_M2" -eq 1 ]; then
            m2_t=${SAMPLE}; m2_n=${NORMAL}; m2_first=${NORMAL}
        else
            m2_t=$(bcftools query -l m2.bb.vcf.gz | grep -Fx "${SAMPLE}" || true)
            m2_n=$(bcftools query -l m2.bb.vcf.gz | grep -Fx "${NORMAL}" || true)
            m2_first=$(bcftools query -l m2.bb.vcf.gz | head -1)
            if [ -z "$m2_t" ] || [ -z "$m2_n" ]; then
                echo "ERROR: samples '${SAMPLE}' / '${NORMAL}' not found in ${MUTECT2}" >&2; exit 1
            fi
            # mask = Mutect2 records that SURVIVED the vote (a record dropped for too few callers
            # must not hide another caller's overlapping call)
            bcftools concat -a part.snv.vcf.gz part.indel.vcf.gz -O u | bcftools view -f PASS,. -O u \
              | bcftools query -f '%CHROM\t%POS0\t%END\n' > m2_mask.bed
            [ -s m2_mask.bed ] && mask=(-T ^m2_mask.bed)
        fi
        : > fill.body
        while IFS=$'\t' read -r bcf caller; do
            [ "$caller" = "mutect2" ] && continue
            # a caller with no selected records contributes nothing (and may lack sample columns)
            [ "$(bcftools view -H "$bcf" | awk 'END { print NR }')" -eq 0 ] && continue
            # callers write TUMOR/NORMAL, or the read-group sample names
            c_t=$(bcftools query -l "$bcf" | awk -v t="${SAMPLE}" 'toupper($0) ~ /^TUMOU?R$/ || $0 == t' | head -1 || true)
            c_n=$(bcftools query -l "$bcf" | awk -v n="${NORMAL}" 'toupper($0) == "NORMAL" || $0 == n' | head -1 || true)
            if [ -z "$c_t" ] || [ -z "$c_n" ]; then
                echo "ERROR: could not resolve tumour/normal samples in $bcf ($caller)" >&2; exit 1
            fi
            bcftools view "${mask[@]}" -O u "$bcf" \
              | bcftools annotate -x INFO -O u - \
              | bcftools annotate -a fill.tsv.gz -h fill_hdr.txt \
                    -c CHROM,POS,REF,ALT,INFO/CALLERS,INFO/NCALLERS,INFO/CONS_SOURCE,INFO/CONS_FILTER -O u - \
              | bcftools view -i "INFO/CONS_SOURCE=\"$caller\"" -O b -o fill_sel.bcf -
            bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%INFO/CONS_FILTER\tCALLERS=%INFO/CALLERS;NCALLERS=%INFO/NCALLERS;CONS_SOURCE=%INFO/CONS_SOURCE\tGT:AD:AF:DP\n' fill_sel.bcf > f.fix
            bcftools query -s "$c_t" -f '[./.:%AD:%AF:%DP]\n' fill_sel.bcf > f.t
            bcftools query -s "$c_n" -f '[./.:%AD:%AF:%DP]\n' fill_sel.bcf > f.n
            if [ "$m2_first" = "$m2_t" ]; then paste f.fix f.t f.n; else paste f.fix f.n f.t; fi >> fill.body
        done < sel.list
        n_fill=$(wc -l < fill.body)
        { if [ "$NO_M2" -eq 1 ]; then cat nom2_header.vcf; else bcftools view -h m2.bb.vcf.gz; fi
          sort -k1,1V -k2,2n fill.body; } \
          | bcftools sort -O z -o part.fill.vcf.gz -
        bcftools index -f -t part.fill.vcf.gz
        parts+=(part.fill.vcf.gz)
    fi

    bcftools concat -a "${parts[@]}" -O u \
      | bcftools sort -O z -o ${PREFIX}.vcf.gz -
    bcftools index -f -t ${PREFIX}.vcf.gz

    {   echo "CONSENSUS SUMMARY: ${SAMPLE}"
        if [ $((n_snv_callers > n_indel_callers ? n_snv_callers : n_indel_callers)) -le 1 ]; then
            echo "mode: single caller (somatic_consensus: false) - its somatic PASS calls, no vote; source in INFO/CONS_SOURCE"
        else
            echo "mode: N-of-M vote over all callers; backbone=${BACKBONE}; provenance in INFO/CALLERS,NCALLERS,CONS_SOURCE,M2_FILTER; filter=${FILTER_MODE}"
        fi
        echo "callers_run_snv=$n_snv_callers  callers_run_indel=$n_indel_callers  eff_snv=$eff_snv  eff_indel=$eff_indel  muse_exclude=${MUSE_EXCLUDE}"
        echo "mutect2 PASS input:   $( [ "$NO_M2" -eq 1 ] && echo 'not run' || bcftools view -H -f PASS,. ${MUTECT2} | wc -l)"
        echo "mutect2-filtered kept by vote (M2_FILTER): $(bcftools query -f '%INFO/M2_FILTER\n' ${PREFIX}.vcf.gz | awk '$1 != "." { n++ } END { print n + 0 }')"
        echo "sites without a Mutect2 record (CONS_SOURCE): $n_fill"
        echo "consensus SNVs:       $(bcftools view -H -v snps   ${PREFIX}.vcf.gz | wc -l)"
        echo "consensus MNPs:       $(bcftools view -H -v mnps   ${PREFIX}.vcf.gz | wc -l)"
        echo "consensus indels:     $(bcftools view -H -v indels ${PREFIX}.vcf.gz | wc -l)"
        echo "total consensus:      $(bcftools view -H ${PREFIX}.vcf.gz | wc -l)"
    } > ${PREFIX}_stats.txt
    
