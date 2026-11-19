#!/usr/bin/awk -f
#------------------------------------------------------------------------------
# Author : Pankaj Chauhan
# Purpose: Somatic/germline pre-classification of VEP-annotated variants (CLASSIFY_VARIANTS).
#          Called by the PAPIO process CLASSIFY_VARIANTS; moved from an inline heredoc.
# Usage  : awk -v AF_CUTOFF=.. -v COSMIC_MIN=.. -v COSMIC_STRONG=.. -v MIN_STARS=.. -v MODE=somatic|germline -v BA1=.. -v GERM_LOW_VAF=.. -v MOSAIC_MIN_DP=.. -f classify_variants.awk variants_raw.tsv
#          reads cgc_lookup.tsv and af_thr.tsv from the working directory
#------------------------------------------------------------------------------
BEGIN {
    OFS = "\t"
    while ((getline line < "cgc_lookup.tsv") > 0) {
        split(line, a, "\t"); role[a[1]] = a[2]; tier[a[1]] = a[3]; germ[a[1]] = a[4]
    }
    while ((getline line < "af_thr.tsv") > 0) {
        split(line, a, "\t"); gba1[a[1]] = a[2]; gbs1[a[1]] = a[3]
    }
}
function clinvar_class(s,   t, i, n, tok, hp, hb) {
    if (s == "" || s == ".") return "other"
    s = tolower(s)
    if (s ~ /conflicting/) return "conflicting"
    gsub(/[\/,&|]+/, " ", s)
    n = split(s, t, /[ ]+/)
    for (i = 1; i <= n; i++) {
        tok = t[i]
        if (tok == "pathogenic" || tok == "likely_pathogenic") hp = 1
        if (tok == "benign"     || tok == "likely_benign")     hb = 1
    }
    if (hp && !hb) return "pathogenic"
    if (hb && !hp) return "benign"
    if (hp && hb)  return "conflicting"
    return "other"
}
# ClinVar ONC keeps its grade: Oncogenic -> oncogenic; Likely_oncogenic or the mixed
# "Oncogenic/Likely_oncogenic" -> likely_oncogenic (conservative).
function clinvar_onc_class(s,   t, i, n, ho, hlo, hb) {
    if (s == "" || s == ".") return "other"
    s = tolower(s)
    if (s ~ /conflict/) return "conflicting"
    gsub(/[\/,&|]+/, " ", s)
    n = split(s, t, /[ ]+/)
    for (i = 1; i <= n; i++) {
        if (t[i] == "oncogenic")                            ho  = 1
        if (t[i] == "likely_oncogenic")                     hlo = 1
        if (t[i] == "benign" || t[i] == "likely_benign")    hb  = 1
    }
    if ((ho || hlo) && hb) return "conflicting"
    if (hlo)               return "likely_oncogenic"
    if (ho)                return "oncogenic"
    if (hb)                return "benign"
    return "other"
}
# ClinVar review status -> stars. VEP --custom treats a two-element comma list as REF,ALT
# and keeps only the second element, so "criteria_provided,_single_submitter" arrives as
# "_single_submitter" (and "...,_conflicting_classifications" as "_conflicting_..."): match
# every element on its own, not only "criteria_provided".
function clinvar_stars(s) {
    if (s == "" || s == ".")            return 0
    s = tolower(s)
    if (s ~ /no_assertion/)             return 0
    if (s ~ /no_classification/)        return 0
    if (s ~ /practice_guideline/)       return 4
    if (s ~ /reviewed_by_expert_panel/) return 3
    if (s ~ /multiple_submitters/)      return 2
    if (s ~ /single_submitter/)         return 1
    if (s ~ /conflicting/)              return 1
    if (s ~ /criteria_provided/)        return 1
    return 0
}
{
    chrom=$1; pos=$2; ref=$3; alt=$4; gene=$5; csq=$6; impact=$7
    af=$8; clin=$9; revstat=$10; cosmic=$11; onc_cv=$12; onc_rev=$13; exon=$14
    cdspos=$15; nmd=$16; ad=$17
    in_cgc   = (gene in role)
    cgc_role = in_cgc         ? role[gene] : "NA"
    cgc_tier = (gene in tier) ? tier[gene] : "NA"
    is_common    = (af     != "" && af     != "." && af+0     >  AF_CUTOFF+0)
    is_recurrent = (cosmic != "" && cosmic != "." && cosmic+0 >= COSMIC_MIN+0)
    high         = (impact == "HIGH")
    moderate     = (impact == "MODERATE")
    cv    = clinvar_class(clin)
    stars = clinvar_stars(revstat)
    if ((cv == "pathogenic" || cv == "benign") && stars < MIN_STARS) cv = "other"
    oc       = clinvar_onc_class(onc_cv)
    oc_stars = clinvar_stars(onc_rev)
    if (oc != "other" && oc != "conflicting" && oc_stars < MIN_STARS) oc = "other"
    lof = (csq ~ /(stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost)/)
    tsg = (tolower(cgc_role) ~ /tsg/)
    onc = (tolower(cgc_role) ~ /oncogene/)
    is_hotspot = (cosmic != "" && cosmic != "." && cosmic+0 >= COSMIC_STRONG+0)
    germ_ok    = ((gene in germ) && germ[gene] != "N")    # CGC GERMLINE=yes, or no GERMLINE column
    # NMD escape (ClinGen SVI PVS1; SOP OVS1 caveat): a nonsense/frameshift in the last exon, or
    # in the last 50 nt of the penultimate exon, of a multi-exon transcript. The Ensembl NMD plugin
    # also flags the first 100 coding bases and single-exon transcripts, which the SVI tree does not
    # downgrade, so its flag is used only for the penultimate-exon case.
    ne = split(exon, ex, "/"); split(cdspos, cp, "-")
    if (ne == 2) { nr = split(ex[1], er, "-"); ex[1] = er[nr] }     # "14-15/16" -> exon 15
    trunc = (csq ~ /(stop_gained|frameshift_variant)/)
    multi = (ne == 2 && ex[2]+0 > 1)
    pen50 = (ex[1]+0 == ex[2]-1 && nmd == "NMD_escaping_variant" && cp[1] ~ /^[0-9]+$/ && cp[1]+0 > 101)
    nmd_escape = trunc && multi && (ex[1]+0 == ex[2]+0 || pen50)
    # Gene-specific germline frequency thresholds (VCEP), else the global ones.
    ba1g = (gene in gba1) ? gba1[gene] : BA1
    bs1g = (gene in gbs1) ? gbs1[gene] : AF_CUTOFF
    # VAF from the sample AD (ref,alt after allele splitting).
    vaf = "."; depth = 0
    if (split(ad, d, ",") == 2 && d[1] ~ /^[0-9]+$/ && d[2] ~ /^[0-9]+$/ && d[1] + d[2] > 0) {
        depth = d[1] + d[2]; vaf = sprintf("%.3f", d[2] / depth) }
    if (MODE == "germline") {
        # ACMG-style vocabulary. ClinVar (star-gated) is a GERMLINE pathogenicity source.
        if      (cv == "pathogenic")                cls = "PATHOGENIC_OR_LIKELY"
        else if (cv == "benign")                    cls = "BENIGN_OR_LIKELY"
        else if (af != "" && af != "." && af+0 > ba1g+0)  cls = "BENIGN_OR_LIKELY"        # BA1
        else if (af != "" && af != "." && af+0 > bs1g+0)  cls = "LIKELY_BENIGN"           # BS1
        else if (lof && tsg && germ_ok)             cls = "LOF_CANCER_GENE_REVIEW"   # LoF mechanism: germline TSG
        else if (high || moderate)                  cls = "VUS"
        else                                        cls = "LOW_IMPACT"        # not assessed; not a benign call
    } else {
        # Somatic ONCOGENICITY, gene-role aware; points follow the ClinGen/CGC/VICC SOP
        # (Horak 2022: >=10 oncogenic, 6-9 likely oncogenic):
        #   TSG null variant (OVS1, 8)                -> MODERATE; + COSMIC >= cosmic_recurrence_min (OM3) -> STRONG
        #     last-exon truncation (NMD escape)       -> one step lower (VUS_CANCER_GENE; MODERATE if recurrent)
        #   no coding/splice impact (LOW/MODIFIER)    -> LOW_IMPACT (not assessed; not LIKELY_BENIGN)
        #   oncogene non-null, COSMIC-recurrent       -> MODERATE; >= cosmic_strong_min (OS3-like) -> STRONG
        #   null variant in an oncogene               -> VUS_CANCER_GENE (LoF is not the mechanism; review)
        #   ClinVar ONC (star-gated) sets a floor     : Oncogenic -> STRONG, Likely oncogenic -> >= MODERATE
        #   ClinVar germline Benign/LB (star-gated)   -> LIKELY_BENIGN unless ClinVar ONC is (likely) oncogenic
        #   ClinVar germline P/LP                     -> shown only, never oncogenicity evidence
        if (oc == "benign" || is_common)            cls = "LIKELY_BENIGN"
        else {
            if      (tsg && lof && nmd_escape)      cls = is_recurrent ? "ONCOGENIC_MODERATE" : "VUS_CANCER_GENE"
            else if (tsg && lof)                    cls = is_recurrent ? "ONCOGENIC_STRONG" : "ONCOGENIC_MODERATE"
            else if (onc && !lof && moderate && is_recurrent) cls = is_hotspot ? "ONCOGENIC_STRONG" : "ONCOGENIC_MODERATE"
            else if (onc && lof)                    cls = "VUS_CANCER_GENE"
            else if (in_cgc && is_recurrent && (high || moderate)) cls = "ONCOGENIC_MODERATE"
            else if (in_cgc && (high || moderate))  cls = "VUS_CANCER_GENE"
            else if (high || moderate)              cls = "VUS"
            else                                    cls = "LOW_IMPACT"    # not assessed; not a benign call
            if      (oc == "oncogenic")             cls = "ONCOGENIC_STRONG"
            else if (oc == "likely_oncogenic")      { if (cls != "ONCOGENIC_STRONG") cls = "ONCOGENIC_MODERATE" }
            else if (cv == "benign")                cls = "LIKELY_BENIGN"
        }
    }
    # Origin flag (does not change the class). Germline: a VAF well below heterozygous suggests
    # mosaicism or clonal haematopoiesis (blood DNA). Paired somatic: absent from the matched normal.
    if (MODE == "germline")
        origin = (vaf == "." || depth < MOSAIC_MIN_DP+0) ? "UNRESOLVED" : (vaf+0 < GERM_LOW_VAF+0 ? "POSSIBLE_MOSAIC_OR_CH" : "GERMLINE")
    else
        origin = "SOMATIC"
    print chrom, pos, ref, alt, gene, csq, impact, (af==""?".":af), (clin==""?".":clin), stars, (cosmic==""?".":cosmic), cgc_role, cgc_tier, cls, (onc_cv==""?".":onc_cv), oc_stars, vaf, origin
}
