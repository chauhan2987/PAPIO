#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Shared JVM-options helper (single source of truth in utils).
include { gatkJavaOpts; pubDir } from './utils'

// ============================================================
//  components.nf  -  All Process Definitions
// ============================================================



// ============================================================
//  STAGE 5A: SV - Manta somatic (tumor-normal)
// ============================================================
process MANTA {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    // scratch = true set globally in nextflow.config

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/sv/manta" }, mode: 'copy'

    input:
    tuple val(meta),
          path(tumor_bam), path(tumor_bai),
          path(normal_bam), path(normal_bai)
    path  ref_fasta
    path  ref_fai
    path  call_regions_gz
    path  call_regions_tbi

    output:
    tuple val(meta),
          path("${meta.sample_id}_manta_somaticSV.vcf.gz"),
          path("${meta.sample_id}_manta_somaticSV.vcf.gz.tbi"),           emit: sv_vcf
    tuple val(meta),
          path("${meta.sample_id}_manta_candidateSmallIndels.vcf.gz"),
          path("${meta.sample_id}_manta_candidateSmallIndels.vcf.gz.tbi"), emit: indel_candidates

    script:
    def call_regions = (meta.assay == 'wes' && call_regions_gz) ? "--callRegions ${call_regions_gz}" : ''
    def exome_flag   = (meta.assay == 'wes') ? '--exome' : ''

    """
    set -euo pipefail

    # Ensure the BAM indexes match Manta's expected string matching patterns (*.bam.bai)
    if [ ! -f "${tumor_bam}.bai" ]; then ln -s ${tumor_bai} "${tumor_bam}.bai"; fi
    if [ ! -f "${normal_bam}.bai" ]; then ln -s ${normal_bai} "${normal_bam}.bai"; fi

    configManta.py \\
        --tumorBam       ${tumor_bam} \\
        --normalBam      ${normal_bam} \\
        --referenceFasta ${ref_fasta} \\
        ${call_regions} ${exome_flag} \\
        --runDir manta_run

    manta_run/runWorkflow.py -m local -j ${task.cpus}

    cp manta_run/results/variants/somaticSV.vcf.gz               ${meta.sample_id}_manta_somaticSV.vcf.gz
    cp manta_run/results/variants/somaticSV.vcf.gz.tbi            ${meta.sample_id}_manta_somaticSV.vcf.gz.tbi
    cp manta_run/results/variants/candidateSmallIndels.vcf.gz      ${meta.sample_id}_manta_candidateSmallIndels.vcf.gz
    cp manta_run/results/variants/candidateSmallIndels.vcf.gz.tbi  ${meta.sample_id}_manta_candidateSmallIndels.vcf.gz.tbi
    """
}


// ============================================================
//  STAGE 5B-SNV/INDEL: SOMATIC - Strelka2
// ============================================================
process STRELKA2 {
    tag { "${meta.sample_id}" }
    label 'process_strelka'

    // scratch = true set globally in nextflow.config

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/strelka2" }, mode: 'copy'

    maxRetries    2

    input:
    tuple val(meta),
          path(tumour_bam), path(tumour_bai),
          path(normal_bam), path(normal_bai)
    path  ref_fasta
    path  ref_fai
    path  call_regions_gz
    path  call_regions_tbi
    path  manta_indel_candidates
    path  manta_indel_candidates_tbi

    output:
    tuple val(meta),
          path("${meta.sample_id}_strelka2_snvs.vcf.gz"),
          path("${meta.sample_id}_strelka2_snvs.vcf.gz.tbi"),
          path("${meta.sample_id}_strelka2_indels.vcf.gz"),
          path("${meta.sample_id}_strelka2_indels.vcf.gz.tbi"), emit: vcf

    script:
    def call_regions = (meta.assay == 'wes' && call_regions_gz) ? "--callRegions ${call_regions_gz}" : ''
    def indel_cands  = manta_indel_candidates ? "--indelCandidates ${manta_indel_candidates}" : ''
    def exome_flag   = (meta.assay == 'wes') ? '--exome' : ''
    def mem_gb       = task.memory ? task.memory.toGiga().intValue() : 60
    def jobs         = task.cpus ? (task.cpus as Integer) : 4
    """
    set -euo pipefail
    tmp_dir="\$(pwd)/strelka_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    # Strelka (Illumina workflow) expects <bam>.bai; our index is <sample>_markdup.bai.
    if [ ! -f "${tumour_bam}.bai" ]; then ln -s ${tumour_bai} "${tumour_bam}.bai"; fi
    if [ ! -f "${normal_bam}.bai" ]; then ln -s ${normal_bai} "${normal_bam}.bai"; fi
    
    configureStrelkaSomaticWorkflow.py \\
        --tumorBam       ${tumour_bam} \\
        --normalBam      ${normal_bam} \\
        --referenceFasta ${ref_fasta} \\
        ${call_regions} ${indel_cands} ${exome_flag} \\
        --runDir strelka_run

    # --memGb (-g): explicit memory budget - matches production pipelines
    #   (UCSF Francis Lab, Janis). Without it, Strelka2 auto-detects available
    #   memory and oversubscribes its -j workers on deep (250-500X) data,
    #   causing silent SIGSEGV (-11) on the heaviest genome segments.
    strelka_run/runWorkflow.py -m local -j ${jobs} -g ${mem_gb}

    cp strelka_run/results/variants/somatic.snvs.vcf.gz       ${meta.sample_id}_strelka2_snvs.vcf.gz
    cp strelka_run/results/variants/somatic.snvs.vcf.gz.tbi   ${meta.sample_id}_strelka2_snvs.vcf.gz.tbi
    cp strelka_run/results/variants/somatic.indels.vcf.gz      ${meta.sample_id}_strelka2_indels.vcf.gz
    cp strelka_run/results/variants/somatic.indels.vcf.gz.tbi  ${meta.sample_id}_strelka2_indels.vcf.gz.tbi
    """
}


// ============================================================
//  STAGE 5C Add STRELKA AF: SOMATIC - Format AF for strelka output
// ============================================================

// Inject standard AF/AD FORMAT fields into Strelka2 somatic VCFs.
// Strelka does not emit AF/AD natively. Per Strelka's user guide:
//   SNVs:   tier1 ref/alt counts come from FORMAT/{A,C,G,T}U (field matching
//           REF/ALT base), tier1 = first comma-delimited value.
//   Indels: tier1 ref/alt counts come from FORMAT/TAR (ref) and FORMAT/TIR (alt).
//   AF = alt_tier1 / (alt_tier1 + ref_tier1)

process TAG_STRELKA_AF {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/strelka2" }, mode: 'copy'

    input:
    tuple val(meta),
          path(strelka_snv),   path(strelka_snv_tbi),
          path(strelka_indel), path(strelka_indel_tbi)

    output:
    tuple val(meta),
          path("${meta.sample_id}_strelka2_snvs.af.vcf.gz"),
          path("${meta.sample_id}_strelka2_snvs.af.vcf.gz.tbi"),
          path("${meta.sample_id}_strelka2_indels.af.vcf.gz"),
          path("${meta.sample_id}_strelka2_indels.af.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    cat > strelka_add_af.py <<'PYEOF'
#!/usr/bin/env python3

import sys
import pysam

def main():
    vcf_in, vcf_out, variant_type = sys.argv[1], sys.argv[2], sys.argv[3]
    assert variant_type in ("snv", "indel"), "variant_type must be 'snv' or 'indel'"

    vcf = pysam.VariantFile(vcf_in)
    vcf.header.formats.add("AD", ".", "Integer",
        "Allelic depths (ref,alt); tier1 counts derived from Strelka U/TAR-TIR fields")
    vcf.header.formats.add("AF", "1", "Float",
        "Allele frequency; tier1 alt/(tier1 alt+tier1 ref), derived per Strelka user guide")

    out = pysam.VariantFile(vcf_out, "w", header=vcf.header)
    base_field = {"A": "AU", "C": "CU", "G": "GU", "T": "TU"}

    n_ok, n_fallback = 0, 0
    for rec in vcf:
        ref = rec.ref.upper()
        alt = rec.alts[0].upper() if rec.alts else None

        for sample_name in rec.samples:
            s = rec.samples[sample_name]
            ref_t1, alt_t1 = 0, 0
            try:
                if variant_type == "snv":
                    ref_field = base_field.get(ref)
                    alt_field = base_field.get(alt) if alt else None
                    if ref_field and ref_field in s and s[ref_field] is not None:
                        ref_t1 = s[ref_field][0]
                    if alt_field and alt_field in s and s[alt_field] is not None:
                        alt_t1 = s[alt_field][0]
                else:
                    if "TAR" in s and s["TAR"] is not None:
                        ref_t1 = s["TAR"][0]
                    if "TIR" in s and s["TIR"] is not None:
                        alt_t1 = s["TIR"][0]
                n_ok += 1
            except Exception:
                n_fallback += 1

            denom = ref_t1 + alt_t1
            af = round(alt_t1 / denom, 4) if denom > 0 else 0.0
            s["AD"] = (int(ref_t1), int(alt_t1))
            s["AF"] = af

        out.write(rec)

    out.close()
    sys.stderr.write(f"strelka_add_af: {n_ok} records tagged OK, {n_fallback} fallback-to-zero\\n")

if __name__ == "__main__":
    main()
PYEOF

    python3 strelka_add_af.py ${strelka_snv}   snvs.af.vcf   snv
    python3 strelka_add_af.py ${strelka_indel} indels.af.vcf indel

    bgzip -c snvs.af.vcf   > ${meta.sample_id}_strelka2_snvs.af.vcf.gz
    bgzip -c indels.af.vcf > ${meta.sample_id}_strelka2_indels.af.vcf.gz
    bcftools index -t ${meta.sample_id}_strelka2_snvs.af.vcf.gz
    bcftools index -t ${meta.sample_id}_strelka2_indels.af.vcf.gz
    """
}

// ============================================================
//  STAGE 5D-SNV/INDEL: SOMATIC - MuSE call (per tumor/normal pair)
// ============================================================
process MUSE_CALL {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    maxRetries    2

    input:
    tuple val(meta),
        path(tumor_bam), path(tumor_bai),
        path(normal_bam), path(normal_bai)
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta), path("${meta.sample_id}.MuSE.txt"), emit: call_output

    script:
    def java_opts = gatkJavaOpts(task)
    """
        set -euo pipefail
        tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_muse_tmp"
        mkdir -p "\$tmp_dir"
        trap 'rm -rf "\$tmp_dir"' EXIT
        # MuSE call fits a Markov substitution model per tumor/normal pair and
        # writes raw candidate SNV calls with their computed tier scores to a
        # flat intermediate file; MUSE_SUMP applies the population-level cutoff
        # against dbSNP in a separate, much lighter step.
        MuSE call \\
            -f ${ref_fasta} \\
            -O ${meta.sample_id} \\
            -n ${task.cpus} \\
            ${tumor_bam} ${normal_bam}
        """
}


// ============================================================
//  STAGE 5E-SNV/INDEL: SOMATIC - MuSE sump (per tumor/normal pair)
// ============================================================

process MUSE_SUMP {
    tag { "${meta.sample_id}" }
    label 'process_low'
    maxRetries    2

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/muse" }, mode: 'copy'

    input:
    tuple val(meta), path(muse_call_txt)
        path  dbsnp_vcf
        path  dbsnp_vcf_tbi

    output:
    tuple val(meta),
          path("${meta.sample_id}.MuSE.vcf.gz"),
          path("${meta.sample_id}.MuSE.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    def mode_flag = (meta.assay == 'wes') ? '-E' : '-G'

    """
    set -euo pipefail
    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_muse_sump_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT
    MuSE sump \\
        -I ${muse_call_txt} \\
        ${mode_flag} \\
        -D ${dbsnp_vcf} \\
        -O ${meta.sample_id}.MuSE.vcf \\
        -n ${task.cpus}
    bgzip ${meta.sample_id}.MuSE.vcf
    bcftools index -t ${meta.sample_id}.MuSE.vcf.gz
    """
}

// ============================================================
//  STAGE 5F Add AF MUSE: SOMATIC - Format AF for MUSE output
// ============================================================
process TAG_MUSE_AF {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/muse" }, mode: 'copy'

    input:
    tuple val(meta), path(muse_vcf), path(muse_tbi)
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta),
          path("${meta.sample_id}_muse.af.vcf.gz"),
          path("${meta.sample_id}_muse.af.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    cat > muse_add_af.py <<'PYEOF'
#!/usr/bin/env python3
# MuSE emits AD as (ref,alt) but no decimal AF. Compute AF = alt/(ref+alt).
import sys, pysam
vcf = pysam.VariantFile(sys.argv[1])
if "AF" not in vcf.header.formats:
    vcf.header.formats.add("AF", "1", "Float", "Allele frequency; computed from MuSE AD (alt/(ref+alt))")
out = pysam.VariantFile(sys.argv[2], "w", header=vcf.header)
for rec in vcf:
    for sn in rec.samples:
        s = rec.samples[sn]
        af = 0.0
        try:
            if "AD" in s and s["AD"] is not None:
                ad = s["AD"]
                if isinstance(ad, (tuple, list)) and len(ad) >= 2:
                    ref_d, alt_d = ad[0] or 0, ad[1] or 0
                    if (ref_d + alt_d) > 0:
                        af = round(alt_d / (ref_d + alt_d), 4)
        except Exception:
            pass
        s["AF"] = af
    out.write(rec)
out.close()
PYEOF

    bcftools norm -m -any -f ${ref_fasta} -O v ${muse_vcf} > split.vcf
    python3 muse_add_af.py split.vcf muse.af.vcf
    bgzip -c muse.af.vcf > ${meta.sample_id}_muse.af.vcf.gz
    bcftools index -t ${meta.sample_id}_muse.af.vcf.gz
    """
}



// ============================================================
//  STAGE 5G-SNV/INDEL: SOMATIC - VARSCAN
// ============================================================
process VARSCAN {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    maxRetries    2

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/varscan" }, mode: 'copy'

    input:
    tuple val(meta),
          path(tumor_bam), path(tumor_bai),
          path(normal_bam), path(normal_bai)
    path  ref_fasta
    path  ref_fai
    path  interval_bed

    output:
    tuple val(meta),
          path("${meta.sample_id}_varscan_snv.vcf.gz"),
          path("${meta.sample_id}_varscan_snv.vcf.gz.tbi"),
          path("${meta.sample_id}_varscan_indel.vcf.gz"),
          path("${meta.sample_id}_varscan_indel.vcf.gz.tbi"), emit: vcf

    script:
    def min_cov     = params.varscan_min_coverage    ?: 8
    def min_var_frq = params.varscan_min_var_freq    ?: 0.05
    def somatic_p   = params.varscan_somatic_pvalue  ?: 0.05
    def region_flag = (meta.assay == 'wes' && interval_bed) ? "-l ${interval_bed}" : ''
    """
    set -euo pipefail

    # samtools mpileup expects <bam>.bai; our index is <sample>_markdup.bai.
    if [ ! -f "${tumor_bam}.bai" ];  then ln -s ${tumor_bai}  "${tumor_bam}.bai";  fi
    if [ ! -f "${normal_bam}.bai" ]; then ln -s ${normal_bai} "${normal_bam}.bai"; fi

    # VarScan somatic reads a COMBINED normal+tumor mpileup on stdin (--mpileup 1).
    # Restrict to capture BED for WES (${region_flag}); q20/Q20 baseline quality.
    samtools mpileup \\
        -f ${ref_fasta} \\
        ${region_flag} \\
        -q 20 -Q 20 \\
        ${normal_bam} ${tumor_bam} \\
    | varscan somatic \\
        /dev/stdin \\
        ${meta.sample_id}_varscan \\
        --mpileup 1 \\
        --min-coverage ${min_cov} \\
        --min-var-freq ${min_var_frq} \\
        --somatic-p-value ${somatic_p} \\
        --output-vcf 1

    # VarScan writes <prefix>.snp.vcf and <prefix>.indel.vcf. Compress+index+rename.
    for kind in snp indel; do
        out=\$( [ "\$kind" = "snp" ] && echo snv || echo indel )
        if [ -s "${meta.sample_id}_varscan.\${kind}.vcf" ]; then
            bgzip -c "${meta.sample_id}_varscan.\${kind}.vcf" > "${meta.sample_id}_varscan_\${out}.vcf.gz"
        else
            # emit a valid empty VCF so downstream staging never fails
            bcftools view -h "${meta.sample_id}_varscan.\${kind}.vcf" 2>/dev/null \\
                | bgzip > "${meta.sample_id}_varscan_\${out}.vcf.gz" || \\
            printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' \\
                | bgzip > "${meta.sample_id}_varscan_\${out}.vcf.gz"
        fi
        bcftools index -t "${meta.sample_id}_varscan_\${out}.vcf.gz"
    done
    """
}

// ============================================================
//  STAGE 5H Add AF VARSCAN: SOMATIC - Format AF for VARSCAN output
// ============================================================
process TAG_VARSCAN_AF {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/varscan" }, mode: 'copy'

    input:
    tuple val(meta),
          path(varscan_snv),   path(varscan_snv_tbi),
          path(varscan_indel), path(varscan_indel_tbi)
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta),
          path("${meta.sample_id}_varscan_snv.af.vcf.gz"),
          path("${meta.sample_id}_varscan_snv.af.vcf.gz.tbi"),
          path("${meta.sample_id}_varscan_indel.af.vcf.gz"),
          path("${meta.sample_id}_varscan_indel.af.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    cat > varscan_add_af.py <<'PYEOF'
#!/usr/bin/env python3
# VarScan: RD=ref depth, AD=alt depth (single value), FREQ="NN.NN%".
# Header AD is redeclared to Number=R by sed upstream so pysam accepts a
# (ref,alt) pair. We only WRITE values here; every read defaults None -> 0
# so no None ever reaches int().
import sys, pysam

def first_num(x):
    # unwrap tuple/list, coerce None -> 0
    if x is None:
        return 0
    if isinstance(x, (tuple, list)):
        x = next((v for v in x if v is not None), 0)
    return x if x is not None else 0

def main():
    vcf = pysam.VariantFile(sys.argv[1])
    if "AF" not in vcf.header.formats:
        vcf.header.formats.add("AF", "1", "Float",
            "Allele frequency (0-1); from VarScan FREQ/AD")
    out = pysam.VariantFile(sys.argv[2], "w", header=vcf.header)
    for rec in vcf:
        for sn in rec.samples:
            s = rec.samples[sn]
            ref_d, alt_d, af = 0, 0, 0.0
            try:
                ref_d = first_num(s.get("RD"))
                # VarScan AD = ALT count; after Number=R reinterpretation it may
                # be a tuple with None padding - take last non-None value.
                ad_raw = s.get("AD")
                if isinstance(ad_raw, (tuple, list)):
                    alt_d = next((v for v in reversed(ad_raw) if v is not None), 0)
                else:
                    alt_d = ad_raw if ad_raw is not None else 0

                fq = s.get("FREQ")
                if fq is not None:
                    fq = fq[0] if isinstance(fq, (tuple, list)) else fq
                    try:
                        af = round(float(str(fq).rstrip("%")) / 100.0, 4)
                    except (ValueError, TypeError):
                        af = 0.0
                elif (ref_d + alt_d) > 0:
                    af = round(alt_d / (ref_d + alt_d), 4)
            except Exception:
                pass
            s["AD"] = (int(ref_d), int(alt_d))
            s["AF"] = af
        out.write(rec)
    out.close()

if __name__ == "__main__":
    main()
PYEOF

    for kind in snv indel; do
        invcf=\$([ "\$kind" = "snv" ] && echo ${varscan_snv} || echo ${varscan_indel})
        # 1. split multiallelics (VarScan single-value RD/AD can't represent multi-alt)
        bcftools norm -m -any -f ${ref_fasta} -O v "\$invcf" > \${kind}.split.vcf
        # 2. fix AD header declaration to Number=R via sed (pysam can't re-add an existing ID)
        sed 's/##FORMAT=<ID=AD,Number=[0-9A-Za-z.]*,Type=Integer/##FORMAT=<ID=AD,Number=R,Type=Integer/' \\
            \${kind}.split.vcf > \${kind}.hdrfix.vcf
        # 3. inject AF + rewrite AD as (ref,alt)
        python3 varscan_add_af.py \${kind}.hdrfix.vcf \${kind}.af.vcf
        bgzip -c \${kind}.af.vcf > ${meta.sample_id}_varscan_\${kind}.af.vcf.gz
        bcftools index -t ${meta.sample_id}_varscan_\${kind}.af.vcf.gz
    done
    """
}

// ============================================================
//  STAGE 5J-SNV/INDEL: SOMATIC - Mutect2 (per interval)
// ============================================================
process MUTECT2 {
    tag { "${meta.sample_id}:${interval_name}" }
    label 'process_variant'

    maxRetries    2

    input:
    tuple val(meta),
          path(tumor_bam), path(tumor_bai),
          path(normal_bam), path(normal_bai),
          val(interval_name), path(interval_file)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  germline_resource
    path  germline_resource_tbi
    path  panel_of_normals
    path  panel_of_normals_tbi

    output:
    tuple val(meta), val(interval_name),
          path("${meta.sample_id}.${interval_name}.vcf.gz"),
          path("${meta.sample_id}.${interval_name}.vcf.gz.tbi"),
          path("${meta.sample_id}.${interval_name}.vcf.gz.stats"), emit: vcf_chunk
    tuple val(meta), val(interval_name),
          path("${meta.sample_id}.${interval_name}.f1r2.tar.gz"), emit: f1r2

    script:
    def java_opts = gatkJavaOpts(task)
    def pon_flag    = panel_of_normals ? "--panel-of-normals ${panel_of_normals}" : ''
    def normal_name = "${meta.patient_id}_N"
    def gr_flag = germline_resource ? "--germline-resource ${germline_resource}" : ''

    """
    set -euo pipefail

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_mutect2_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    # --f1r2-tar-gz collects read-orientation counts per interval; these are
    # gathered per-sample by LearnReadOrientationModel to model FFPE C>T/G>A
    # deamination artifacts (orientation bias).
    gatk --java-options "${java_opts}" Mutect2 \\
        -R ${ref_fasta} \\
        -I ${tumor_bam} -I ${normal_bam} \\
        -normal ${normal_name} \\
        ${gr_flag} \\
        ${pon_flag} \\
        -L ${interval_file} \\
        --f1r2-tar-gz ${meta.sample_id}.${interval_name}.f1r2.tar.gz \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}.${interval_name}.vcf.gz
    """
}

// Use later
// input:
//     tuple val(meta), path(tumor_bam), path(tumor_bai), val(interval_name), path(interval)
//     path(normal_bam)          // real BAM = paired; [] = tumor-only (truthiness, like pon_flag)
//     path(normal_bai)
//     path  ref_fasta
//     path  ref_fai
//     path  ref_dict
//     path  panel_of_normals
//     path  panel_of_normals_tbi
//     path  germline_resource
//     path  germline_resource_tbi
//     // ... outputs unchanged ...
//     script:
//     def java_opts = gatkJavaOpts(task)
//     def pon_flag  = panel_of_normals ? "--panel-of-normals ${panel_of_normals}" : ''
//     def gr_flag   = germline_resource ? "--germline-resource ${germline_resource}" : ''
//     def normal_in = normal_bam ? "-I ${normal_bam} -normal ${meta.patient_id}_N" : ''   // <- the whole change

//     """
//     set -euo pipefail
//     gatk --java-options "${java_opts}" Mutect2 \\
//         -R ${ref_fasta} \\
//         -I ${tumor_bam} ${normal_in} \\
//         ${pon_flag} \\
//         ${gr_flag} \\
//         -L ${interval} \\
//         --f1r2-tar-gz ${meta.sample_id}.${interval_name}.f1r2.tar.gz \\
//         -O ${meta.sample_id}.${interval_name}.unfiltered.vcf.gz
//     """

// ============================================================
//  STAGE 5K-GATHER: Merge + Filter Mutect2 interval chunks
// ============================================================
process MERGE_MUTECT2 {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/mutect2" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf_chunks), path(tbi_chunks), path(stats_chunks)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  contamination_table
    path  tumor_segmentation
    path  ob_priors

    output:
    tuple val(meta),
          path("${meta.sample_id}_mutect2.vcf.gz"),
          path("${meta.sample_id}_mutect2.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    def vcf_inputs   = vcf_chunks.collect   { v -> "-I ${v}" }.join(' \\\n        ')
    def stats_inputs = stats_chunks.collect { v -> "-stats ${v}" }.join(' \\\n        ')
    def contam_flag  = contamination_table ? "--contamination-table ${contamination_table}" : ''
    def seg_flag     = tumor_segmentation  ? "--tumor-segmentation ${tumor_segmentation}" : ''
    def ob_flag      = ob_priors           ? "--orientation-bias-artifact-priors ${ob_priors}" : ''

    """
    set -euo pipefail

    gatk --java-options "${java_opts}" MergeVcfs \\
        ${vcf_inputs} -O ${meta.sample_id}_merged.vcf.gz

    gatk --java-options "${java_opts}" MergeMutectStats \\
        ${stats_inputs} -O ${meta.sample_id}_merged.stats

    gatk --java-options "${java_opts}" FilterMutectCalls \\
        -R ${ref_fasta} \\
        -V ${meta.sample_id}_merged.vcf.gz \\
        --stats ${meta.sample_id}_merged.stats \\
        ${contam_flag} \\
        ${seg_flag} \\
        ${ob_flag} \\
        -O ${meta.sample_id}_mutect2.vcf.gz
    """
}


// ============================================================
//  STAGE 5O: Variant-type-aware ensemble consensus (N-of-M callers)
// ============================================================

//     consensus() {
//         local out="\$1"; shift
//         local minreq="\$1"; shift
//         local files=("\$@")
//         local n=\${#files[@]}
//         if [ "\$n" -eq 0 ]; then
//             bcftools view -h ${mutect2_vcf} -O z -o "\$out"
//             bcftools index -t "\$out"
//             return
//         fi
//         local eff=\$minreq
//         [ "\$n" -lt "\$minreq" ] && eff=\$n
//         if [ "\$n" -eq 1 ]; then
//             cp "\${files[0]}" "\$out"
//         elif [ "\$eff" -le 1 ]; then
//             bcftools concat -a -D "\${files[@]}" -O u | bcftools sort -O z -o "\$out" -
//         else
//             local d=isec_tmp_\$(basename "\$out" .vcf.gz)
//             rm -rf "\$d"
//             bcftools isec -n +\${eff} -c none -w1 -p "\$d" -O z "\${files[@]}"
//             cp "\$d"/0000.vcf.gz "\$out"
//         fi
//         bcftools index -t "\$out"
//     }

//     snv_files=()
//     indel_files=()

//     if [ -s "${mutect2_vcf}" ]; then
//         prep "${mutect2_vcf}" mutect2 pass
//         snv_files+=(mutect2.snv.vcf.gz)
//         indel_files+=(mutect2.indel.vcf.gz)
//     fi
//     if [ -s "${strelka_snv}" ]; then
//         prep "${strelka_snv}" strelka_s pass
//         snv_files+=(strelka_s.snv.vcf.gz)
//     fi
//     if [ -s "${strelka_indel}" ]; then
//         prep "${strelka_indel}" strelka_i pass
//         indel_files+=(strelka_i.indel.vcf.gz)
//     fi
//     if [ -s "${muse_vcf}" ]; then
//         prep "${muse_vcf}" muse muse
//         snv_files+=(muse.snv.vcf.gz)
//     fi
//     if [ -s "${varscan_snv}" ]; then
//         prep "${varscan_snv}" varscan_s pass
//         snv_files+=(varscan_s.snv.vcf.gz)
//     fi
//     if [ -s "${varscan_indel}" ]; then
//         prep "${varscan_indel}" varscan_i pass
//         indel_files+=(varscan_i.indel.vcf.gz)
//     fi

//     consensus consensus_snv.vcf.gz   ${min_snv}   \${snv_files[@]+"\${snv_files[@]}"}
//     consensus consensus_indel.vcf.gz ${min_indel} \${indel_files[@]+"\${indel_files[@]}"}

//     bcftools concat -a -D consensus_snv.vcf.gz consensus_indel.vcf.gz -O u \
//       | bcftools sort -O z -o ${meta.sample_id}_consensus.vcf.gz -
//     bcftools index -t ${meta.sample_id}_consensus.vcf.gz

//     {   echo "CONSENSUS SUMMARY: ${meta.sample_id}"
//         echo "min_snv_callers=${min_snv} min_indel_callers=${min_indel} muse_exclude=${muse_exclude}"
//         echo "snv-capable callers present: \${#snv_files[@]}"
//         echo "indel-capable callers present: \${#indel_files[@]}"
//         echo "consensus SNVs:   \$(bcftools view -H consensus_snv.vcf.gz   | wc -l)"
//         echo "consensus indels: \$(bcftools view -H consensus_indel.vcf.gz | wc -l)"
//         echo "total consensus:  \$(bcftools view -H ${meta.sample_id}_consensus.vcf.gz | wc -l)"
//     } > ${meta.sample_id}_consensus_stats.txt
//     """
// }


// ============================================================
process ENSEMBLE_CONSENSUS {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-aggregated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/callers_consensus" }, mode: 'copy'

    input:
    tuple val(meta),
          path(mutect2_vcf),    path(mutect2_tbi),
          path(strelka_snv),    path(strelka_snv_tbi),
          path(strelka_indel),  path(strelka_indel_tbi),
          path(muse_vcf),       path(muse_tbi),
          path(varscan_snv),    path(varscan_snv_tbi),
          path(varscan_indel),  path(varscan_indel_tbi)
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta),
          path("${meta.sample_id}_consensus.vcf.gz"),
          path("${meta.sample_id}_consensus.vcf.gz.tbi"), emit: vcf
    path  "${meta.sample_id}_consensus_stats.txt",        emit: stats

    script:
    def min_snv      = params.consensus_min_snv_callers   ?: 2
    def min_indel    = params.consensus_min_indel_callers ?: 2
    def muse_exclude = params.muse_exclude_tier           ?: 'Tier5'
    def filter_mode  = params.consensus_filter_mode       ?: 'remove'   // 'remove' | 'tag'
    """
    set -euo pipefail

    # Emit per-variant caller support (CHROM POS REF ALT CALLER) for one caller VCF.
    # PASS/somatic-filtered and left-aligned+split so alleles match across callers.
    emit_support() {
        local invcf="\$1" caller="\$2" mode="\$3" sel
        case "\$mode" in
          muse)    sel=(bcftools view -e "FILTER==\\"${muse_exclude}\\"" -O u "\$invcf") ;;
          varscan)
              # capture the header first: 'grep -q' closing the pipe early under
              # pipefail (large contig header) would make this test spuriously false.
              vs_hdr=\$(bcftools view -h "\$invcf")
              if grep -q 'ID=SS,' <<< "\$vs_hdr"; then
                  sel=(bcftools view -i 'INFO/SS==2' -f PASS,. -O u "\$invcf")
              else
                  sel=(bcftools view -f PASS,. -O u "\$invcf")
              fi ;;
          *)       sel=(bcftools view -f PASS,. -O u "\$invcf") ;;
        esac
        "\${sel[@]}" | bcftools norm -m -any -f ${ref_fasta} -O u - \
          | bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t\${caller}\\n"
    }

    : > support.raw
    [ -s "${mutect2_vcf}" ]   && emit_support "${mutect2_vcf}"   mutect2  pass    >> support.raw
    [ -s "${strelka_snv}" ]   && emit_support "${strelka_snv}"   strelka2 pass    >> support.raw
    [ -s "${strelka_indel}" ] && emit_support "${strelka_indel}" strelka2 pass    >> support.raw
    [ -s "${muse_vcf}" ]      && emit_support "${muse_vcf}"      muse     muse    >> support.raw
    [ -s "${varscan_snv}" ]   && emit_support "${varscan_snv}"   varscan  varscan >> support.raw
    [ -s "${varscan_indel}" ] && emit_support "${varscan_indel}" varscan  varscan >> support.raw

    # aggregate - CHROM POS REF ALT CALLERS(csv) NCALLERS(distinct)
    sort -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 support.raw | awk -F'\\t' 'BEGIN{OFS="\\t"}
      { k=\$1 FS \$2 FS \$3 FS \$4
        if(k!=prev && prev!=""){ print pc,pp,pr,pa,cs,nc; cs="";nc=0 }
        if(\$5!="" && index(","cs",", ","\$5",")==0){ cs=(cs==""?\$5:cs","\$5); nc++ }
        prev=k;pc=\$1;pp=\$2;pr=\$3;pa=\$4 }
      END{ if(prev!="") print pc,pp,pr,pa,cs,nc }' > agg.tsv

    # Mutect2 MNPs (equal-length REF/ALT > 1 bp) are kept as ONE record (correct codon-level
    # annotation). Other callers report the constituent SNVs, so an MNP's vote = callers that
    # support EVERY constituent SNV (intersection) + mutect2. Deterministic caller order.
    awk -F'\\t' 'BEGIN { OFS = "\\t"; no = split("mutect2,strelka2,muse,varscan", ord, ",") }
      { k = \$1 FS \$2 FS \$3 FS \$4; cs[k] = \$5; row[NR] = \$0 }
      END {
        for (r = 1; r <= NR; r++) {
          split(row[r], f, "\\t")
          if (length(f[3]) > 1 && length(f[3]) == length(f[4]) && index("," f[5] ",", ",mutect2,") > 0) {
            for (o in inter) delete inter[o]; first = 1
            for (i = 1; i <= length(f[3]); i++) {
              rb = substr(f[3], i, 1); ab = substr(f[4], i, 1); if (rb == ab) continue
              for (o in cur) delete cur[o]
              m = split(cs[f[1] FS (f[2] + i - 1) FS rb FS ab], c, ",")
              for (j = 1; j <= m; j++) if (c[j] != "") cur[c[j]] = 1
              if (first) { for (o in cur) inter[o] = 1; first = 0 }
              else       { for (o in inter) if (!(o in cur)) delete inter[o] }
            }
            inter["mutect2"] = 1; out = ""; nc = 0
            for (x = 1; x <= no; x++) if (ord[x] in inter) { out = (out == "" ? ord[x] : out "," ord[x]); nc++ }
            print f[1], f[2], f[3], f[4], out, nc
          } else print row[r]
        }
      }' agg.tsv > agg.mnp.tsv && mv agg.mnp.tsv agg.tsv
    bgzip -f agg.tsv
    tabix -f -s1 -b2 -e2 agg.tsv.gz

    # distinct callers present (clamp thresholds so we never require more than exist)
    n_callers=\$(cut -f5 support.raw | sort -u | grep -c . || true)
    eff_snv=${min_snv};   [ "\$n_callers" -lt "\$eff_snv" ]   && eff_snv=\$n_callers
    eff_indel=${min_indel}; [ "\$n_callers" -lt "\$eff_indel" ] && eff_indel=\$n_callers
    [ "\$eff_snv"   -lt 1 ] && eff_snv=1
    [ "\$eff_indel" -lt 1 ] && eff_indel=1

    cat > cons_hdr.txt <<'HDR'
##INFO=<ID=CALLERS,Number=.,Type=String,Description="Somatic callers supporting this variant">
##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description="Number of somatic callers supporting this variant">
HDR

    # Backbone = Mutect2 (carries _T/_N samples + FORMAT). Normalize the SAME way,
    # keep PASS, annotate with the cross-caller vote, then keep records meeting the
    # per-type N-of-M threshold. Output stays a Mutect2-genotyped VCF.
    bcftools view -f PASS,. -O u ${mutect2_vcf} \
      | bcftools norm -m -any -f ${ref_fasta} -O z -o m2.norm.vcf.gz -
    bcftools index -f -t m2.norm.vcf.gz

    bcftools annotate -a agg.tsv.gz -h cons_hdr.txt \
        -c CHROM,POS,REF,ALT,INFO/CALLERS,INFO/NCALLERS \
        m2.norm.vcf.gz -O z -o m2.annot.vcf.gz
    bcftools index -f -t m2.annot.vcf.gz

    # split by type, apply per-type threshold, recombine (all Mutect2 - concat-safe)
    if [ "${filter_mode}" = "tag" ]; then
        bcftools view -v snps,mnps -O u m2.annot.vcf.gz | bcftools filter -m + -s LOWCONSENSUS -e "INFO/NCALLERS < \$eff_snv"   -O z -o part.snv.vcf.gz
        bcftools view -v indels -O u m2.annot.vcf.gz | bcftools filter -m + -s LOWCONSENSUS -e "INFO/NCALLERS < \$eff_indel" -O z -o part.indel.vcf.gz
    else
        bcftools view -v snps,mnps -i "INFO/NCALLERS >= \$eff_snv"   -O z -o part.snv.vcf.gz   m2.annot.vcf.gz
        bcftools view -v indels -i "INFO/NCALLERS >= \$eff_indel" -O z -o part.indel.vcf.gz m2.annot.vcf.gz
    fi
    bcftools index -f -t part.snv.vcf.gz
    bcftools index -f -t part.indel.vcf.gz

    bcftools concat -a part.snv.vcf.gz part.indel.vcf.gz -O u \
      | bcftools sort -O z -o ${meta.sample_id}_consensus.vcf.gz -
    bcftools index -f -t ${meta.sample_id}_consensus.vcf.gz

    {   echo "CONSENSUS SUMMARY: ${meta.sample_id}"
        echo "mode: N-of-M vote over all callers; Mutect2-backbone records; provenance in INFO/CALLERS,NCALLERS; filter=${filter_mode}"
        echo "callers_present=\$n_callers  eff_snv=\$eff_snv  eff_indel=\$eff_indel  muse_exclude=${muse_exclude}"
        echo "mutect2 PASS input:   \$(bcftools view -H -f PASS,. ${mutect2_vcf} | wc -l)"
        echo "consensus SNVs:       \$(bcftools view -H -v snps   ${meta.sample_id}_consensus.vcf.gz | wc -l)"
        echo "consensus MNPs:       \$(bcftools view -H -v mnps   ${meta.sample_id}_consensus.vcf.gz | wc -l)"
        echo "consensus indels:     \$(bcftools view -H -v indels ${meta.sample_id}_consensus.vcf.gz | wc -l)"
        echo "total consensus:      \$(bcftools view -H ${meta.sample_id}_consensus.vcf.gz | wc -l)"
    } > ${meta.sample_id}_consensus_stats.txt
    """
}


// ============================================================
//  STAGE 5P: SOMATIC variant-level filtering (VAF / depth / tumor-in-normal)
//
//  Soft-filters Mutect2/ensemble somatic SNV+indel calls on tumor VAF,
//  alt-read support, total depth, and presence in the matched normal
//  (tumor-in-normal). Tumor/normal sample columns are resolved BY NAME
//  (_T/_N suffix), never by assumed position, since Mutect2 may emit the
//  two samples in either column order. Thresholds are params (clinical
//  defaults). Mode 'tag' (default) keeps all variants and labels failures;
//  mode 'remove' drops them.
// ============================================================
process FILTER_SOMATIC_VCF {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-filtered) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/callers_consensus" }, mode: 'copy'


    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta),
          path("${meta.sample_id}_filtered.vcf.gz"),
          path("${meta.sample_id}_filtered.vcf.gz.tbi"), emit: vcf
    path  "${meta.sample_id}_filter_stats.txt",          emit: stats

    script:
    def min_t_vaf = params.filter_min_tumor_vaf    ?: 0.05
    def min_t_ad  = params.filter_min_tumor_alt_ad ?: 3
    def min_t_dp  = params.filter_min_tumor_dp     ?: 20
    def max_n_vaf = (params.filter_max_normal_vaf != null) ? params.filter_max_normal_vaf : 0.01   // 0 is a valid strict setting
    // Normal evidence needs >= this many alt reads as well as VAF: 1 read in a 30x normal is 3.3% VAF
    // (sequencing noise), while 1 read in an 800x normal passes the VAF test anyway.
    def min_n_alt = (params.filter_min_normal_alt_fail != null) ? params.filter_min_normal_alt_fail : 2
    def mode      = params.somatic_filter_mode     ?: 'tag'
    """
    set -euo pipefail

    # Resolve tumor/normal sample COLUMN INDEX by name suffix (_T/_N). Mutect2
    # may write the two samples in either order, so we key by NAME not position.
    samples="\$(bcftools query -l ${vcf})"
    tumor_idx=""; normal_idx=""; i=0
    while IFS= read -r s; do
        case "\$s" in
            *_T) tumor_idx=\$i ;;
            *_N) normal_idx=\$i ;;
        esac
        i=\$((i+1))
    done <<EOF
\$samples
EOF

    if [ -z "\$tumor_idx" ] || [ -z "\$normal_idx" ]; then
        echo "ERROR: could not resolve tumor(_T)/normal(_N) sample columns in ${vcf}" >&2
        echo "samples found: \$samples" >&2
        exit 1
    fi

    # Fail-condition: tumor VAF/altAD/DP below thresholds OR alt present in normal.
    # Scalar per-sample indexing (FMT/TAG[idx]) so '||' joins resolved values.
    # null-guard every term in the expression but be aware of the residual risk this leaves
    EXPR="(FMT/AF[\${tumor_idx}:0]!='.' && FMT/AF[\${tumor_idx}:0] < ${min_t_vaf}) || (FMT/AD[\${tumor_idx}:1]!='.' && FMT/AD[\${tumor_idx}:1] < ${min_t_ad}) || (FMT/DP[\${tumor_idx}]!='.' && FMT/DP[\${tumor_idx}] < ${min_t_dp}) || (FMT/AF[\${normal_idx}:0]!='.' && FMT/AF[\${normal_idx}:0] > ${max_n_vaf} && FMT/AD[\${normal_idx}:1]!='.' && FMT/AD[\${normal_idx}:1] >= ${min_n_alt})"
    if [ "${mode}" = "remove" ]; then
        bcftools view -e "\$EXPR" -O z -o ${meta.sample_id}_filtered.vcf.gz ${vcf}
    else
        bcftools filter -m + -s LOWSUPPORT -e "\$EXPR" -O z -o ${meta.sample_id}_filtered.vcf.gz ${vcf}
    fi
    bcftools index -t ${meta.sample_id}_filtered.vcf.gz

    {   echo "FILTER SUMMARY: ${meta.sample_id}"
        echo "mode: ${mode}"
        echo "thresholds: tumor_vaf>=${min_t_vaf} tumor_alt_ad>=${min_t_ad} tumor_dp>=${min_t_dp} normal_vaf<=${max_n_vaf} (fail only with normal_alt_ad>=${min_n_alt})"
        echo "tumor_col_idx: \$tumor_idx  normal_col_idx: \$normal_idx"
        echo "input variants:       \$(bcftools view -H ${vcf} | wc -l)"
        echo "output total records: \$(bcftools view -H ${meta.sample_id}_filtered.vcf.gz | wc -l)"
        echo "output PASS records:  \$(bcftools view -H -f PASS ${meta.sample_id}_filtered.vcf.gz | wc -l)"
    } > ${meta.sample_id}_filter_stats.txt
    """
}


// ============================================================
//  STAGE 5R-LOH: SEQUENZA - allele-specific CN + LOH + purity/ploidy
//  (Replaces FACETS) - PARALLEL version
//
//  --parallel runs chromosomes concurrently, writing ONE FILE PER CHROMOSOME
//  named <output-prefix><chromosome>, then merges them before binning.
//  VERIFY on first run that per-chrom files match the pattern in the merge
//  loop (f="<prefix><chrom>"); adjust if the actual suffix differs.
// ============================================================

process SEQUENZA {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    maxRetries    2

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/cnv/sequenza" }, mode: 'copy'

    input:
    tuple val(meta),
          path(tumour_bam), path(tumour_bai),
          path(normal_bam), path(normal_bai)
    path  ref_fasta
    path  ref_fai
    path  sequenza_gc_wiggle

    output:
    tuple val(meta),
          path("${meta.sample_id}_segments.txt"),
          path("${meta.sample_id}_confints_CP.txt"),
          path("${meta.sample_id}_alternative_solutions.txt"),
          path("${meta.sample_id}_genome_view.pdf"),
          path("${meta.sample_id}_CP_contours.pdf"),
          path("${meta.sample_id}_chromosome_view.pdf"),
          path("${meta.sample_id}_sequenza_cp_table.RData"),
          path("${meta.sample_id}_sequenza_extract.RData"),
          path("${meta.sample_id}_mutations.txt"),
          path("${meta.sample_id}_sequenza_log.txt"), emit: sequenza
    tuple val(meta), path("${meta.sample_id}_fit_qc.tsv"), emit: fit_qc

    script:
    def window   = params.sequenza_window   ?: 50
    def gamma    = params.sequenza_gamma    ?: 80
    def kmin     = params.sequenza_kmin     ?: 10
    def assembly = params.sequenza_assembly ?: 'hg38'
    // Map your 0 or 1 gender encoding to sequenzas female flag.
    // ADJUST THIS LINE to your convention (see note below).
    def female_flag = (meta.sex as String == '0') ? 'TRUE' : 'FALSE'
    def chroms   = params.sequenza_chroms   ?: 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'
    // Cap parallel bam2seqz workers (each holds deep-coverage pileup in memory).
    def slots    = Math.min((task.cpus ?: 4) as Integer, (params.sequenza_parallel ?: 4) as Integer)
    """
    set -euo pipefail

    trap '
        rm -f ${meta.sample_id}.seqz.gz*
        rm -f ${meta.sample_id}_*.seqz.gz
        rm -f ${meta.sample_id}.merged.seqz.gz
        rm -f ${meta.sample_id}.binned.seqz.gz
        rm -f run_sequenza.R
        rm -f ${meta.sample_id}_CN_bars.pdf
        rm -f ${meta.sample_id}_chromosome_depths.pdf
        rm -f ${meta.sample_id}_gc_plots.pdf
        rm -f ${meta.sample_id}_model_fit.pdf
        rm -f ${meta.sample_id}_alternative_fit.pdf
    ' EXIT

    # Step 1: BAMs + GC wiggle - raw seqz, parallelized per chromosome.
    # --parallel writes one file per chromosome named <prefix>_<chrom>.seqz.gz
    # -C restricts to chr1-22,X,Y (skips d1.vd1 decoy/viral/HLA contigs).
    sequenza-utils bam2seqz \\
        --normal   ${normal_bam} \\
        --tumor    ${tumour_bam} \\
        --fasta    ${ref_fasta} \\
        -gc        ${sequenza_gc_wiggle} \\
        -C         ${chroms} \\
        --parallel ${slots} \\
        --output   ${meta.sample_id}.seqz.gz \\
        -q 20 -N 30

    # Step 1b: merge per-chromosome seqz in chromosome order, one header.
    first=1
    for chr in ${chroms}; do
        f="${meta.sample_id}_\${chr}.seqz.gz"
        if [ ! -s "\$f" ]; then
            echo "WARN: no seqz for \${chr} (\$f) - skipping" >&2
            continue
        fi
        if [ "\$first" = "1" ]; then
            zcat "\$f"
            first=0
        else
            zcat "\$f" | awk 'NR > 1 || \$1 != "chromosome"'
        fi
    done | bgzip > ${meta.sample_id}.merged.seqz.gz

    # Fail loudly if the merge is empty (prevents ZeroDivisionError in binning)
    if [ \$(zcat ${meta.sample_id}.merged.seqz.gz 2>/dev/null | head -2 | wc -l) -lt 2 ]; then
        echo "ERROR: merged seqz empty - no per-chromosome files matched ${meta.sample_id}_<chr>.seqz.gz" >&2
        ls -la ${meta.sample_id}*.seqz.gz >&2 || true
        exit 1
    fi

    # Step 2: bin merged seqz to reduce memory for R analysis
    sequenza-utils seqz_binning \\
        --seqz   ${meta.sample_id}.merged.seqz.gz \\
        --window ${window} \\
        -o       ${meta.sample_id}.binned.seqz.gz


    # Step 3: Sequenza R analysis
    cat > run_sequenza.R <<'RSCRIPT'
    suppressMessages(library(sequenza))
    suppressMessages(library(iotools))

    # sequenza 3.0.0 calls chunk.apply(parallel=) in four internal functions;
    # iotools >=0.3 renamed that argument to CH.PARALLEL. Patch all four in
    # sequenza's own namespace. Verified: complete set, all patch cleanly.
    # No-op if iotools is old (still has 'parallel').
    if ("CH.PARALLEL" %in% names(formals(iotools::chunk.apply))) {
        .targets <- c("gc.sample.stats", "read.seqz", "read.seqz.chr", "sequenza.extract")
        for (.t in .targets) {
            .orig <- get(.t, envir = asNamespace("sequenza"))
            .src  <- deparse(.orig)
            .src  <- gsub("parallel = parallel", "CH.PARALLEL = parallel", .src, fixed = TRUE)
            .fn   <- eval(parse(text = .src))
            environment(.fn) <- asNamespace("sequenza")
            assignInNamespace(.t, .fn, ns = "sequenza")
        }
        cat("Patched", length(.targets), "sequenza functions to use CH.PARALLEL.\\n")
    }

    args      <- commandArgs(trailingOnly = TRUE)
    seqz_file <- args[1]
    prefix    <- args[2]
    assembly  <- args[3]
    gamma     <- as.integer(args[4])
    kmin      <- as.integer(args[5])

    set.seed(1234)

    seqz_data <- sequenza.extract(
        seqz_file,
        assembly = assembly,
        gamma    = gamma,
        kmin     = kmin,
        verbose  = TRUE
    )

    CP <- sequenza.fit(
        seqz_data,
        cellularity = seq(0.1, 1, 0.01),
        ploidy      = seq(1, 7, 0.1)
    )

    sequenza.results(
        sequenza.extract = seqz_data,
        cp.table         = CP,
        sample.id        = prefix,
        out.dir          = ".",
        female           = ${female_flag}
    )
    RSCRIPT

    Rscript --vanilla run_sequenza.R \\
        ${meta.sample_id}.binned.seqz.gz \\
        ${meta.sample_id} \\
        ${assembly} \\
        ${gamma} \\
        ${kmin}
    # Fit ambiguity QC: top-2 purity/ploidy solutions by SLPP. Reported, not judged -
    # review cases where the two solutions are close (common in hypodiploid/low-purity).
    awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[\$i] = i
                           if (!("cellularity" in h) || !("ploidy" in h) || !("SLPP" in h)) exit
                           next }
                 { print \$(h["cellularity"]) "\\t" \$(h["ploidy"]) "\\t" \$(h["SLPP"]) }' \\
        ${meta.sample_id}_alternative_solutions.txt \\
        | sort -t "\$(printf '\\t')" -k3,3gr | head -n 2 > fit_top2.tmp || true
    {   printf 'rank\\tcellularity\\tploidy\\tSLPP\\n'
        awk -F'\\t' '{ print NR "\\t" \$0 }' fit_top2.tmp
        awk -F'\\t' 'NR == 1 { c1 = \$1; p1 = \$2; l1 = \$3 } NR == 2 { c2 = \$1; p2 = \$2; l2 = \$3 }
                     END { if (NR >= 2) printf "delta_1_minus_2\\t%.3f\\t%.2f\\t%.3f\\n", c1 - c2, p1 - p2, l1 - l2
                           else print "delta_1_minus_2\\tNA\\tNA\\tNA" }' fit_top2.tmp
    } > ${meta.sample_id}_fit_qc.tsv
    rm -f fit_top2.tmp
    """
}

// ============================================================
//  STAGE 5R_A-LOH: SEQUENZA - gene level CN + LOH + purity/ploidy
// ============================================================

process SEQUENZA_GENE_LOH {
    tag { "${meta.sample_id}" }
    label 'process_low'

    publishDir { "${pubDir('variant_calling', true)}/cnv/sequenza" }, mode: 'copy'

    input:
    tuple val(meta), path(segments), path(confints)   // SEQUENZA.out.sequenza[1] segments, [2] confints_CP
    path  gene_bed
    path  arm_bed, stageAs: 'arm_bed/*'  // params.arm_bed (chr,start,end,arm e.g. 1p) or [] -> no arm table

    output:
    tuple val(meta), path("${meta.sample_id}_gene_loh.tsv"),   emit: gene_loh
    tuple val(meta), path("${meta.sample_id}_cn_summary.tsv"), emit: cn_summary
    tuple val(meta), path("${meta.sample_id}_arm_cn.tsv"),     optional: true, emit: arm_cn

    script:
    def arm_min = (params.arm_call_min_frac != null) ? params.arm_call_min_frac : 0.5
    """
    set -euo pipefail
    cat > gene_loh.R <<'GENELOH'
    args <- commandArgs(trailingOnly = TRUE)
    seg <- read.delim(args[1], sep="\\t", header=TRUE, stringsAsFactors=FALSE)
    req <- c("chromosome","start.pos","end.pos","CNt","A","B")
    if (!all(req %in% colnames(seg))) stop(sprintf("segments missing: %s", paste(setdiff(req,colnames(seg)),collapse=",")), call.=FALSE)
    bed <- read.delim(args[2], sep="\\t", header=FALSE, stringsAsFactors=FALSE)
    if (ncol(bed) < 5) stop("gene BED needs 5 cols: chr,start,end,symbol,ENSG", call.=FALSE)
    bed <- bed[,1:5]; colnames(bed) <- c("chr","start","end","gene","gene_id")
    nc <- function(x){ x <- sub("^chr","",as.character(x)); x[x %in% c("M","MT")] <- "MT"; x }
    seg\$chromosome <- nc(seg\$chromosome); bed\$chr <- nc(bed\$chr)
    bed\$start <- bed\$start + 1L
    for (cc in c("CNt","A","B","start.pos","end.pos")) seg[[cc]] <- suppressWarnings(as.numeric(seg[[cc]]))

    # Tumour ploidy = length-weighted mean autosomal CNt; expected CN baseline = round(ploidy).
    auto_re <- "^([1-9]|1[0-9]|2[0-2])\$"
    seglen  <- seg\$end.pos - seg\$start.pos + 1
    ok      <- grepl(auto_re, seg\$chromosome) & !is.na(seg\$CNt)
    ploidy  <- if (any(ok)) sum(seglen[ok] * seg\$CNt[ok]) / sum(seglen[ok]) else NA_real_
    base    <- if (is.na(ploidy)) 2 else floor(ploidy + 0.5)

    clf <- function(cn,b){ if(is.na(cn)&&is.na(b))return("NO_COVERAGE"); if(is.na(b))return("UNKNOWN"); if(!is.na(cn)&&cn==0)return("HOMDEL"); if(b==0)return("LOH"); "HET" }
    sbc <- split(seq_len(nrow(seg)), seg\$chromosome)
    rc <- rep(NA_real_, nrow(bed)); ra <- rc; rb <- rc
    for (i in seq_len(nrow(bed))) {
      idx <- sbc[[ bed\$chr[i] ]]; if (is.null(idx)) next
      ov <- pmin(bed\$end[i], seg\$end.pos[idx]) - pmax(bed\$start[i], seg\$start.pos[idx]) + 1
      hit <- ov > 0; if (!any(hit)) next
      j <- idx[hit][ which.max(ov[hit]) ]; rc[i] <- seg\$CNt[j]; ra[i] <- seg\$A[j]; rb[i] <- seg\$B[j]
    }
    out <- data.frame(Gene_ID=bed\$gene_id, Hugo_Symbol=bed\$gene, Chromosome=bed\$chr, Start=bed\$start, End=bed\$end,
                      Total_CN=rc, Major_Allele=ra, Minor_Allele=rb, stringsAsFactors=FALSE)
    out\$LOH_Status <- mapply(clf, out\$Total_CN, out\$Minor_Allele)
    # CN-neutral relative to tumour ploidy (not CN==2): a tetraploid A=4,B=0 is CN-neutral LOH.
    out\$LOH_Detail <- ifelse(out\$LOH_Status != "LOH", out\$LOH_Status,
                       ifelse(out\$Total_CN == base, "CN_neutral_LOH",
                       ifelse(out\$Total_CN <  base, "LOH_with_loss", "LOH_with_gain")))
    # Males: chrX/chrY are single-copy by nature - B==0 there is expected, not LOH.
    if (identical(args[10], "1")) {
      hemi <- out\$Chromosome %in% c("X", "Y") & out\$LOH_Status == "LOH"
      out\$LOH_Status[hemi] <- "HEMIZYGOUS"
      out\$LOH_Detail[hemi] <- "HEMIZYGOUS_EXPECTED"
    }
    write.table(out, args[3], sep="\\t", quote=FALSE, row.names=FALSE)

    # Genome-level summary: ploidy, fraction genome altered (CNt != baseline), LOH fraction (autosomes).
    L    <- sum(seglen[ok])
    fga  <- if (L > 0) sum(seglen[ok & seg\$CNt != base]) / L else NA_real_
    okb  <- ok & !is.na(seg\$B)
    floh <- if (sum(seglen[okb]) > 0) sum(seglen[okb & seg\$B == 0]) / sum(seglen[okb]) else NA_real_
    # Sequenza's own fit (confints_CP.txt, Sequenza 3.x: row 2 = point estimate).
    cp <- tryCatch(read.delim(args[9], sep="\\t", header=TRUE), error=function(e) NULL)
    s_cell <- if (!is.null(cp) && nrow(cp) >= 2 && "cellularity" %in% colnames(cp)) cp\$cellularity[2] else NA
    s_ploi <- if (!is.null(cp) && nrow(cp) >= 2 && "ploidy.estimate" %in% colnames(cp)) cp\$ploidy.estimate[2] else NA
    summ <- data.frame(metric = c("sample_id","tumour_ploidy_weighted","cn_baseline","autosomal_bp_assessed",
                                  "fraction_genome_altered","fraction_genome_loh",
                                  "sequenza_cellularity","sequenza_ploidy_estimate"),
                       value  = c(args[4], round(ploidy, 3), base, format(L, scientific=FALSE), round(fga, 4), round(floh, 4),
                                  s_cell, s_ploi))
    write.table(summ, args[5], sep="\\t", quote=FALSE, row.names=FALSE)

    # Arm-level calls: fraction of each autosomal arm lost / gained vs baseline, and LOH fraction.
    if (length(args) >= 6 && nzchar(args[6])) {
      minf <- as.numeric(args[7])
      arm  <- read.delim(args[6], sep="\\t", header=FALSE, stringsAsFactors=FALSE, comment.char="#")
      if (ncol(arm) < 4) stop("arm BED needs 4 cols: chr,start,end,arm", call.=FALSE)
      arm <- arm[,1:4]; colnames(arm) <- c("chr","start","end","arm")
      arm\$chr <- nc(arm\$chr); arm\$start <- arm\$start + 1L
      rows <- list()
      for (i in seq_len(nrow(arm))) {
        if (!grepl(auto_re, arm\$chr[i])) next
        idx <- sbc[[ arm\$chr[i] ]]; if (is.null(idx)) next
        ov <- pmin(arm\$end[i], seg\$end.pos[idx]) - pmax(arm\$start[i], seg\$start.pos[idx]) + 1
        h  <- ov > 0 & !is.na(seg\$CNt[idx]); if (!any(h)) next
        o  <- ov[h]; cn <- seg\$CNt[idx][h]; b <- seg\$B[idx][h]; cov <- sum(o)
        lf <- sum(o[cn < base]) / cov; gf <- sum(o[cn > base]) / cov; hf <- sum(o[!is.na(b) & b == 0]) / cov
        call <- if (lf >= minf) "LOSS" else if (gf >= minf) "GAIN" else "NEUTRAL"
        rows[[length(rows) + 1]] <- data.frame(arm=arm\$arm[i], chromosome=arm\$chr[i], covered_bp=format(cov, scientific=FALSE),
                                               loss_frac=round(lf,4), gain_frac=round(gf,4), loh_frac=round(hf,4), call=call)
      }
      tab <- if (length(rows)) do.call(rbind, rows) else
             data.frame(arm=character(), chromosome=character(), covered_bp=numeric(), loss_frac=numeric(),
                        gain_frac=numeric(), loh_frac=numeric(), call=character())
      write.table(tab, args[8], sep="\\t", quote=FALSE, row.names=FALSE)
    }
    GENELOH
    Rscript --vanilla gene_loh.R ${segments} ${gene_bed} ${meta.sample_id}_gene_loh.tsv \\
        ${meta.sample_id} ${meta.sample_id}_cn_summary.tsv \\
        "${arm_bed}" ${arm_min} ${meta.sample_id}_arm_cn.tsv ${confints} "${meta.sex ?: 'NA'}"
    rm -f gene_loh.R
    """
}

// ============================================================
//  STAGE 7: CONSOLIDATE_PATIENT_REPORT 
// ============================================================
process CONSOLIDATE_PATIENT_REPORT {
    tag { patient_id }
    label 'process_low'

    // -- 07 classification (analysis-separated) ------------------------
    publishDir { "${pubDir('classification', true)}/reports" }, mode: 'copy'

    input:
    tuple val(patient_id), val(sample_id), val(sex),
        path(classification, stageAs: 'class.tsv'),
        path(cnvkit_cns,     stageAs: 'cnv.cns'),
        path(sequenza_seg,   stageAs: 'seq_segments.txt'),
        path(sv_annot,       stageAs: 'sv_annot.tsv'),
        path(msi,            stageAs: 'msi_summary'),
        path(tmb,            stageAs: 'tmb.txt'),
        path(contam,         stageAs: 'contamination.table')

    output:
    tuple val(patient_id), path("${patient_id}_consolidated.tsv"), emit: report

    script:
    // Optional absolute-CN threshold for the AMP label; unset -> no AMP label (GAIN only).
    def amp_min = (params.cn_amp_min_cn != null) ? params.cn_amp_min_cn : ''
    def sex_lab = (sex?.toString() == '0') ? 'female' : (sex?.toString() == '1') ? 'male' : 'unknown'
    """
    set -uo pipefail
    REPORT=${patient_id}_consolidated.tsv

    {
      echo "## ============================================================"
      echo "## CONSOLIDATED VARIANT REPORT - PATIENT: ${patient_id}   (tumor sample: ${sample_id})"
      echo "## Generated: \$(date '+%Y-%m-%d %H:%M:%S')"
      echo "## NOTE: computational summary for review; not a final clinical tier."
      echo "## ============================================================"
    } > \$REPORT

    printf '\\n## [1] SOMATIC SNV/INDEL  (07_classification; PASS calls only)\\n' >> \$REPORT
    if [ -s class.tsv ]; then grep -v '^#' class.tsv >> \$REPORT; else echo "(no classification output)" >> \$REPORT; fi

    # CN calls relative to the EXPECTED copy number: autosomes = tumour ploidy baseline,
    # chrX/chrY sex-aware (male X/Y = baseline/2; female Y skipped; unknown sex -> X/Y skipped).
    # Source: Sequenza (purity/ploidy-fitted) when present, else CNVkit (NOT purity-adjusted).
    cat > cn_calls.awk <<'AWK'
BEGIN { FS = OFS = "\\t" }
function isauto(c) { return c ~ /^(chr)?([1-9]|1[0-9]|2[0-2])\$/ }
function isx(c)    { return c ~ /^(chr)?X\$/ }
function isy(c)    { return c ~ /^(chr)?Y\$/ }
# pass 1 (sequenza only): tumour ploidy = length-weighted mean autosomal CNt
NR == FNR && src == "sequenza" {
    if (FNR == 1) { for (i = 1; i <= NF; i++) h[\$i] = i; next }
    cn = \$(h["CNt"]); if (cn == "" || cn == "NA") next
    if (isauto(\$(h["chromosome"]))) { L = \$(h["end.pos"]) - \$(h["start.pos"]) + 1; sw += L * cn; sl += L }
    next
}
FNR == 1 {
    for (i = 1; i <= NF; i++) h[\$i] = i
    if (src == "sequenza") { ploidy = (sl > 0) ? sw / sl : 2; base = int(ploidy + 0.5)
                             cc = "chromosome"; sc = "start.pos"; ec = "end.pos"; nc = "CNt"
                             note = sprintf("tumour_ploidy=%.2f (length-weighted mean autosomal CNt)", ploidy) }
    else                   { base = 2; cc = "chromosome"; sc = "start"; ec = "end"; nc = "cn"
                             note = "NOT purity/ploidy-adjusted" }
    printf "## cn_source=%s  %s  autosomal_baseline=%d  sex=%s  amp_min_cn=%s\\n", src, note, base, sexlab, (ampmin == "" ? "unset" : ampmin)
    print "CHROM", "START", "END", "GENE", "LOG2", "COPY_NUMBER", "CALL"
    next
}
{
    c = \$(h[cc]); cn = \$(h[nc]); if (cn == "" || cn == "NA") next; cn += 0
    if (isauto(c))    ex = base
    else if (isx(c)) { if (sex == "1") ex = base / 2; else if (sex == "0") ex = base; else next }
    else if (isy(c)) { if (sex == "1") ex = base / 2; else next }
    else next
    if (cn == 0)                              call = "HOMDEL"
    else if (ampmin != "" && cn >= ampmin + 0) call = "AMP"
    else if (cn > ex)                         call = "GAIN"
    else if (cn < ex)                         call = "LOSS"
    else next
    if (src == "sequenza") { dr = \$(h["depth.ratio"]); lg = (dr + 0 > 0) ? sprintf("%.3f", log(dr) / log(2)) : "."; g = "." }
    else                   { lg = \$(h["log2"]); g = \$(h["gene"]) }
    print c, \$(h[sc]), \$(h[ec]), g, lg, cn, call
}
AWK
    if [ -s seq_segments.txt ]; then
      printf '\\n## [2] COPY NUMBER  (Sequenza segments; purity/ploidy-fitted; calls vs sex-aware expected CN)\\n' >> \$REPORT
      awk -v src=sequenza -v sex="${sex}" -v sexlab="${sex_lab}" -v ampmin="${amp_min}" -f cn_calls.awk seq_segments.txt seq_segments.txt >> \$REPORT
    elif [ -s cnv.cns ]; then
      printf '\\n## [2] COPY NUMBER  (CNVkit call.cns; NOT purity-adjusted; calls vs sex-aware expected CN)\\n' >> \$REPORT
      awk -v src=cnvkit -v sex="${sex}" -v sexlab="${sex_lab}" -v ampmin="${amp_min}" -f cn_calls.awk cnv.cns >> \$REPORT
    else
      printf '\\n## [2] COPY NUMBER\\n(no Sequenza or CNVkit output)\\n' >> \$REPORT
    fi

    printf '\\n## [3] LOH  (Sequenza; minor allele B==0)\\n' >> \$REPORT
    if [ -s seq_segments.txt ]; then
      awk -F'\\t' 'BEGIN{OFS="\\t"; print "CHROM","START","END","CNt","A","B","STATUS"}
        NR==1{ for(i=1;i<=NF;i++) h[\$i]=i; next }
        { b=\$(h["B"]); if(b=="")next;
          if(b+0==0 && \$(h["CNt"])!="") print \$(h["chromosome"]),\$(h["start.pos"]),\$(h["end.pos"]),\$(h["CNt"]),\$(h["A"]),b,"LOH" }' seq_segments.txt >> \$REPORT
    else echo "(no Sequenza output)" >> \$REPORT; fi

    printf '\\n## [4] STRUCTURAL VARIANTS  (Manta + AnnotSV; one row/SV, PASS. ACMG = germline framework, advisory)\\n' >> \$REPORT
    if [ -s sv_annot.tsv ] && [ \$(grep -vc '^AnnotSV_ID' sv_annot.tsv) -gt 0 ]; then
      awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++)h[\$i]=i;print"CHROM","START","END","SVTYPE","GENE(S)","ACMG";next}
        (("Annotation_mode" in h)?\$(h["Annotation_mode"]):"full")!="full"{next}
        (("FILTER" in h)&&\$(h["FILTER"])!="PASS"&&\$(h["FILTER"])!="."){next}
        {g=(("Gene_name" in h)?\$(h["Gene_name"]):".");g=(g==""?".":g)
         a=(("ACMG_class" in h)?\$(h["ACMG_class"]):".");a=(a==""||a=="NA"?".":a)
         print \$(h["SV_chrom"]),\$(h["SV_start"]),\$(h["SV_end"]),\$(h["SV_type"]),g,a}' OFS='\\t' sv_annot.tsv >> \$REPORT
    else echo "(no annotated SV output)" >> \$REPORT; fi

    printf '\\n## [5] MSI  (MSIsensor-pro)\\n' >> \$REPORT
    if [ -s msi_summary ]; then cat msi_summary >> \$REPORT; else echo "(no MSI output)" >> \$REPORT; fi

    printf '\\n## [6] TMB\\n' >> \$REPORT
    if [ -s tmb.txt ]; then cat tmb.txt >> \$REPORT; else echo "(no TMB output)" >> \$REPORT; fi

    printf '\\n## [7] TUMOUR CROSS-SAMPLE CONTAMINATION  (GATK CalculateContamination; threshold %s)\\n' "${params.qc_max_contamination ?: 0.03}" >> \$REPORT
    if [ -s contamination.table ]; then
        awk -F'\\t' -v max="${params.qc_max_contamination ?: 0.03}" 'NR == 1 { next }
            { printf "sample\\t%s\\ncontamination\\t%s\\nerror\\t%s\\nstatus\\t%s\\n", \$1, \$2, \$3, (\$2 + 0 > max + 0 ? "WARN" : "PASS") }' \\
            contamination.table >> \$REPORT
    else echo "(no contamination estimate - pileup_variants not set)" >> \$REPORT; fi
    """
}