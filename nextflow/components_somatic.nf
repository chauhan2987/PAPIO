#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Shared JVM-options helper (single source of truth in utils).
include { gatkJavaOpts; pubDir } from './utils'

// ============================================================
//  components_somatic.nf  -  Somatic (tumor-normal) process definitions
// ============================================================


// ============================================================
//  STAGE 5 variant_calling: MANTA - somatic SVs and candidate indels - Manta
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
//  STAGE 5 variant_calling: STRELKA2 - somatic SNVs/indels - Strelka2
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
//  STAGE 5 variant_calling: TAG_STRELKA_AF - Strelka2 AD/AF fields
//
//  Inject standard AF/AD FORMAT fields into Strelka2 somatic VCFs.
//  Strelka does not emit AF/AD natively. Per Strelka's user guide:
//    SNVs:   tier1 ref/alt counts come from FORMAT/{A,C,G,T}U (field matching
//            REF/ALT base), tier1 = first comma-delimited value.
//    Indels: tier1 ref/alt counts come from FORMAT/TAR (ref) and FORMAT/TIR (alt).
//    AF = alt_tier1 / (alt_tier1 + ref_tier1)
// ============================================================

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


    strelka_add_af.py ${strelka_snv}   snvs.af.vcf   snv
    strelka_add_af.py ${strelka_indel} indels.af.vcf indel

    bgzip -c snvs.af.vcf   > ${meta.sample_id}_strelka2_snvs.af.vcf.gz
    bgzip -c indels.af.vcf > ${meta.sample_id}_strelka2_indels.af.vcf.gz
    bcftools index -t ${meta.sample_id}_strelka2_snvs.af.vcf.gz
    bcftools index -t ${meta.sample_id}_strelka2_indels.af.vcf.gz
    """
}


// ============================================================
//  STAGE 5 variant_calling: MUSE_CALL - somatic SNVs - MuSE call
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
    """
        set -euo pipefail
        tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_muse_tmp.XXXXXX")
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
//  STAGE 5 variant_calling: MUSE_SUMP - somatic SNVs - MuSE sump
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
    def mode_flag = (meta.assay == 'wes') ? '-E' : '-G'

    """
    set -euo pipefail
    tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_muse_sump_tmp.XXXXXX")
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
//  STAGE 5 variant_calling: TAG_MUSE_AF - MuSE AF field
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


    bcftools norm -m -any -f ${ref_fasta} -O v ${muse_vcf} > split.vcf
    muse_add_af.py split.vcf muse.af.vcf
    bgzip -c muse.af.vcf > ${meta.sample_id}_muse.af.vcf.gz
    bcftools index -t ${meta.sample_id}_muse.af.vcf.gz
    """
}


// ============================================================
//  STAGE 5 variant_calling: VARSCAN - somatic SNVs/indels - VarScan2 (per shard on WGS)
//
//  WES: one shard ('all') restricted to the capture BED - identical to the
//  former single-task run. WGS: one task per scatter shard, because a single
//  samtools mpileup | varscan pipe is single-threaded and takes days on WGS.
// ============================================================

process VARSCAN {
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
    path  interval_bed

    output:
    tuple val(meta),
          path("${meta.sample_id}.${interval_name}.varscan.snp.vcf"),
          path("${meta.sample_id}.${interval_name}.varscan.indel.vcf"), emit: vcf_chunk

    script:
    def min_cov     = params.varscan_min_coverage    ?: 8
    def min_var_frq = params.varscan_min_var_freq    ?: 0.05
    def somatic_p   = params.varscan_somatic_pvalue  ?: 0.05
    def prefix      = "${meta.sample_id}.${interval_name}.varscan"
    // Region: WGS shard (interval_list or BED) | WES capture BED | none.
    def region_prep = ''
    def region_flag = ''
    if (interval_file) {
        // a BED shard is used as is (copying it to region.bed fails when it is already named so)
        region_prep = interval_file.name.endsWith('.bed')
            ? ''
            : "grep -v '^@' ${interval_file} | awk 'BEGIN{OFS=\"\\t\"} {print \$1, \$2 - 1, \$3}' > region.bed"
        region_flag = interval_file.name.endsWith('.bed') ? "-l ${interval_file}" : "-l region.bed"
    } else if (meta.assay == 'wes' && interval_bed) {
        region_flag = "-l ${interval_bed}"
    }
    """
    set -euo pipefail

    # samtools mpileup expects <bam>.bai; our index is <sample>_markdup.bai.
    if [ ! -f "${tumor_bam}.bai" ];  then ln -s ${tumor_bai}  "${tumor_bam}.bai";  fi
    if [ ! -f "${normal_bam}.bai" ]; then ln -s ${normal_bai} "${normal_bam}.bai"; fi
    ${region_prep}

    # VarScan somatic reads a COMBINED normal+tumor mpileup on stdin (--mpileup 1).
    # -B (no BAQ), as the VarScan documentation prescribes: BAQ lowers base qualities next
    # to indels and would hide real somatic SNVs and indel-adjacent alleles from VarScan.
    samtools mpileup \\
        -f ${ref_fasta} \\
        ${region_flag} \\
        -B -q 20 -Q 20 \\
        ${normal_bam} ${tumor_bam} \\
    | varscan somatic \\
        /dev/stdin \\
        ${prefix} \\
        --mpileup 1 \\
        --min-coverage ${min_cov} \\
        --min-var-freq ${min_var_frq} \\
        --somatic-p-value ${somatic_p} \\
        --output-vcf 1

    # VarScan writes <prefix>.snp.vcf / <prefix>.indel.vcf; guarantee both exist.
    for kind in snp indel; do
        [ -e "${prefix}.\${kind}.vcf" ] || : > "${prefix}.\${kind}.vcf"
    done
    """
}


// ============================================================
//  STAGE 5 variant_calling: VARSCAN_MERGE - VarScan2 shards merged in reference order
// ============================================================

process VARSCAN_MERGE {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/varscan" }, mode: 'copy'

    input:
    tuple val(meta), path(snv_chunks), path(indel_chunks)
    path  ref_fai

    output:
    tuple val(meta),
          path("${meta.sample_id}_varscan_snv.vcf.gz"),
          path("${meta.sample_id}_varscan_snv.vcf.gz.tbi"),
          path("${meta.sample_id}_varscan_indel.vcf.gz"),
          path("${meta.sample_id}_varscan_indel.vcf.gz.tbi"), emit: vcf

    script:
    def snv_list   = (snv_chunks   instanceof List ? snv_chunks   : [snv_chunks]).join(' ')
    def indel_list = (indel_chunks instanceof List ? indel_chunks : [indel_chunks]).join(' ')
    """
    gather_vcf_shards.sh ${meta.sample_id}_varscan_snv.vcf.gz   ${ref_fai} ${snv_list}
    gather_vcf_shards.sh ${meta.sample_id}_varscan_indel.vcf.gz ${ref_fai} ${indel_list}
    """
}


// ============================================================
//  STAGE 5 variant_calling: TAG_VARSCAN_AF - VarScan2 AD/AF fields
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
    varscan_tag_af.sh ${meta.sample_id} ${ref_fasta} ${varscan_snv} ${varscan_indel}
    """
}


// ============================================================
//  STAGE 5 variant_calling: MUTECT2 - somatic SNVs/indels - Mutect2 (per shard)
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
    def normal_name = meta.normal_id
    def gr_flag = germline_resource ? "--germline-resource ${germline_resource}" : ''

    """
    set -euo pipefail

    tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_mutect2_tmp.XXXXXX")
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


// ============================================================
//  STAGE 5 variant_calling: MERGE_MUTECT2 - Mutect2 shards merged and FilterMutectCalls
// ============================================================

process MERGE_MUTECT2 {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/snv/mutect2" }, mode: 'copy'

    input:
    // arity '1..*': always a list, so the -I / -stats builders below also work for a
    // single-shard scatter (a lone Path would be iterated over its path components).
    tuple val(meta), path(vcf_chunks, arity: '1..*'), path(tbi_chunks, arity: '1..*'), path(stats_chunks, arity: '1..*')
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
//  STAGE 5 variant_calling: ENSEMBLE_CONSENSUS - caller consensus (Mutect2-anchored or N-of-M)
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
    def filter_mode  = params.consensus_filter_mode       ?: 'tag'      // 'tag' | 'remove'
    def backbone     = params.consensus_backbone          ?: 'mutect2'  // 'mutect2' | 'union'
    """
    ensemble_consensus.sh --sample ${meta.sample_id} --normal ${meta.normal_id} --fasta ${ref_fasta} \\
        --mutect2 ${mutect2_vcf} \\
        ${strelka_snv   ? "--strelka-snv ${strelka_snv}"     : ''} ${strelka_indel ? "--strelka-indel ${strelka_indel}" : ''} \\
        ${muse_vcf      ? "--muse ${muse_vcf}"               : ''} \\
        ${varscan_snv   ? "--varscan-snv ${varscan_snv}"     : ''} ${varscan_indel ? "--varscan-indel ${varscan_indel}" : ''} \\
        --min-snv ${min_snv} --min-indel ${min_indel} --muse-exclude ${muse_exclude} \\
        --filter-mode ${filter_mode} --backbone ${backbone}
    """
}


// ============================================================
//  STAGE 5 variant_calling: FILTER_SOMATIC_VCF - somatic filters (VAF, depth, tumor-in-normal)
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
    def min_t_vaf = (params.filter_min_tumor_vaf    != null) ? params.filter_min_tumor_vaf    : 0.05   // 0 disables the test
    def min_t_ad  = (params.filter_min_tumor_alt_ad != null) ? params.filter_min_tumor_alt_ad : 3
    def min_t_dp  = (params.filter_min_tumor_dp     != null) ? params.filter_min_tumor_dp     : 20
    def max_n_vaf = (params.filter_max_normal_vaf != null) ? params.filter_max_normal_vaf : 0.01   // 0 is a valid strict setting
    // Normal evidence needs >= this many alt reads as well as VAF: 1 read in a 30x normal is 3.3% VAF
    // (sequencing noise), while 1 read in an 800x normal passes the VAF test anyway.
    def min_n_alt = (params.filter_min_normal_alt_fail != null) ? params.filter_min_normal_alt_fail : 2
    def mode      = params.somatic_filter_mode     ?: 'tag'
    """
    filter_somatic_vcf.sh --vcf ${vcf} --tumor ${meta.sample_id} --normal ${meta.normal_id} --mode ${mode} \\
        --min-tumor-vaf ${min_t_vaf} --min-tumor-ad ${min_t_ad} --min-tumor-dp ${min_t_dp} \\
        --max-normal-vaf ${max_n_vaf} --min-normal-alt ${min_n_alt}
    """
}


// ============================================================
//  STAGE 5 variant_calling: SEQUENZA - allele-specific copy number, purity and ploidy - Sequenza
//
//  sequenza-utils bam2seqz --parallel writes one seqz file per chromosome
//  (<prefix>_<chrom>.seqz.gz); they are merged in chromosome order before
//  binning and fitting.
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
    def seg_filt = (params.sequenza_segment_filter != null) ? params.sequenza_segment_filter : 3000000
    def assembly = params.sequenza_assembly ?: 'hg38'
    // Samplesheet sex (0 = female, 1 = male); NA is inferred from the normal's chrX depth in
    // bin/sequenza.sh. Sequenza leaves chrX/chrY out of the purity/ploidy fit for males.
    def female_flag = (meta.sex?.toString() == '0') ? 'TRUE' : (meta.sex?.toString() == '1') ? 'FALSE' : 'NA'
    def sex_fmin    = params.sex_female_min_xratio ?: 0.80
    def sex_mmax    = params.sex_male_max_xratio   ?: 0.65
    // Primary chromosomes as named in the reference (chr1 or 1 style); a fixed chr-prefixed
    // list matched nothing on references named 1..22,X,Y.
    def chroms_param = params.sequenza_chroms ?: ''
    // Cap parallel bam2seqz workers (each holds deep-coverage pileup in memory).
    def slots    = Math.min((task.cpus ?: 4) as Integer, (params.sequenza_parallel ?: 4) as Integer)
    """
    sequenza.sh --sample ${meta.sample_id} --tumor ${tumour_bam} --normal ${normal_bam} \\
        --fasta ${ref_fasta} --fai ${ref_fai} --gc-wiggle ${sequenza_gc_wiggle} --parallel ${slots} \\
        --window ${window} --gamma ${gamma} --kmin ${kmin} --female ${female_flag} \\
        --female-min-xratio ${sex_fmin} --male-max-xratio ${sex_mmax} \\
        --segment-filter ${seg_filt} --assembly ${assembly} ${chroms_param ? "--chroms \"${chroms_param}\"" : ''}
    """
}


// ============================================================
//  STAGE 5 variant_calling: SEQUENZA_GENE_LOH - gene- and arm-level copy number and LOH
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
    Rscript --vanilla "\$(command -v gene_loh.R)" ${segments} ${gene_bed} ${meta.sample_id}_gene_loh.tsv \\
        ${meta.sample_id} ${meta.sample_id}_cn_summary.tsv \\
        "${arm_bed}" ${arm_min} ${meta.sample_id}_arm_cn.tsv ${confints} "${meta.sex ?: 'NA'}"
    """
}


// ============================================================
//  STAGE 7 classification: CONSOLIDATE_PATIENT_REPORT - consolidated somatic report per tumor
// ============================================================

process CONSOLIDATE_PATIENT_REPORT {
    tag { patient_id }
    label 'process_low'

    // -- 07 classification (analysis-separated) ------------------------
    publishDir { "${pubDir('classification', true)}/reports" }, mode: 'copy'

    input:
    // patient_id = report name: the patient for a single tumour, the tumour sample_id when a
    // patient has several (samplesheet 'sample' column) so reports never overwrite each other.
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
    patient_report.sh --patient ${patient_id} --sample ${sample_id} --sex "${sex}" --sex-label ${sex_lab} \\
        --max-contam ${params.qc_max_contamination ?: 0.03} ${amp_min != '' ? "--amp-min ${amp_min}" : ''}
    """
}