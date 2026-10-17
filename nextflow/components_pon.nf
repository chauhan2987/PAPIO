#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
//  components_pon.nf  -   Standalone Mutect2 (somatic) Panel-of-Normals (PON) builder
// ============================================================

include { gatkJavaOpts; pubDir } from './utils'


// ============================================================
//  STAGE P1 (scatter): Mutect2 tumor-only mode per normal, per interval
//  NOTE: --max-mnp-distance 0 is NOT optional - GATK's own docs flag this
//  as a required workaround for a GenomicsDBImport bug on MNP records.
// ============================================================
process PON_MUTECT2 {
    tag { "${meta.sample_id}:${interval_name}" }
    label 'process_variant'

    maxRetries 2

    input:
    tuple val(meta), path(normal_bam), path(normal_bai),
          val(interval_name), path(interval_file)
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(meta), val(interval_name),
          path("${meta.sample_id}.${interval_name}.pon.vcf.gz"),
          path("${meta.sample_id}.${interval_name}.pon.vcf.gz.tbi"), emit: vcf_chunk

    script:
    def java_opts = gatkJavaOpts(task)
    """
    set -euo pipefail

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_${interval_name}_ponm2_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" Mutect2 \\
        -R ${ref_fasta} \\
        -I ${normal_bam} \\
        -L ${interval_file} \\
        --max-mnp-distance 0 \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}.${interval_name}.pon.vcf.gz
    """
}


// ============================================================
//  STAGE P2 (gather): merge one normal's scattered chunks into one VCF
//  No FilterMutectCalls here - GATK's PON recipe feeds RAW tumor-only calls
//  straight into GenomicsDBImport; filtering is not part of this workflow.
// ============================================================
process PON_MERGE_NORMAL {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // -- 00 reference / pon build artifacts -----------------------------
    publishDir { "${params.outdir}/00_reference/pon/per_normal" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf_chunks), path(tbi_chunks)

    output:
    tuple val(meta),
          path("${meta.sample_id}_for_pon.vcf.gz"),
          path("${meta.sample_id}_for_pon.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts   = gatkJavaOpts(task)
    def vcf_inputs  = vcf_chunks.collect { v -> "-I ${v}" }.join(' \\\n        ')
    """
    set -euo pipefail

    gatk --java-options "${java_opts}" MergeVcfs \\
        ${vcf_inputs} -O ${meta.sample_id}_for_pon.vcf.gz
    """
}


// ============================================================
//  STAGE P3: GenomicsDBImport - combine ALL normals' VCFs into one DB
//  Scattered by the SAME interval shards as PON_MUTECT2 (keeps DB build
//  parallel and memory-bounded; --batch-size caps concurrent sample readers,
//  which is what actually controls GenomicsDBImport's memory footprint on
//  large cohorts - NOT task memory alone).
// ============================================================
process PON_GENOMICSDB_IMPORT {
    tag { "${interval_name}" }
    label 'process_genomicsdb'

    maxRetries 2

    input:
    tuple val(interval_name), path(interval_file)
    path  normal_vcfs   // all per-normal *_for_pon.vcf.gz, collected
    path  normal_tbis   // all per-normal *_for_pon.vcf.gz.tbi, collected
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(interval_name), path("pon_db_${interval_name}"), emit: db

    script:
    def java_opts  = gatkJavaOpts(task)
    def vcf_inputs = normal_vcfs.collect { v -> "-V ${v}" }.join(' \\\n        ')
    def batch_size = params.pon_genomicsdb_batch_size ?: 50
    """
    set -euo pipefail

    tmp_dir="\$PWD/pon_gdb_${interval_name}_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" GenomicsDBImport \\
        -R ${ref_fasta} \\
        -L ${interval_file} \\
        --merge-input-intervals \\
        --genomicsdb-workspace-path pon_db_${interval_name} \\
        --batch-size ${batch_size} \\
        --tmp-dir "\$tmp_dir" \\
        ${vcf_inputs}
    """
}


// ============================================================
//  STAGE P4: CreateSomaticPanelOfNormals per interval shard
//  --germline-resource (gnomAD) lets the tool exclude likely-germline sites
//  (--max-germline-probability, default 0.5) - real germline filtering still
//  belongs to Mutect2's own --germline-resource model at call time; this
//  just keeps the PON itself focused on technical artifacts, per GATK docs.
// ============================================================
process PON_CREATE_SHARD {
    tag { "${interval_name}" }
    label 'process_medium'

    maxRetries 2

    input:
    tuple val(interval_name), path(pon_db)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  germline_resource
    path  germline_resource_tbi

    output:
    tuple val(interval_name),
          path("pon_${interval_name}.vcf.gz"),
          path("pon_${interval_name}.vcf.gz.tbi"), emit: vcf_chunk

    script:
    def java_opts  = gatkJavaOpts(task)
    def gr_flag    = germline_resource ? "--germline-resource ${germline_resource}" : ''
    def min_sample = params.pon_min_sample_count ?: 2   // GATK default: site must recur in >=2 normals
    """
    set -euo pipefail

    tmp_dir="\$PWD/pon_create_${interval_name}_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" CreateSomaticPanelOfNormals \\
        -R ${ref_fasta} \\
        -V gendb://${pon_db} \\
        --min-sample-count ${min_sample} \\
        ${gr_flag} \\
        --tmp-dir "\$tmp_dir" \\
        -O pon_${interval_name}.vcf.gz
    """
}


// ============================================================
//  STAGE P5 (gather): merge shard PONs into the final panel_of_normals.vcf.gz
// ============================================================

process PON_MERGE_FINAL {
    tag { "final_pon" }
    label 'process_medium'

    scratch false

    publishDir { "${params.outdir}/00_reference/final_pon" }, mode: 'copy'

    input:
    path pon_vcf_chunks
    path pon_tbi_chunks
    path ref_dict

    output:
    path "panel_of_normals.vcf.gz",     emit: vcf
    path "panel_of_normals.vcf.gz.tbi", emit: tbi

    script:
    def java_opts  = gatkJavaOpts(task)
    def vcf_inputs = pon_vcf_chunks.collect { v -> "-I ${v}" }.join(' \\\n        ')
    """
    set -euo pipefail

    gatk --java-options "${java_opts}" SortVcf \\
        ${vcf_inputs} \\
        --SEQUENCE_DICTIONARY ${ref_dict} \\
        -O panel_of_normals.vcf.gz
    """
}