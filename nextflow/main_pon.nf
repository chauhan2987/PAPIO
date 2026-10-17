#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
//  main_pon.nf  -  Standalone Mutect2 (somatic) Panel-of-Normals (PON) builder
//
//  Deliberately NOT wired into main.nf: building a PON is a one-time (or
//  occasional) reference-curation task, not a per-patient clinical run - it
//  has its own samplesheet shape (normal BAMs only, no patient/tumor pairing)
//  and shouldn't share main's assay/stage-gating machinery. Run it once,
//  publish inhouse_pon.vcf.gz, then point params.panel_of_normals at it
//  for every future somatic MUTECT2 run.
//
//  Run:
//    nextflow run main_pon.nf -c pon.config -params-file params_pon.yaml\
//        -profile production
//
//  Input samplesheet (simple, 3 columns):
//    sample_id,bam,bai
//    S001_N,/path/S001_N.bam,/path/S001_N.bam.bai
//    S002_N,/path/S002_N.bam,/path/S002_N.bam.bai
// ============================================================

include { gatkJavaOpts } from './utils'
include {
    PON_MUTECT2; PON_MERGE_NORMAL;
    PON_GENOMICSDB_IMPORT; PON_CREATE_SHARD; PON_MERGE_FINAL
} from './components_pon'

include { toFileChannel; toOptionalFileChannel; stripChunkMeta } from './workflows_dna'


// ============================================================
//  GENERATE_PON_META_INFO - minimal meta factory for normal-only input
// ============================================================
workflow GENERATE_PON_META_INFO {

    take:
    samplesheet

    main:
    def REQUIRED_COLUMNS = ['sample_id', 'bam', 'bai']

    ch_normals = channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            REQUIRED_COLUMNS.each { col ->
                if (!row.containsKey(col) || row[col] == null || row[col] == '') {
                    error "PON samplesheet row missing required column '${col}': ${row}"
                }
            }
            def meta = [sample_id: row.sample_id]
            def bam  = file(row.bam, checkIfExists: true)
            def bai  = file(row.bai, checkIfExists: true)
            tuple(meta, bam, bai)
        }

    emit:
    normals = ch_normals
}


// ============================================================
//  PON_WORKFLOW - entry sub-workflow
// ============================================================
workflow PON_WORKFLOW {

    main:
    if (!params.pon_samplesheet)
        error "ERROR: params.pon_samplesheet is required (CSV: sample_id,bam,bai - normal samples only)."
    if (!params.scatter_intervals)
        error "ERROR: params.scatter_intervals is required (dir of interval-list shards, same as the somatic Mutect2 arm)."

    // -- Reference channels (value - broadcast to all tasks) ---
    ch_ref_fasta = toFileChannel(params.ref_fasta)
    ch_ref_fai   = toFileChannel(params.ref_fai)
    ch_ref_dict  = toFileChannel(params.ref_dict)

    ch_germline_resource     = toOptionalFileChannel(params.germline_resource)
    ch_germline_resource_tbi = params.germline_resource ? toFileChannel("${params.germline_resource}.tbi") : channel.value([])

    // interval shards - SAME directory the somatic MUTECT2 scatter uses
    ch_intervals = channel.fromPath("${params.scatter_intervals}/*", checkIfExists: true)
        .map { f -> tuple(f.baseName, f) }

    // -- normals - intervals cross product (scatter) ----------------------
    ch_normals = GENERATE_PON_META_INFO(params.pon_samplesheet).normals

    ch_scatter_in = ch_normals.combine(ch_intervals)
        .map { meta, bam, bai, iname, ifile -> tuple(meta, bam, bai, iname, ifile) }

    // -- STAGE P1: Mutect2 tumor-only per normal per interval -------------
    PON_MUTECT2(ch_scatter_in, ch_ref_fasta, ch_ref_fai, ch_ref_dict)

    // -- STAGE P2: gather chunks -> one VCF per normal ---------------------
    ch_per_normal_chunks = PON_MUTECT2.out.vcf_chunk
        .map { meta, iname, vcf, tbi -> tuple(meta.sample_id, meta, vcf, tbi) }
        .groupTuple()
        .map { sid, metas, vcfs, tbis -> tuple(metas[0], vcfs, tbis) }

    PON_MERGE_NORMAL(ch_per_normal_chunks)

    // -- STAGE P3: GenomicsDBImport (all normals, scattered by interval) --
    ch_all_normal_vcfs = PON_MERGE_NORMAL.out.vcf.map { meta, vcf, tbi -> vcf }.collect()
    ch_all_normal_tbis = PON_MERGE_NORMAL.out.vcf.map { meta, vcf, tbi -> tbi }.collect()

    PON_GENOMICSDB_IMPORT(ch_intervals, ch_all_normal_vcfs, ch_all_normal_tbis,
                           ch_ref_fasta, ch_ref_fai, ch_ref_dict)

    // -- STAGE P4: CreateSomaticPanelOfNormals per shard -------------------
    PON_CREATE_SHARD(PON_GENOMICSDB_IMPORT.out.db, ch_ref_fasta, ch_ref_fai, ch_ref_dict, ch_germline_resource, ch_germline_resource_tbi)

    // -- STAGE P5: merge shard PONs -> final panel_of_normals.vcf.gz ------
    ch_pon_vcfs = PON_CREATE_SHARD.out.vcf_chunk.map { iname, vcf, tbi -> vcf }.collect()
    ch_pon_tbis = PON_CREATE_SHARD.out.vcf_chunk.map { iname, vcf, tbi -> tbi }.collect()

    PON_MERGE_FINAL(ch_pon_vcfs, ch_pon_tbis, ch_ref_dict)
}


workflow {
    PON_WORKFLOW()
}
