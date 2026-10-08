#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
//  utils.nf  -  Channel Factory + Helper Sub-Workflows
// ============================================================



// ============================================================
//  SHARED CONSTANTS  -  THE single source of truth
//
//  Previously stageRankMap()'s contents were declared TWICE:
//  stageRank() in main_v3.nf AND a local STAGE_RANK map inside
//  DNA_WORKFLOW. Declared once here; both import it. Insert a stage in
//  one place and validation + execution stay in lockstep.
// ============================================================
def stageRankMap() {
    [ 'mapping': 1, 'markdup': 2, 'metrics': 3,
      'bqsr': 4, 'variant_calling': 5, 'annotation': 6,
      'classification': 7 ]
}

def validStages()     { ['mapping', 'markdup', 'metrics', 'bqsr', 'variant_calling', 'annotation', 'classification'] }
def validAssays()     { ['wes', 'wgs', 'rna'] }
def validAnalyses()   { ['germline', 'somatic', 'tumor_only', 'rna_fusion'] }
def validCallers()    { ['haplotypecaller', 'deepvariant'] }
def validDedupModes() { ['umi_consensus', 'umi_markdup', 'markdup', 'external'] }

// Single stage-window predicate (ranks from the shared stageRankMap()), so
// no gate hard-codes a rank literal. Identical mechanism to the DNA arm.
def stageRuns(String s) {
    def r = stageRankMap()[s]
    stageRankMap()[params.start_at] <= r && r <= stageRankMap()[params.stop_at]
}

// ============================================================
//  JVM / SCRATCH HELPERS  -  per-process boilerplate, declared once
//
//  Identical content to the formulas/blocks previously inlined in every
//  GATK/Picard/fgbio process body. A JVM-flag or scratch-policy change
//  is now a one-line edit here instead of a ~20-site sweep.
// ============================================================

// Reserve 2 GB for non-heap (Java metaspace + native); floor at 4 GB.
// Byte-identical to the old inline `Math.max(task.memory.toGiga() - 2, 4)`.
def jvmHeapGb(task) {
    Math.max(task.memory.toGiga() - 2, 4)
}

// Standard GATK java-options string (unchanged content).
def gatkJavaOpts(task) {
    "-Xmx${jvmHeapGb(task)}g -XX:+UseParallelGC"
}

// Scratch-dir shell prelude. Emits the EXACT bash hand-written in
// GENOTYPE_GVCFS / BASE_RECALIBRATOR / MUSE_CALL / etc.: a per-task
// TMPDIR under node-local scratch, created and auto-cleaned on EXIT.
// `tag` keeps each process's descriptive suffix (e.g. 'ggt', 'bqsr',
// 'muse') so work dirs stay self-documenting and collision-free.
// The executed command is identical to the previous inline versions
// (bash ignores the leading indentation of the spliced lines).
def scratchPrelude(meta, String tag) {
    """\
    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_${tag}_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT"""
}


// Bases to hard-trim at read start in the 'markdup' path: the leading run of
// NON-template segments (M=UMI, S=linker, B=barcode) in the fgbio read structure,
// up to the first template (T). A null/absent structure (standard non-UMI library)
// yields 0 - nothing is trimmed. This is the SINGLE source of truth for the trim;
// there is no separate umi_skip_bases to contradict it.
int umiLeadingSkip(String rs) {
    if (!rs) return 0
    def toks = (rs =~ /(?:\d+|\+)[A-Za-z]/).findAll()
    return toks
        .takeWhile { t -> t[-1].toUpperCase() in ['M','S','B'] && t[0..-2] != '+' }
        .sum(0) { t -> t[0..-2].toInteger() }
}


// execution ranks from stageRankMap() + layout-only pre-stage names,
// so a stage-less output (manifest) is numbered without a hard-coded literal
// and stageRankMap() stays the pure --start_at/--stop_at contract.
def stageDir(String stage) {
    def rank = (stageRankMap() + ['manifest': 0])[stage]
    if (rank == null) error "stageDir: unknown stage '${stage}'"
    "${String.format('%02d', rank)}_${stage}"
}

// analysisAware=true segregates by params.analysis (somatic|germline|tumor_only)
// - variant-calling and everything downstream. false = shared preprocessing/QC.
def pubDir(String stage, boolean analysisAware = false) {
    def base = "${params.outdir}/${stageDir(stage)}"
    // RNA outputs always land under rna_fusion/, whatever params.analysis defaults to.
    def analysis = (params.assay == 'rna') ? 'rna_fusion' : params.analysis
    analysisAware ? "${base}/${analysis}" : base
}
// ============================================================
//  GENERATE_META_INFO
//
//  Emits: tuple(meta Map, [fastq_1, fastq_2]) - one per chunk.
//
//  meta fields:
//    Chunk-level (unique per row):  uid, rg
//    Sample-level (shared):         sample_id, patient_id, type,
//                                   sex, batch, assay
//
//  Grouping note: downstream grouping keys on the SCALAR
//  meta.sample_id string (in workflows.nf), NOT on the map.
//  So adding fields to meta later cannot break scatter-gather.
// ============================================================
workflow GENERATE_META_INFO {

    take:
    samplesheet

    main:
    // Strict syntax: constants are LOCAL to the workflow, not top-level.
    def REQUIRED_COLUMNS = [
        'patient_id', 'sample_type', 'sex', 'batch',
        'flowcell', 'run_num', 'lane_num', 'chunk',
        'read1', 'read2'
    ]
    // DT/FFPE = tumour ; DN/BL/FF = normal
    def VALID_SAMPLE_TYPES = ['DT', 'DN']

    ch_meta = channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->

            REQUIRED_COLUMNS.each { col ->
                if (!row.containsKey(col) || row[col] == null || row[col] == '') {
                    error "Samplesheet row missing required column '${col}': ${row}"
                }
            }

            if (!VALID_SAMPLE_TYPES.contains(row.sample_type)) {
                error "Invalid sample_type '${row.sample_type}' for patient ${row.patient_id}. " +
                      "Must be one of: ${VALID_SAMPLE_TYPES.join(', ')}"
            }

            def suffix = (row.sample_type == 'DT')  ? 'T' :
                        (row.sample_type == 'DN')  ? 'N' :
                        (row.sample_type == 'RN')  ? 'R' : 'X'
            def sample_id = "${row.patient_id}_${suffix}"
            def uid       = "${sample_id}_L${row.lane_num}_C${row.chunk}"

            def pu         = "${row.flowcell}.${row.run_num}.${row.lane_num}"
            def read_group = "@RG\\tID:${uid}\\tSM:${sample_id}\\tLB:${row.batch}\\tPL:ILLUMINA\\tPU:${pu}"

            def fastq_1 = file("${params.fastq_root}/${row.read1}", checkIfExists: true)
            def fastq_2 = file("${params.fastq_root}/${row.read2}", checkIfExists: true)

            // sex: pipeline convention 0 = female, 1 = male (Sequenza female flag, sex checks).
            // Accepts 0/1, F/M, female/male and PLINK 2 = female; blank/NA/unknown -> 'NA'.
            def sx  = (row.sex ?: '').toString().trim().toLowerCase()
            def sex = sx in ['0', 'f', 'female', '2'] ? '0' : sx in ['1', 'm', 'male'] ? '1' :
                      sx in ['', 'na', 'nan', 'unknown', 'u', '-9'] ? 'NA' : null
            if (sex == null)
                error "Invalid sex '${row.sex}' for ${sample_id}: use 0/F/female, 1/M/male, or NA"

            def meta = [
                // chunk-level (unique per row) - stripped before grouping
                uid        : uid,
                rg         : read_group,
                // sample-level (shared across all chunks of a sample)
                sample_id  : sample_id,
                patient_id : row.patient_id,
                type       : row.sample_type,
                sex        : sex,
                batch      : row.batch,
                assay      : params.assay
            ]

            return tuple(meta, [fastq_1, fastq_2])
        }

    emit:
    reads = ch_meta
}



// ============================================================
//  GENERATE_RNA_META_INFO
//
//  Emits one tuple(meta, reads) per samplesheet row (chunk). Same samplesheet
//  schema as the DNA arm, but sample_type is 'RN' (or 'RT') -> suffix 'R', and
//  read2 is OPTIONAL (single-end supported; paired-end recommended). Grouping
//  to sample level happens in RNA_WORKFLOW on the scalar meta.sample_id.
// ============================================================
workflow GENERATE_RNA_META_INFO {

    take:
    samplesheet

    main:
    def REQUIRED_COLUMNS = [
        'patient_id', 'sample_type', 'sex', 'batch',
        'flowcell', 'run_num', 'lane_num', 'chunk', 'read1'
    ]
    def VALID_SAMPLE_TYPES = ['RN', 'RT']   // RNA sample; both -> 'R'

    ch_meta = channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->

            REQUIRED_COLUMNS.each { col ->
                if (!row.containsKey(col) || row[col] == null || row[col] == '') {
                    error "RNA samplesheet row missing required column '${col}': ${row}"
                }
            }
            if (!VALID_SAMPLE_TYPES.contains(row.sample_type)) {
                error "Invalid RNA sample_type '${row.sample_type}' for patient ${row.patient_id}. " +
                      "Must be one of: ${VALID_SAMPLE_TYPES.join(', ')}"
            }

            def sample_id = "${row.patient_id}_R"
            def uid       = "${sample_id}_L${row.lane_num}_C${row.chunk}"
            def pu        = "${row.flowcell}.${row.run_num}.${row.lane_num}"
            def read_group = "@RG\\tID:${uid}\\tSM:${sample_id}\\tLB:${row.batch}\\tPL:ILLUMINA\\tPU:${pu}"

            def has_r2  = row.containsKey('read2') && row.read2?.toString()?.trim()
            def fastq_1 = file("${params.fastq_root}/${row.read1}", checkIfExists: true)
            def reads   = [fastq_1]
            if (has_r2) reads << file("${params.fastq_root}/${row.read2}", checkIfExists: true)

            def meta = [
                uid         : uid,
                rg          : read_group,
                sample_id   : sample_id,
                patient_id  : row.patient_id,
                type        : row.sample_type,
                sex         : row.sex,
                batch       : row.batch,
                assay       : params.assay,
                single_end  : !has_r2
            ]
            return tuple(meta, reads)
        }

    emit:
    reads = ch_meta
}

