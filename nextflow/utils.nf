#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
//  utils.nf  -  Channel Factory + Helper Sub-Workflows
// ============================================================



// ============================================================
//  SHARED CONSTANTS  -  THE single source of truth
//
//  Previously stageRankMap()'s contents were declared TWICE:
//  stageRank() in main.nf AND a local STAGE_RANK map inside
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
// rna_strandedness as a lower-case word. A params YAML reads an unquoted yes/no as a
// boolean (YAML 1.1), so true/false are taken as Arriba's 'yes'/'no'.
def rnaStrandedness() {
    def v = params.rna_strandedness
    if (v instanceof Boolean) return v ? 'yes' : 'no'
    return (v == null || v.toString().trim() == '') ? 'auto' : v.toString().trim().toLowerCase()
}
// Somatic callers: with somatic_consensus a caller runs when its use_<caller> is on; without it
// only somatic_caller runs. somaticCallSet() names the final call set: consensus | single_<caller>.
def somaticCaller()            { (params.somatic_caller ?: 'mutect2').toString().toLowerCase() }
def somaticCallerOn(String c)  { params.somatic_consensus ? (params["use_${c}".toString()] as boolean) : somaticCaller() == c }
def somaticCallSet()           { params.somatic_consensus ? 'consensus' : "single_${somaticCaller()}" }
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

// Reserve 2 GB for non-heap (Java metaspace + native). On small allocations (< 8 GB) keep
// the heap at 75% of the task instead: the old 4 GB floor gave a heap as large as, or larger
// than, a 4 GB task, which the scheduler then kills. Unchanged for tasks of 8 GB or more.
def jvmHeapGb(task) {
    def gb = task.memory.toGiga()
    Math.max((gb - 2) as int, Math.max(1, (gb * 0.75d) as int))
}

// Standard GATK java-options string (unchanged content).
def gatkJavaOpts(task) {
    "-Xmx${jvmHeapGb(task)}g -XX:+UseParallelGC"
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

// Samplesheet FASTQ path: absolute paths are used as given, relative ones are
// resolved against params.fastq_root.
def fastqPath(String p) {
    def s = p.trim()
    if (s.startsWith('/'))
        return file(s, checkIfExists: true)
    if (!params.fastq_root)
        error "ERROR: samplesheet FASTQ '${s}' is a relative path but params.fastq_root is not set"
    return file("${params.fastq_root}/${s}", checkIfExists: true)
}

// Natural order for shard names ('shard_2' < 'shard_10'); plain string order
// puts 'shard_10' first when the names are not zero-padded.
def naturalCompare(String a, String b) {
    def ta = a.findAll(/\d+|\D+/)
    def tb = b.findAll(/\d+|\D+/)
    def n  = Math.min(ta.size(), tb.size())
    def first_diff = (0..<n).collect { i ->
        def x = ta[i]
        def y = tb[i]
        (x.isNumber() && y.isNumber()) ? (x as BigInteger) <=> (y as BigInteger) : x <=> y
    }.find { c -> c != 0 }
    return first_diff != null ? first_diff : (ta.size() <=> tb.size())
}

// analysisAware=true segregates by params.analysis (somatic|germline|tumor_only)
// - variant-calling and everything downstream. false = shared preprocessing/QC.
def pubDir(String stage, boolean analysisAware = false) {
    def base = "${params.outdir}/${stageDir(stage)}"
    // RNA outputs always land under rna_fusion/, whatever params.analysis defaults to.
    def analysis = (params.assay == 'rna') ? 'rna_fusion' : params.analysis
    analysisAware ? "${base}/${analysis}" : base
}
// Optional samplesheet column 'sample': distinguishes several samples of one
// patient and type (multi-region / longitudinal tumours, repeat RNA). Absent or
// blank -> the classic <patient>_<T|N|R> ID, so existing samplesheets are unchanged.
// No underscore allowed: the swap check recovers patient_id by stripping the last
// '_' token of sample_id.
def sampleLabel(Map row) {
    def lab = (row.containsKey('sample') && row.sample != null) ? row.sample.toString().trim() : ''
    if (lab && !(lab ==~ /[A-Za-z0-9-]+/))
        error "Invalid sample '${row.sample}' for patient ${row.patient_id}: use letters, digits or '-' only"
    return lab
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
            def label     = sampleLabel(row)
            def sample_id = label ? "${row.patient_id}_${label}" : "${row.patient_id}_${suffix}"
            def uid       = "${sample_id}_L${row.lane_num}_C${row.chunk}"

            def pu         = "${row.flowcell}.${row.run_num}.${row.lane_num}"
            def read_group = "@RG\\tID:${uid}\\tSM:${sample_id}\\tLB:${row.batch}\\tPL:ILLUMINA\\tPU:${pu}"

            def fastq_1 = fastqPath(row.read1)
            def fastq_2 = fastqPath(row.read2)

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

            def label     = sampleLabel(row)
            def sample_id = label ? "${row.patient_id}_${label}" : "${row.patient_id}_R"
            def uid       = "${sample_id}_L${row.lane_num}_C${row.chunk}"
            def pu        = "${row.flowcell}.${row.run_num}.${row.lane_num}"
            def read_group = "@RG\\tID:${uid}\\tSM:${sample_id}\\tLB:${row.batch}\\tPL:ILLUMINA\\tPU:${pu}"

            def has_r2  = row.containsKey('read2') && row.read2?.toString()?.trim()
            def fastq_1 = fastqPath(row.read1)
            def reads   = [fastq_1]
            if (has_r2) reads << fastqPath(row.read2)

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

