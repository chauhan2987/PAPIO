#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Shared JVM-options helper.
include { gatkJavaOpts; pubDir; umiLeadingSkip} from './utils'
// ============================================================
//  components_common.nf  -  Processes shared by the DNA arms (and some by RNA)
// ============================================================


// ============================================================
//  STAGE 0 setup: MANIFEST_ROW - run manifest - one audit row per sample
//
//  Records each sample's identity, T/N type, sex, assay, batch, input
//  file, size, and (optionally) a checksum, consolidated into one
//  run_manifest.tsv. Checksum level is a param: 'none' (size only, light,
//  default), 'header' (BAM header md5 - light, provenance), or 'full'
//  (whole-file md5 - audit-grade chain-of-custody, expensive on large BAMs).
// ============================================================

process MANIFEST_ROW {
    tag { "${meta.uid}" }
    label 'process_low'

    input:
    tuple val(meta), path(input_file)

    output:
    path "${meta.uid}_manifest_row.tsv", emit: row

    script:
    def cmode = params.manifest_checksum ?: 'none'

    """
    set -euo pipefail

    sz=\$(stat -c%s "${input_file}" 2>/dev/null || echo NA)

    case "${cmode}" in
        full)
            cks=\$(md5sum "${input_file}" | awk '{print \$1}')
            cmode_out="md5:full"
            ;;
        header)
            if samtools view -H "${input_file}" >/dev/null 2>&1; then
                cks=\$(samtools view -H "${input_file}" | md5sum | awk '{print \$1}')
                cmode_out="md5:header"
            else
                cks="NA"; cmode_out="header_unavailable"
            fi
            ;;
        *)
            cks="NA"; cmode_out="none"
            ;;
    esac

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        "${meta.uid}" "${meta.patient_id}" "${meta.type}" \
        "${meta.sex}" "${meta.assay}" "${meta.batch}" \
        "${input_file}" "\$sz" "\$cks" "\$cmode_out" \
        > ${meta.uid}_manifest_row.tsv
    """
}


// ============================================================
//  STAGE 0 setup: WRITE_MANIFEST - run manifest - write the cohort table
// ============================================================

process WRITE_MANIFEST {
    label 'process_low'

    // -- 00 manifest -----------------------------------------
    publishDir { "${pubDir('manifest', true)}" }, mode: 'copy'


    input:
    path rows

    output:
    path "run_manifest.tsv", emit: manifest

    script:
    """
    set -euo pipefail
    {
        printf 'uid\\tpatient_id\\ttype\\tsex\\tassay\\tbatch\\tinput_file\\tsize_bytes\\tchecksum\\tchecksum_mode\\n'
        cat ${rows}
    } > run_manifest.tsv
    """
}


// ============================================================
//  SHARED: MULTIQC - MultiQC report (called by several stages)
//
//   Multisatge: MultiQC Aggregation
// ============================================================

process MULTIQC {
    tag { "${label}" }
    label 'process_low'

    // MultiQC aggregates ACROSS all samples of a run - it is not a per-stage
    // -- 0x multiqc  -----------------------------------------
    publishDir { "${params.outdir}/0x_multiqc/${params.assay == 'rna' ? 'rna_fusion' : params.analysis}" }, mode: 'copy'

    input:
    tuple val(label), path(files, arity: '1..*')

    output:
    path "${label}_multiqc_report.html", emit: report

    script:
    """
    set -euo pipefail
    multiqc ${files.join(' ')} --filename ${label}_multiqc_report.html --force
    """
}


// ============================================================
//  SHARED: CRAM_CONVERT - BAM to CRAM archive (markdup, bqsr and RNA mapping)
//
//  Multisatge CRAM CORVERT: CRAM CONVERSION
// ============================================================

process CRAM_CONVERT {
    tag { "${label}" }

    label 'process_medium'

    // -- 0x cram save ----------------------------------------
    publishDir { "${params.outdir}/0x_crams/${label}" }, mode: 'copy', pattern: '*.cram*'


    input:
    val(label)
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta),
          path("${meta.sample_id}.cram"),
          path("${meta.sample_id}.cram.crai"), emit: cram_crai_pair

    script:
    """
    set -euo pipefail
    samtools view -@ ${task.cpus} -C -T ${ref_fasta} -o ${meta.sample_id}.cram ${bam}
    samtools index -@ ${task.cpus} ${meta.sample_id}.cram
    """
}


// ============================================================
//  STAGE 1 mapping: PREPROCESS_AND_ALIGN_DNA - fastp trimming, bwa-mem2 alignment, fixmate (dedup_mode markdup)
//
//  Architecture: fastp -> bwa-mem2 -> samtools fixmate
//
//  OPTIMIZATION: Thread budget explicitly partitioned so total
//  threads never exceeds task.cpus (prevents LSF over-subscription).
//  fastp(2) + fixmate(2) + bwa(remainder) = task.cpus exactly.
// ============================================================

process PREPROCESS_AND_ALIGN_DNA {
    tag { "${meta.uid}" }
    label 'process_high_cpu'

    // scratch = true set globally in nextflow.config
    maxRetries    2
    // memory        { task.attempt == 1 ? '80 GB' : task.attempt == 2 ? '85 GB' : '95 GB' }

    input:
    tuple val(meta), path(reads)
    path  ref_fasta                 // GRCh38 primary assembly FASTA
    path  ref_index                 // All bwa-mem2 index files (.0123, .bwt.2bit.64, etc.)

    output:
    tuple val(meta), path("${meta.uid}.fixmate.bam"), emit: chunk_bams
    path "${meta.uid}_fastp.json",                    emit: json_report
    path "${meta.uid}_fastp.html",                    emit: html_report

    script:
    def trim_flag       = params.skip_trimming ? '--disable_adapter_trimming' : '--detect_adapter_for_pe'
    def fastp_threads   = 2
    def fixmate_threads = 2
    def bwa_threads     = Math.max(task.cpus - (fastp_threads + fixmate_threads), 1)

    // Hard-trim the inline UMI+linker prefix ONLY when this library actually
    // carries an inline UMI (umi_read_structure explicitly set). A plain
    // non-UMI markdup library leaves umi_read_structure null -> no trim.
    // R2 may carry a different structure (params.umi_read_structure_r2; null = same as R1).
    def rs_r2      = params.umi_read_structure_r2 ?: params.umi_read_structure
    def skip_n     = params.umi_read_structure ? umiLeadingSkip(params.umi_read_structure) : 0
    def skip_n2    = rs_r2 ? umiLeadingSkip(rs_r2) : 0
    def front_trim = ((skip_n > 0) ? "--trim_front1 ${skip_n} " : '') + ((skip_n2 > 0) ? "--trim_front2 ${skip_n2}" : '')

    """
    set -euo pipefail

    fastp \\
        --in1   ${reads[0]} --in2   ${reads[1]} \\
        --stdout --thread ${fastp_threads} \\
        ${trim_flag} --trim_poly_g ${front_trim} \\
        --json  ${meta.uid}_fastp.json \\
        --html  ${meta.uid}_fastp.html \\
    | \\
    bwa-mem2 mem \\
        -t ${bwa_threads} -p -K 100000000 -Y \\
        -R "${meta.rg}" ${ref_fasta} - \\
    | \\
    samtools fixmate -@ ${fixmate_threads} -m -u - ${meta.uid}.fixmate.bam
    """
}


// ============================================================
//  STAGE 1 mapping: ALIGN_DNA_UMI - UMI extraction (fgbio) and bwa-mem2 alignment (dedup_mode umi_*)
//
//  fgbio FastqToBam reads the 3M2S+T structure on BOTH mates:
//    3M = 3 bp molecular barcode  - RX tag
//    2S = 2 bp linker             - discarded
//    +T = remaining template      - genomic sequence kept
//  The dual barcodes are stored as "AAA-BBB" in RX. We then map with
//  bwa-mem2 and ZipperBams re-attaches RX + the read group onto the
//  aligned records (bwa strips tags otherwise).
// ============================================================

process ALIGN_DNA_UMI {
    tag { "${meta.uid}" }
    label 'process_high_cpu'
    maxRetries    2

    input:
    tuple val(meta), path(reads)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  ref_index

    output:
    tuple val(meta), path("${meta.uid}.umi.mapped.bam"), emit: chunk_bams

    script:
    // Heaps and sort buffer scale with the task's memory (capped at the tuned values) so a
    // smaller allocation does not over-commit: FastqToBam runs alone; ZipperBams, bwa-mem2
    // (index-bound) and the 2-thread sort run together.
    def mem_gb       = task.memory ? task.memory.toGiga() : 48
    def fastq_heap   = Math.max(2, Math.min(8, (mem_gb * 0.5) as int))
    def zipper_heap  = Math.max(1, Math.min(4, (mem_gb * 0.08) as int))
    def sort_mem_mb  = Math.max(256, Math.min(4096, (mem_gb * 1024 * 0.1 / 2) as int))

    // Dynamically calculate threads to avoid over-subscribing the cluster node
    def bwa_threads  = Math.max((task.cpus - 4).toInteger(), 4)

    """
    set -euo pipefail

    tmp_dir="\$(pwd)/fgbio_tmp"
    mkdir -p "\${tmp_dir}"
    trap 'rm -rf "\${tmp_dir}"' EXIT

    # STEP 1: Generate the unmapped BAM containing UMI structural molecular tags
    fgbio --tmp-dir "\${tmp_dir}" -Xmx${fastq_heap}g FastqToBam \\
        --input ${reads[0]} ${reads[1]} \\
        --read-structures ${params.umi_read_structure} ${params.umi_read_structure_r2 ?: params.umi_read_structure} \\
        --sample ${meta.sample_id} \\
        --library ${meta.batch} \\
        --read-group-id ${meta.uid} \\
        --platform ILLUMINA \\
        --output ./chunk_unmapped.bam

    # STEP 2: Stream through bwa-mem2 -> ZipperBams -> samtools sort
    # This guarantees the file saved to disk is 100% coordinate sorted with valid headers.
    samtools fastq ./chunk_unmapped.bam \\
    | bwa-mem2 mem -t ${bwa_threads} -p -K 100000000 -Y ${ref_fasta} - \\
    | fgbio --tmp-dir "\${tmp_dir}" --compression 0 -Xmx${zipper_heap}g ZipperBams \\
        --unmapped ./chunk_unmapped.bam \\
        --ref ${ref_fasta} \\
    | samtools sort -@ 2 -m ${sort_mem_mb}M -T "\${tmp_dir}"/sam_sort_tmp -o ${meta.uid}.umi.mapped.bam -

    # Clean up intermediate unmapped data
    rm -f ./chunk_unmapped.bam
    rm -rf "\${tmp_dir}"
    """
}


// ============================================================
//  STAGE 2 markdup: MERGE_AND_MARK_DUPLICATES - merge chunks, sort and mark duplicates (samtools markdup)
// ============================================================

process MERGE_AND_MARK_DUPLICATES {
    tag { "${meta.sample_id}" }

    label 'process_markdup'

    maxRetries    2

    // -- 02 markdup (shared; dedup modes mutually exclusive) ------------
    publishDir { pubDir('markdup') }, mode: 'copy', pattern: '*_dup_metrics.json'

    input:
    tuple val(meta), path(chunk_bams)

    output:
    tuple val(meta),
          path("${meta.sample_id}_markdup.bam"),
          path("${meta.sample_id}_markdup.bai"), emit: bam_bai_pair
    // JSON, not text: MultiQC 1.35's samtools/markdup_txt search pattern is the
    // literal substring "^COMMAND:", which never occurs in samtools' text stats,
    // so text stats are silently undetected. markdup_json matches correctly.
    path  "${meta.sample_id}_dup_metrics.json",      emit: metrics

    script:
    // NOTE: samtools sort -m is PER-THREAD. Total sort buffer = sort_mem * task.cpus.
    // The 0.8 factor caps that at ~80% of the request, leaving headroom for the
    // merge + markdup stages running concurrently in the same pipe.
    def sort_mem = "${(task.memory.toGiga() * 0.8 / task.cpus).toInteger()}G"

    """
    set -euo pipefail

    samtools merge -@ ${task.cpus} -u - ${chunk_bams} \\
    | \\
    samtools sort -@ ${task.cpus} -m ${sort_mem} -u - \\
    | \\
    samtools markdup \\
        -@ ${task.cpus} -d 2500 -S --json \\
        -f ${meta.sample_id}_dup_metrics.json \\
        - ${meta.sample_id}_markdup.bam

    samtools index -@ ${task.cpus} ${meta.sample_id}_markdup.bam ${meta.sample_id}_markdup.bai
    """
}


// ============================================================
//  STAGE 2 markdup: MERGE_AND_CONSENSUS - merge chunks and call UMI consensus reads (dedup_mode umi_consensus)
//
//  fgbio single-strand consensus, the ideal FFPE path:
//    GroupReadsByUmi (adjacency)  - families by UMI + position
//    CallMolecularConsensusReads  - one error-corrected read/family
//    re-align consensus reads     - bwa-mem2 + ZipperBams
//    FilterConsensusReads         - drop low-support / high-error consensus
//  Consensus REPLACES duplicate marking - families are already collapsed.
// ============================================================

process MERGE_AND_CONSENSUS {
    tag { "${meta.sample_id}" }
    label 'process_high_memory' // 150 GB RAM / 14 CPU

    maxRetries    2

    // -- 02 markdup (shared; dedup modes mutually exclusive) ------------
    publishDir { pubDir('markdup') }, mode: 'copy' , pattern: '*_family_sizes.txt'

    input:
    tuple val(meta), path(chunk_bams)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  ref_index

    output:
    tuple val(meta), path("${meta.sample_id}_consensus.bam"), path("${meta.sample_id}_consensus.bai"), emit: bam_bai_pair
    path  "${meta.sample_id}_family_sizes.txt", emit: metrics

    script:
    def bwa_threads      = Math.max(task.cpus - 2, 1)

    // Tuned for 150 GB; every buffer scales down with a smaller allocation (factor f) so
    // the concurrent JVMs + sort never exceed task.memory.
    def mem_gb           = task.memory ? task.memory.toGiga() : 150
    def f                = Math.min(1.0d, mem_gb / 150.0d)
    def ncpu             = (task.cpus ?: 1) as int
    def phase1_sort_mem  = "${Math.max(256, Math.min(5120, (mem_gb * 1024 * 0.5 / ncpu) as int))}M"   // sort runs alone
    def group_heap       = Math.max(2, (24 * f) as int)
    def call_heap        = Math.max(2, (16 * f) as int)
    def filter_heap      = Math.max(2, (16 * f) as int)
    def zipper_heap      = Math.max(1, (8 * f) as int)
    def phase4_sort_mem  = "${Math.max(128, Math.min(1500, (mem_gb * 1024 * 0.15 / ncpu) as int))}M"

    // Consensus parameters
    def min_reads    = params.consensus_min_reads ?: 1
    def err_pre      = params.consensus_error_rate_pre  ?: 45
    def err_post     = params.consensus_error_rate_post ?: 40
    def max_read_err = params.consensus_max_read_error ?: 0.05
    def min_base_q   = params.consensus_min_base_quality ?: 20

    """
    set -euo pipefail

    # task-local scratch (JVM and samtools spill files)
    tmp_dir="./mer_con_tmp"
    mkdir -p "\$tmp_dir"

    export TMPDIR="\$tmp_dir"

    # intermediates are removed on any exit
    trap 'rm -rf "\$tmp_dir" ./sorted_template.bam ./grouped_transient.bam ./unmapped_consensus.bam' EXIT

    # =========================================================================
    # STEP 1: INITIAL TEMPLATE COORDINATE SORT (STANDALONE)
    # =========================================================================
    echo "[STEP 1] Running TemplateCoordinate sort..."
    samtools merge -@ ${task.cpus} -u - ${chunk_bams} \\
    | samtools sort -@ ${task.cpus} --template-coordinate -m ${phase1_sort_mem} -T "\$tmp_dir/sam_spill" -O BAM -o ./sorted_template.bam -

    # =========================================================================
    # STEP 2: GROUP READS BY UMI (ISOLATED PHASE TO PREVENT PIPE DEADLOCKS)
    # =========================================================================
    echo "[STEP 2] Running GroupReadsByUmi..."
    fgbio --tmp-dir "\$tmp_dir" --compression 0 -Xmx${group_heap}g GroupReadsByUmi \\
        --input ./sorted_template.bam \\
        --strategy adjacency \\
        --raw-tag RX --assign-tag MI \\
        --min-map-q 20 \\
        --family-size-histogram ${meta.sample_id}_family_sizes.txt \\
        --output ./grouped_transient.bam

    # free disk before the next phase
    rm -f ./sorted_template.bam

    # =========================================================================
    # STEP 3: STREAMED CONSENSUS CALLING & FILTERING
    # =========================================================================
    echo "[STEP 3] Running Call -> Filter stream..."
    fgbio --tmp-dir "\$tmp_dir" --compression 0 -Xmx${call_heap}g CallMolecularConsensusReads \\
        --input ./grouped_transient.bam \\
        --min-reads ${min_reads} \\
        --error-rate-pre-umi ${err_pre} \\
        --error-rate-post-umi ${err_post} \\
        --output /dev/stdout \\
    | fgbio --tmp-dir "\$tmp_dir" --compression 0 -Xmx${filter_heap}g FilterConsensusReads \\
        --input /dev/stdin \\
        --ref ${ref_fasta} \\
        --min-reads ${min_reads} \\
        --max-read-error-rate ${max_read_err} \\
        --min-base-quality ${min_base_q} \\
        --output ./unmapped_consensus.bam

    rm -f ./grouped_transient.bam

    # =========================================================================
    # STEP 4: CONCURRENT ALIGNMENT & FINAL COORDINATE SORT
    # =========================================================================
    echo "[STEP 4] Re-align consensus + ZipperBams + coordinate sort..."
    samtools fastq ./unmapped_consensus.bam \\
    | bwa-mem2 mem -t ${bwa_threads} -p -K 100000000 -Y ${ref_fasta} - \\
    | fgbio --tmp-dir "\$tmp_dir" --compression 1 -Xmx${zipper_heap}g ZipperBams \\
        --unmapped ./unmapped_consensus.bam \\
        --ref ${ref_fasta} \\
        --tags-to-reverse Consensus \\
        --tags-to-revcomp Consensus \\
        --output /dev/stdout \\
    | samtools sort -@ ${task.cpus} -m ${phase4_sort_mem} -T "\$tmp_dir/final_spill" -o ${meta.sample_id}_consensus.bam -

    samtools index -@ ${task.cpus} ${meta.sample_id}_consensus.bam ${meta.sample_id}_consensus.bai
    """
}


// ============================================================
//  STAGE 2 markdup: MERGE_AND_MARKDUP_UMI - merge chunks and UMI-aware duplicate marking (dedup_mode umi_markdup)
//
//  Lighter than consensus: uses the UMI to dedup more accurately than
//  coordinates alone (Picard MarkDuplicates BARCODE_TAG=RX), but does
//  NOT build consensus / error-correct. Better than plain markdup,
//  cheaper than consensus. Good when family sizes are too small for
//  consensus to pay off (low coverage).
// ============================================================

process MERGE_AND_MARKDUP_UMI {
    tag { "${meta.sample_id}" }
    label 'process_markdup'

    maxRetries    2

    // -- 02 markdup (shared; dedup modes mutually exclusive) ------------
    publishDir { pubDir('markdup') }, mode: 'copy', pattern: '*_dup_metrics.txt'

    input:
    tuple val(meta), path(chunk_bams)
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(meta), path("${meta.sample_id}_markdup.bam"), path("${meta.sample_id}_markdup.bai"), emit: bam_bai_pair
    path  "${meta.sample_id}_dup_metrics.txt", emit: metrics

    script:
    def parallel_gc_threads = Math.max((task.cpus / 2).toInteger(), 2)
    // Reserve 15% (min 10G) as OS/JVM-overhead cushion so LSF doesn't evict.
    // Non-heap reserve: 15% (at least 10 GB on large nodes, but never more than a quarter of
    // the task, so a small allocation still leaves a positive heap).
    def task_gb       = task.memory.toGiga()
    def system_buffer = Math.max(task_gb * 0.15d, Math.min(10.0d, task_gb * 0.25d))
    def mem_gb        = Math.max((task_gb - system_buffer) as int, 1)
    // def java_opts     = "-Xmx${mem_gb}g -XX:ParallelGCThreads=${parallel_gc_threads} -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_read_samtools=true"
    def java_opts = "-Xmx${mem_gb}g -XX:ParallelGCThreads=${parallel_gc_threads}"
    def bam_list      = chunk_bams instanceof List ? chunk_bams : [chunk_bams]
    def input_list    = bam_list.collect { b -> "-I ${b}" }.join(' ')
    """
    set -euo pipefail

    tmp_dir="\$(pwd)/mdumi_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" MarkDuplicates \\
        ${input_list} \\
        -O ${meta.sample_id}_markdup.bam \\
        -M ${meta.sample_id}_dup_metrics.txt \\
        --BARCODE_TAG RX \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        --SORTING_COLLECTION_SIZE_RATIO 0.5 \\
        --TMP_DIR "\$tmp_dir"

    samtools index -@ ${task.cpus} ${meta.sample_id}_markdup.bam ${meta.sample_id}_markdup.bai
    """
}


// ============================================================
//  STAGE 3 metrics: MOSDEPTH - coverage - mosdepth
// ============================================================

process MOSDEPTH {
    tag { "${meta.sample_id}" }
    label 'process_qc'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/mosdepth" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  interval_bed
    path  coverage_report_bed, stageAs: 'report_bed/*'   // params.coverage_report_bed (reportable genes/exons, col4 = gene) or []

    output:
    tuple val(meta), path("${meta.sample_id}.mosdepth.summary.txt"),     emit: summary
    tuple val(meta), path("${meta.sample_id}.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(meta), path("${meta.sample_id}.mosdepth.region.dist.txt"), optional: true, emit: region_dist
    tuple val(meta), path("${meta.sample_id}.regions.bed.gz"),           optional: true, emit: regions
    tuple val(meta), path("${meta.sample_id}.thresholds.bed.gz"),        optional: true, emit: thresholds
    tuple val(meta), path("${meta.sample_id}.gene_coverage.tsv"),        optional: true, emit: gene_coverage
    tuple val(meta), path("${meta.sample_id}.quantized.bed.gz"),         emit: quantized   // depth classes (callable territory for TMB)

    script:
    // WES: per-target regions (--by capture BED) enable --thresholds + --quantize.
    // WGS: no capture regions. mosdepth --thresholds REQUIRES --by, and --fast-mode
    //      is incompatible with --quantize - so WGS uses windowed --by instead, which
    //      keeps thresholds/quantize working and gives genome coverage in windows.
    def is_wgs = (meta.assay == 'wgs')
    def by_flag = is_wgs
        ? "--by ${params.mosdepth_wgs_window ?: 500}"
        : ((meta.assay == 'wes' && interval_bed) ? "--by ${interval_bed}" : '')
    // --thresholds is a TARGET-capture metric (% of target bases >= X). Meaningful
    // for WES only. On WGS windows it measures "% of arbitrary windows crossing X",
    // which is not comparable to the WES gate - so WGS coverage QC uses summary +
    // global.dist instead. --thresholds also REQUIRES --by, hence the WES+BED guard.
    def thresholds_flag = (!is_wgs && by_flag) ? "--thresholds 1,10,20,30,50,100" : ''
    // Per-gene coverage of reportable genes (WES + WGS): a negative result is only
    // reportable where the gene is covered. MAPQ-filtered (usable) depth.
    def gene_mapq = (params.coverage_report_min_mapq != null) ? params.coverage_report_min_mapq : 20
    def gene_cmd  = !coverage_report_bed ? '' :
        "gene_coverage.sh --sample ${meta.sample_id} --bam ${bam} --bed ${coverage_report_bed} --mapq ${gene_mapq} --threads ${task.cpus}"
    """
    set -euo pipefail
    # --no-per-base: the per-base track (several GB for WGS) is never used; --quantize writes
    # the merged depth classes instead
    mosdepth \\
        --threads ${task.cpus} --no-per-base \\
        --quantize 0:1:10:20:30:50:100: \\
        ${thresholds_flag} \\
        ${by_flag} \\
        ${meta.sample_id} ${bam}
    ${gene_cmd}
    """
}


// ============================================================
//  STAGE 3 metrics: INFER_SEX - sex inference from X/autosome coverage
// ============================================================

process INFER_SEX {
    tag { meta.sample_id }
    label 'process_low'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/sex" }, mode: 'copy'

    input:
    tuple val(meta), path(mosdepth_summary)

    output:
    tuple val(meta), path("${meta.sample_id}_sex.tsv"), emit: sex

    script:
    def female_min = params.sex_female_min_xratio ?: 0.80
    def male_max   = params.sex_male_max_xratio   ?: 0.65
    def female_max_y = (params.sex_female_max_yratio != null) ? params.sex_female_max_yratio : 0.10
    def male_min_y   = (params.sex_male_min_yratio   != null) ? params.sex_male_min_yratio   : 0.20
    """
    infer_sex.sh --sample ${meta.sample_id} --summary ${mosdepth_summary} --expected "${meta.sex}" \\
        --type ${meta.type} --female-min ${female_min} --male-max ${male_max} \\
        --female-max-y ${female_max_y} --male-min-y ${male_min_y}
    """
}


// ============================================================
//  STAGE 3 metrics: SAMTOOLS_STATS - alignment statistics - samtools stats
//
//  Captures: mapping rate, error rate, average read length,
//  insert size distribution, mismatch rate.
//  mosdepth captures depth; samtools stats captures alignment
//  quality. Both are needed for a complete QC picture.
//  MultiQC parses samtools stats natively.
// ============================================================

process SAMTOOLS_STATS {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/samtools" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta), path("${meta.sample_id}.stats"), emit: stats

    script:
    """
    set -euo pipefail
    samtools stats \\
        --threads   ${task.cpus} \\
        --reference ${ref_fasta} \\
        ${bam} > ${meta.sample_id}.stats
    """
}


// ============================================================
//  STAGE 3 metrics: METRICS - Picard CollectMultipleMetrics and CollectHsMetrics
// ============================================================

process METRICS {
    tag { "${meta.sample_id}" }

    label 'process_qc'

    maxRetries    2

    //   Picard CollectHsMetrics memory-maps the BAM, so its VIRTUAL size scales with BAM depth
    memory {
        def in_gb = (bam.size() ?: 0) / 1e9
        (in_gb > params.deep_input_gb ? 150 : 16) * 1.GB * task.attempt
    }
    time {
        // Deep tier: CollectMultipleMetrics (+GcBias +SequencingArtifact) and
        // CollectHsMetrics run serially; a deep WES tumor exceeded 150 min
        // (killed ~40 min into HsMetrics, chr14). 300 min sizes attempt 1 to the
        // real job; the time-derived queue routes it to medium.
        def in_gb = (bam.size() ?: 0) / 1e9
        (in_gb > params.deep_input_gb ? 300 : 60) * 1.min * task.attempt
    }

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/picard" }, mode: 'copy', pattern: '*{txt,metrics}'

    input:
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  bait_interval_list
    path  target_interval_list

    output:
    tuple val(meta), path("${meta.sample_id}.alignment_summary_metrics"), emit: alignment
    tuple val(meta), path("${meta.sample_id}.insert_size_metrics"),        emit: insert_size
    tuple val(meta), path("${meta.sample_id}.hs_metrics"), optional: true, emit: hs_metrics
    // GC bias + sequencing-artifact (FFPE deamination / OxoG) metrics; params.picard_extra_metrics
    tuple val(meta), path("${meta.sample_id}.{gc_bias.summary_metrics,pre_adapter_summary_metrics,bait_bias_summary_metrics}"),
          optional: true, emit: extra_metrics

    script:
    def java_opts = gatkJavaOpts(task)
    def run_hs    = (meta.assay == 'wes') && bait_interval_list && target_interval_list
    // Same single BAM pass as the existing programs (no extra traversal).
    def extra_programs = params.picard_extra_metrics
        ? "--PROGRAM CollectGcBiasMetrics --PROGRAM CollectSequencingArtifactMetrics" : ''


    """
    set -euo pipefail

    tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_metrics_tmp.XXXXXX")
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" CollectMultipleMetrics \\
        -I ${bam} -R ${ref_fasta} -O ${meta.sample_id} \\
        --PROGRAM CollectAlignmentSummaryMetrics \\
        --PROGRAM CollectInsertSizeMetrics \\
        ${extra_programs} \\
        --VALIDATION_STRINGENCY SILENT --QUIET true \\
        --TMP_DIR "\$tmp_dir"
    """ +
    (run_hs ? """

    gatk --java-options "${java_opts}" CollectHsMetrics \\
        -I ${bam} -R ${ref_fasta} \\
        -BAIT_INTERVALS ${bait_interval_list} -TARGET_INTERVALS ${target_interval_list} \\
        -O ${meta.sample_id}.hs_metrics \\
        --VALIDATION_STRINGENCY SILENT --TMP_DIR "\$tmp_dir"
    """ : "")
}


// ============================================================
//  STAGE 3 metrics: VERIFYBAMID2 - contamination of normals - VerifyBamID2
// ============================================================

process VERIFYBAMID2 {
    tag { "${meta.sample_id}" }

    label 'process_medium'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/verifybamid" }, mode: 'copy', pattern: '*.selfSM'

    input:
    tuple val(meta), path(bam), path(bai)
    path  verifybamid_db
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta), path("${meta.sample_id}.selfSM"), emit: contamination

    script:
    """
    set -euo pipefail
    VerifyBamID \\
        --SVDPrefix ${verifybamid_db}/${params.verifybamid_svd} \\
        --Reference ${ref_fasta} \\
        --BamFile   ${bam} \\
        --Output    ${meta.sample_id} \\
        --NumThread ${task.cpus} \\
        --DisableSanityCheck
    """
}


// ============================================================
//  STAGE 3 metrics: SOMALIER_EXTRACT - sample-identity sites - somalier extract
//
//  somalier extract fingerprints each BAM/CRAM at known sites;
//  somalier relate then checks that tumor and normal from the same
//  patient_id actually match (relatedness ~1) and flags swaps/mislabels.
//  Critical for cohorts of many tumor/normal pairs.
// ============================================================

process SOMALIER_EXTRACT {
    tag { "${meta.sample_id}" }
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai
    path  somalier_sites

    output:
    path "extracted/${meta.sample_id}.somalier", emit: counts
    tuple val(meta), path("extracted/${meta.sample_id}.somalier"), emit: counts_meta

    script:
    """
    set -euo pipefail
    mkdir -p extracted
    somalier extract \\
        --sites ${somalier_sites} \\
        --fasta ${ref_fasta} \\
        --sample-prefix "" \\
        -d extracted \\
        ${bam}
    """
}


// ============================================================
//  STAGE 3 metrics: SOMALIER_RELATE - relatedness - somalier relate
// ============================================================

process SOMALIER_RELATE {
    tag   "cohort"
    label 'process_low'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/somalier" }, mode: 'copy'

    input:
    path somalier_counts   // all per-sample .somalier files

    output:
    path "somalier.html",          emit: report
    path "somalier.pairs.tsv",     emit: pairs
    path "somalier.samples.tsv",   emit: samples

    script:
    """
    set -euo pipefail
    somalier relate \\
        --output-prefix somalier \\
        ${somalier_counts}
    """
}


// ============================================================
//  STAGE 3 metrics: SOMALIER_ANCESTRY - genetic ancestry - somalier ancestry
//
//  Normals only, projected onto the 1000 Genomes somalier profiles.
//  Needs the same somalier sites file as SOMALIER_EXTRACT.
// ============================================================

process SOMALIER_ANCESTRY {
    tag   "cohort"
    label 'process_low'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/somalier" }, mode: 'copy'

    input:
    path normal_counts                  // <pid>_N.somalier files
    path labels                         // ancestry-labels-1kg.tsv
    path ref_dir                        // 1kg-somalier/ (one .somalier per 1000 Genomes sample)

    output:
    path "somalier-ancestry.somalier-ancestry.tsv",  emit: ancestry
    path "somalier-ancestry.somalier-ancestry.html", emit: report

    script:
    """
    set -euo pipefail
    somalier ancestry \\
        --labels ${labels} \\
        --output-prefix somalier-ancestry \\
        ${ref_dir}/*.somalier ++ ${normal_counts}
    """
}


// ============================================================
//  STAGE 3 metrics: CHECK_SAMPLE_SWAPS - tumor-normal pairing and sample-swap verdicts
// ============================================================

process CHECK_SAMPLE_SWAPS {
    tag { "swap_check" }
    label 'process_low'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/somalier" }, mode: 'copy'

    input:
    path pairs_tsv

    output:
    path "sample_swap_report.tsv", emit: report

    script:
    def min_rel   = params.somalier_min_relatedness ?: 0.90   // same-patient pairs should exceed this
    def max_unrel = params.somalier_max_unrelated   ?: 0.20   // unrelated pairs should stay below this
    // Same person => opposite homozygotes (IBS0) ~0 even with tumour LOH (LOH turns het->hom,
    // never hom->opposite hom). Low relatedness with IBS0 fraction <= this is LOH, not a swap.
    def max_ibs0  = (params.somalier_max_ibs0_frac != null) ? params.somalier_max_ibs0_frac : 0.02
    def analysis  = params.analysis ?: 'germline'
    def dup_min   = (params.somalier_dup_min_relatedness != null) ? params.somalier_dup_min_relatedness : 0.80
    // Somatic/tumor-only: a patient's samples must match, other patients must not. Germline:
    // relatives may exist, so only identity-level matches between sample IDs are flagged.
    """
    check_sample_swaps.sh --pairs ${pairs_tsv} --analysis ${analysis} \\
        --min-rel ${min_rel} --max-unrel ${max_unrel} --max-ibs0 ${max_ibs0} --dup-min ${dup_min}
    """

}


// ============================================================
//  STAGE 3 metrics: QC_GATE - soft QC gate per sample (reports only, never stops a run)
//
//  Evaluates per-sample coverage (mosdepth) and contamination
//  (VerifyBamID2 FREEMIX) against thresholds, emitting PASS/WARN/UNKNOWN
//  per metric. Always exits 0 - flags are advisory; you decide whether to
//  keep/remove a flagged sample. T/N concordance is assessed separately in
//  somalier.pairs.tsv (already produced). Shared QC (somatic + germline).
// ============================================================

process QC_GATE {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/gate" }, mode: 'copy'

    input:
    tuple val(meta), path(mosdepth_summary), path(selfsm),
          path(hs_metrics, stageAs: 'hs_metrics.txt'),
          path(samtools_stats, stageAs: 'samtools.stats'),
          path(global_dist,    stageAs: 'global.dist.txt')       // mosdepth global.dist (WGS breadth)

    output:
    tuple val(meta), path("${meta.sample_id}_qc_gate.tsv"), emit: report

    script:
    // Coverage floor by assay x sample type. WES: qc_min_mean_coverage for both types unless
    // qc_min_mean_coverage_normal is set. WGS: separate tumour / normal floors.
    def is_normal     = (meta.type == 'DN')
    def min_cov       = (meta.assay == 'wgs')
        ? (is_normal ? (params.qc_wgs_min_mean_coverage_normal ?: 30) : (params.qc_wgs_min_mean_coverage_tumor ?: 60))
        : ((is_normal && params.qc_min_mean_coverage_normal) ? params.qc_min_mean_coverage_normal
                                                             : (params.qc_min_mean_coverage ?: 100))
    def max_contam    = params.qc_max_contamination   ?: 0.03
    def min_breadth   = params.qc_min_pct_target_30x  ?: 0.90
    def min_breadth20 = params.qc_min_pct_target_20x  ?: 0.90
    // Breadth gate by role: tumours at 30x; normals / germline samples (WES floor 50x mean) at
    // 20x, the usual germline exome breadth metric. The other breadth row is reported as INFO.
    def gate20        = is_normal
    def min_ontarget  = params.qc_min_on_target       ?: 0.50
    def min_mapping   = params.qc_min_mapping_rate    ?: 0.95
    def max_fold80    = params.qc_max_fold80          ?: 3.0
    // WGS breadth: share of the autosomes (chromosome length, N gaps included) at >= 30x / 20x
    def is_wgs        = (meta.assay == 'wgs')
    def wgs_b30       = params.qc_wgs_min_pct_genome_30x ?: 0.80
    def wgs_b20       = params.qc_wgs_min_pct_genome_20x ?: 0.80
    """
    qc_gate.sh --sample ${meta.sample_id} --summary ${mosdepth_summary} --selfsm "${selfsm}" \\
        --hs hs_metrics.txt --stats samtools.stats \\
        --min-cov ${min_cov} --max-contam ${max_contam} \\
        --min-breadth30 ${is_wgs ? wgs_b30 : min_breadth} --min-breadth20 ${is_wgs ? wgs_b20 : min_breadth20} --gate20 ${gate20} \\
        ${is_wgs ? "--global-dist global.dist.txt" : ''} \\
        --min-ontarget ${min_ontarget} --min-mapping ${min_mapping} --max-fold80 ${max_fold80}
    """
}


// ============================================================
//  STAGE 3 metrics: QC_GATE_SUMMARY - soft QC gate - cohort summary
// ============================================================

process QC_GATE_SUMMARY {
    label 'process_low'

    // -- 03 metrics / QC (shared) --------------------------------------
    publishDir { "${pubDir('metrics')}/gate" }, mode: 'copy'

    input:
    path gate_tsvs

    output:
    path "cohort_qc_summary.tsv", emit: summary

    script:
    """
    qc_gate_summary.sh ${gate_tsvs}
    """
}


// ============================================================
//  STAGE 4 bqsr: BASE_RECALIBRATOR - BaseRecalibrator (per shard on WGS)
// ============================================================

process BASE_RECALIBRATOR {
    tag { "${meta.sample_id}:${interval_name}" }
    label 'process_medium'

    maxRetries    2

    input:
    tuple val(meta), path(bam), path(bai), val(interval_name), path(interval_file)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  dbsnp
    path  dbsnp_tbi
    path  mills
    path  mills_tbi
    path  known_indels
    path  known_indels_tbi
    path  interval_bed

    output:
    tuple val(meta), path("${meta.sample_id}.${interval_name}.recal.table"), emit: recal_table

    script:
    // WGS: the scatter shard; WES: the capture BED (single 'all' shard).
    def interval_flag = interval_file ? "-L ${interval_file}"
                      : ((meta.assay == 'wes' && interval_bed) ? "-L ${interval_bed}" : '')
    def java_opts = gatkJavaOpts(task)

    """
    set -euo pipefail

    tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_${interval_name}_bqsr_tmp.XXXXXX")
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" BaseRecalibrator \\
        -R ${ref_fasta} -I ${bam} \\
        --known-sites ${dbsnp} \\
        --known-sites ${mills} \\
        --known-sites ${known_indels} \\
        ${interval_flag} \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}.${interval_name}.recal.table
    """
}


// ============================================================
//  STAGE 4 bqsr: GATHER_BQSR_REPORTS - GatherBQSRReports - one table per sample
// ============================================================

process GATHER_BQSR_REPORTS {
    tag { "${meta.sample_id}" }
    label 'process_low'

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("${meta.sample_id}_recal.table"), emit: recal_table

    script:
    def java_opts = gatkJavaOpts(task)
    def inputs    = (tables instanceof List ? tables : [tables]).collect { t -> "-I ${t}" }.join(' ')
    """
    set -euo pipefail
    gatk --java-options "${java_opts}" GatherBQSRReports ${inputs} -O ${meta.sample_id}_recal.table
    """
}


// ============================================================
//  STAGE 4 bqsr: APPLY_BQSR - ApplyBQSR
// ============================================================

process APPLY_BQSR {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    maxRetries    2

    // Whole-BAM pass (per-shard ApplyBQSR would duplicate reads spanning shard
    // boundaries). Single-threaded; WGS needs far more wall time than WES. WGS is not
    // scaled by attempt: a retry must still fit the medium queue's 24 h limit.
    time { meta.assay == 'wgs' ? 20.h : 4.h * task.attempt }

    // -- 04 bqsr (shared) ----------------------------------------------
    publishDir { "${pubDir('bqsr')}/bqsr_bams" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  interval_bed

    output:
    // GATK ApplyBQSR produces .bai not .bam.bai
    tuple val(meta),
          path("${meta.sample_id}_bqsr.bam"),
          path("${meta.sample_id}_bqsr.bai"), emit: bam_bai_pair

    script:
    def java_opts = gatkJavaOpts(task)

    """
    set -euo pipefail

    tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_applybqsr_tmp.XXXXXX")
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" ApplyBQSR \\
        -R ${ref_fasta} -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \\
        --create-output-bam-index true \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}_bqsr.bam
    """
}


// ============================================================
//  STAGE 5 variant_calling: GETPILEUPSUMMARIES - GetPileupSummaries (tumor and normal)
//
//  Pileup at common biallelic germline SNP sites (params.pileup_variants),
//  restricted to the capture BED for WES. Feeds CalculateContamination.
//  One process, called twice (tumor BAM, normal BAM).
// ============================================================

process GETPILEUPSUMMARIES {
    tag { "${meta.sample_id}:${label}" }
    label 'process_medium'

    maxRetries    2

    input:
    tuple val(meta), val(label), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  pileup_variants
    path  pileup_variants_tbi
    path  interval_bed

    output:
    tuple val(meta), val(label), path("${meta.sample_id}.${label}.pileups.table"), emit: table

    script:
    def java_opts = gatkJavaOpts(task)
    def interval_args = (meta.assay == 'wes' && interval_bed)
        ? "-L ${interval_bed} -L ${pileup_variants} --interval-set-rule INTERSECTION"
        : "-L ${pileup_variants}"

    """
    set -euo pipefail

    tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_${label}_pile_tmp.XXXXXX")
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" GetPileupSummaries \\
        -I ${bam} \\
        -V ${pileup_variants} \\
        ${interval_args} \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}.${label}.pileups.table
    """
}


// ============================================================
//  STAGE 5 variant_calling: CALCULATE_CONTAMINATION - CalculateContamination (tumor vs matched normal)
//
//  Produces the contamination table + tumor segmentation consumed by
//  FilterMutectCalls. Matched-normal pileup sharpens the estimate.
// ============================================================

process CALCULATE_CONTAMINATION {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (contamination check) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/contamination" }, mode: 'copy', pattern: '*.{contamination,segments}.table'

    input:
    tuple val(meta), path(tumor_pileup), path(normal_pileup)

    output:
    tuple val(meta),
          path("${meta.sample_id}.contamination.table"),
          path("${meta.sample_id}.segments.table"), emit: contamination

    script:
    def mem_gb    = Math.max(task.memory.toGiga() - 1, 2)
    def java_opts = "-Xmx${mem_gb}g"

    """
    set -euo pipefail

    gatk --java-options "${java_opts}" CalculateContamination \\
        -I ${tumor_pileup} \\
        --matched-normal ${normal_pileup} \\
        --tumor-segmentation ${meta.sample_id}.segments.table \\
        -O ${meta.sample_id}.contamination.table
    """
}

// Use this later for compliant with tumor only cases


// ============================================================
//  STAGE 5 variant_calling: LEARN_READ_ORIENTATION_MODEL - LearnReadOrientationModel
//
//  Gathers all per-interval f1r2 tarballs for a sample into one
//  orientation-bias prior model. FFPE deamination produces strand-biased
//  C>T/G>A calls; --orientation-bias-artifact-priors lets FilterMutectCalls
//  flag them. This is the call-level complement to UMI consensus.
// ============================================================

process LEARN_READ_ORIENTATION_MODEL {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    maxRetries    2

    input:
    tuple val(meta), path(f1r2_tarballs, arity: '1..*')   // always a list (a single shard is still one -I)

    output:
    tuple val(meta), path("${meta.sample_id}.read-orientation-model.tar.gz"), emit: priors

    script:
    def java_opts = gatkJavaOpts(task)
    def inputs    = f1r2_tarballs.collect { v -> "-I ${v}" }.join(' \\\n        ')

    """
    set -euo pipefail

    gatk --java-options "${java_opts}" LearnReadOrientationModel \\
        ${inputs} \\
        -O ${meta.sample_id}.read-orientation-model.tar.gz
    """
}


// ============================================================
//  STAGE 5 variant_calling: MUTATIONAL_SIGNATURES - mutational signatures - COSMIC fitting
// ============================================================

process MUTATIONAL_SIGNATURES {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // -- 05 signatures (analysis-separated) ----------------------------
    publishDir { "${pubDir('variant_calling', true)}/signatures" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.sample_id}_signatures"), emit: signatures

    script:
    def genome   = params.sig_genome_build   ?: 'GRCh38'
    def cosmic_v = params.sig_cosmic_version  ?: 3.5
    // WGS must use genome-wide COSMIC signatures regardless of params.sig_exome.
    def exome    = (meta.assay == 'wgs') ? 'False'
                 : ((params.sig_exome == null ? true : params.sig_exome) ? 'True' : 'False')
    def min_snv  = params.sig_min_snv ?: 50
    """
    mutational_signatures.sh ${meta.sample_id} ${vcf} ${genome} ${cosmic_v} ${exome} ${min_snv}
    """
}


// ============================================================
//  STAGE 5 variant_calling: CNVKIT - copy number - CNVkit
//
//  FIX vs v5.0: matched normal BAM now passed and used
//  explicitly. Previously channel.value([]) was passed,
//  causing CNVkit to fall back to PON-only mode which
//  produces less accurate log2 ratios.
// ============================================================

process CNVKIT {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/cnv/cnvkit" }, mode: 'copy'

    input:
    tuple val(meta),
          path(tumor_bam), path(tumor_bai),
          path(normal_bam), path(normal_bai)
    path  ref_fasta
    path  interval_bed
    path  cnvkit_pon

    output:
    tuple val(meta),
          path("${meta.sample_id}.cns"),
          path("${meta.sample_id}.cnr"), emit: cnv

    script:
    // Argument set (PoN -> --reference only; else a reference from the normal) and the
    // tumor-only / WES-BED / WGS-PoN guards are in bin/run_cnvkit.sh.
    def normal_args = normal_bam ? "--normal ${normal_bam} --normal-bai ${normal_bai}" : ''
    def pon_arg     = cnvkit_pon ? "--pon ${cnvkit_pon}" : ''
    def bed_arg     = interval_bed ? "--intervals ${interval_bed}" : ''
    """
    run_cnvkit.sh --sample ${meta.sample_id} --tumor ${tumor_bam} --tumor-bai ${tumor_bai} \\
        ${normal_args} ${pon_arg} ${bed_arg} \\
        --assay ${meta.assay ?: 'wes'} --threads ${task.cpus} --fasta ${ref_fasta}
    """
}


// ============================================================
//  STAGE 5 variant_calling: CNVKIT_CALL - CNVkit integer copy number (purity from Sequenza when available)
//
//  -m clonal + --purity rescales log2 for tumour content before integer CN.
//  Without Sequenza: CNVkit default thresholds (as before). Sex passed when known.
//  Limitation: CNVkit centres log2 on the genome median, i.e. assumes the
//  median is diploid; for near-triploid/tetraploid tumours Sequenza CN is primary.
// ============================================================

process CNVKIT_CALL {
    tag { "${meta.sample_id}" }
    label 'process_low'

    publishDir { "${pubDir('variant_calling', true)}/cnv/cnvkit" }, mode: 'copy', pattern: '*-call.cns'

    input:
    tuple val(meta), path(cns), path(confints, stageAs: 'confints/*')

    output:
    tuple val(meta), path("${meta.sample_id}-call.cns"), emit: call_cns

    script:
    def sex_lab = (meta.sex?.toString() == '0') ? 'female' : (meta.sex?.toString() == '1') ? 'male' : 'unknown'
    // [] when Sequenza did not run; a single staged file otherwise
    def ci_list = (confints instanceof List) ? confints : [confints]
    def ci_file = ci_list ? ci_list[0].toString() : ''
    """
    cnvkit_call.sh ${meta.sample_id} ${cns} ${sex_lab} ${ci_file}
    """
}


// ============================================================
//  STAGE 5 variant_calling: MSISENSOR_PRO - microsatellite instability - MSIsensor-pro
// ============================================================

process MSISENSOR_PRO {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/msi" }, mode: 'copy'

    input:
    // normal_bam/bai carried IN the tuple; pass [] (empty) for tumor-only.
    tuple val(meta), path(tumour_bam), path(tumour_bai)
    path(normal_bam)
    path(normal_bai)
    path  msisensor_list   // scan list (paired) OR baseline dir (tumor-only)
    path  interval_bed     // WES: restrict evaluated sites to capture; [] or WGS = genome-wide

    output:
    tuple val(meta), path("${meta.sample_id}_msi*"), emit: msi

    script:
    def region_flag = (meta.assay == 'wes' && interval_bed) ? "-e ${interval_bed}" : ''
    // per-site minimum depth: MSIsensor-pro recommends 20 for WXS, 15 for WGS
    def min_cov     = (meta.assay == 'wgs') ? 15 : 20
    if (normal_bam) {
        // paired - msi submodule, msisensor_ref = scan list
        """
        set -euo pipefail
        msisensor-pro msi -d ${msisensor_list} \\
            -t ${tumour_bam} -n ${normal_bam} ${region_flag} \\
            -o ${meta.sample_id}_msi -b ${task.cpus} -c ${min_cov}
        """
    } else {
        // tumor-only - pro submodule, msisensor_ref = pre-built baseline
        """
        set -euo pipefail
        msisensor-pro pro -d ${msisensor_list} \\
            -t ${tumour_bam} ${region_flag} \\
            -o ${meta.sample_id}_msi -b ${task.cpus} -c ${min_cov}
        """
    }
}


// ============================================================
//  STAGE 5 variant_calling: MSISENSOR_GENE_MAP - unstable microsatellite sites mapped to genes
//
//  Gene mapping of MSIsensor-pro unstable microsatellite sites.
//  Reads the paired-mode _unstable file (_somatic in older releases); base R (runs in an R
//  container, NOT the MSIsensor image). Keyed on Ensembl gene ID.
// ============================================================

process MSISENSOR_GENE_MAP {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/msi" }, mode: 'copy'

    input:
    tuple val(meta), path(somatic)     // ${sample}_msi_unstable (or _somatic), extracted at the call site
    path  gene_bed                     // ensembl115_genes.bed (chr,start,end,symbol,ENSG)

    output:
    tuple val(meta), path("${meta.sample_id}_msi_genes.tsv"), emit: msi_genes

    script:
    """
    set -euo pipefail
    Rscript --vanilla "\$(command -v msi_gene_map.R)" ${somatic} ${gene_bed} ${meta.sample_id}_msi_genes.tsv
    """
}


// ============================================================
//  STAGE 5 variant_calling: CLASSIFY_MSI - MSI status
//
//  Parses the MSIsensor-pro main output (<sample>_msi): line 2, column 3
//  is the MSI score (% unstable sites). Applies a configurable cutoff
//  (default 20% per MSIsensor-pro recommendation; tumor-type dependent).
// ============================================================

process CLASSIFY_MSI {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 classification (analysis-separated) ------------------------
    publishDir { "${pubDir('variant_calling', true)}/msi" }, mode: 'copy'

    input:
    tuple val(meta), path(msi_files)

    output:
    tuple val(meta), path("${meta.sample_id}_msi_status.txt"), emit: status

    script:
    // Paired 'msi' mode = MSIsensor score (% unstable sites): >= high MSI-H, < low MSS, in
    // between MSI-INDETERMINATE (Middha et al. 2017; Niu et al. 2014). Tumour-only 'pro'
    // mode = msisensor-pro score, one cut-off (developers recommend >= 20%).
    def tumor_only = (params.analysis == 'tumor_only')
    def cut_high   = tumor_only ? params.msi_cutoff_tumor_only : params.msi_cutoff_paired_high
    def cut_low    = tumor_only ? params.msi_cutoff_tumor_only : params.msi_cutoff_paired_low
    def min_sites = (params.msi_min_sites != null) ? params.msi_min_sites : 50
    """
    classify_msi.sh --sample ${meta.sample_id} --msi ${meta.sample_id}_msi \\
        --high ${cut_high} --low ${cut_low} --min-sites ${min_sites}
    """
}


// ============================================================
//  STAGE 6 annotation: ANNOTATE_VCF - variant annotation - VEP
// ============================================================

process ANNOTATE_VCF {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 06 annotation (analysis-separated) ------------------------
    publishDir { pubDir('annotation', true) }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  vep_cache_dir
    path  ref_fasta
    path  ref_fai
    path  vep_plugins_dir
    path  clinvar
    path  clinvar_tbi
    path  cosmic
    path  cosmic_tbi
    path  gnomad_faf
    path  gnomad_faf_tbi

    output:
    tuple val(meta),
          path("${meta.sample_id}_vep.vcf.gz"),
          path("${meta.sample_id}_vep.vcf.gz.tbi"),
          path("${meta.sample_id}_vep_summary.html"), emit: annotated_vcf

    script:
    // Optional flags: empty when the param is null (toOptionalFileChannel passes []),
    // matching the pon_flag / xl_flag idiom used by MUTECT2. dbNSFP/SpliceAI plugins
    // are intentionally NOT wired here (no data files); add explicit --plugin lines
    // later if those resources are installed.
    def cache_ver = params.vep_cache_version ?: 115
    def plugins_flag = vep_plugins_dir ? "--dir_plugins ${vep_plugins_dir}" : ''
    def syn_path = "${vep_cache_dir}/homo_sapiens/${cache_ver}_GRCh38/chr_synonyms.txt"
    // ONC/ONCREVSTAT = ClinVar oncogenicity (somatic) classification; empty in releases before 2024
    def clinvar_flag = clinvar ? "--custom ${clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,ONC,ONCREVSTAT" : ''
    // COSMIC is a somatic-mutation catalog - apply only in somatic contexts, not germline.
    def cosmic_flag = (cosmic && params.analysis in ['somatic', 'tumor_only'])
        ? "--custom ${cosmic},COSMIC,vcf,exact,0,GENOME_SCREEN_SAMPLE_COUNT,TIER" : ''
    // gnomAD v4.1 grpmax filtering AF (95% CI; ClinGen-recommended for BA1/BS1): README Appendix A, step 7.6
    def faf_flag = gnomad_faf ? "--custom ${gnomad_faf},gnomAD,vcf,exact,0,fafmax_faf95_max_joint" : ''
    """
    set -euo pipefail
    syn_flag=""
    if [ -s "${syn_path}" ]; then syn_flag="--synonyms ${syn_path}"; fi
    # Ensembl NMD plugin (no data files): flags truncations predicted to escape NMD.
    nmd_flag=""
    if [ -n "${vep_plugins_dir}" ] && [ -s "${vep_plugins_dir}/NMD.pm" ]; then
        nmd_flag="--plugin NMD"
    else
        echo "WARN: NMD.pm not in vep_plugins_dir; NMD escape judged from the last exon only" >&2
    fi

    # One ALT per record: per-allele consequence, ClinVar/COSMIC exact match and the
    # classification table all need biallelic records (germline VCFs arrive unsplit).
    bcftools norm --threads ${task.cpus} -m -any -f ${ref_fasta} -Oz -o split.vcf.gz ${vcf}

    vep \\
        --input_file  split.vcf.gz \\
        --output_file ${meta.sample_id}_vep.vcf.gz \\
        --format vcf --vcf --compress_output bgzip \\
        --fork ${task.cpus} \\
        --species homo_sapiens --assembly GRCh38 \\
        --offline --cache \\
        --dir_cache     ${vep_cache_dir} \\
        --cache_version ${cache_ver} \\
        \$syn_flag \\
        --fasta ${ref_fasta} \\
        --everything \\
        --flag_pick_allele_gene \\
        --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \\
        ${plugins_flag} \\
        \$nmd_flag \\
        ${clinvar_flag} \\
        ${cosmic_flag} \\
        ${faf_flag} \\
        --stats_file  ${meta.sample_id}_vep_summary.html \\
        --buffer_size 5000 --no_escape

    tabix -p vcf ${meta.sample_id}_vep.vcf.gz
    """
}


// ============================================================
//  STAGE 6 annotation: COMPUTE_TMB - tumor mutational burden
//
//  Counts PASS nonsynonymous-coding somatic variants (worst consequence
//  per variant, via bcftools +split-vep) and divides by the non-overlapping
//  covered megabases of the capture territory. Matched-normal calling has
//  already removed germline (the gold-standard TMB method). Cutoff default
//  10 mut/Mb (FDA/KEYNOTE-158). Splice inclusion optional.
// ============================================================

process COMPUTE_TMB {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 06 annotation (analysis-separated) ------------------------
    publishDir { "${pubDir('annotation', true)}/tmb" }, mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf), path(annotated_tbi),
          path(tumor_depth,  stageAs: 'depth_tumor/*'),     // mosdepth quantized BED of the tumour, or []
          path(normal_depth, stageAs: 'depth_normal/*')     // ... of the matched normal, or []
    path  interval_bed
    path  coding_bed, stageAs: 'coding_bed/*'   // params.tmb_coding_bed (CDS BED) or [] -> capture-only territory

    output:
    tuple val(meta), path("${meta.sample_id}_tmb.txt"), emit: tmb

    script:
    // WGS reports NA (the interval BED is the CNVkit access BED, not a coding territory).
    def cutoff      = params.tmb_cutoff         ?: 10
    def incl_splice = params.tmb_include_splice ?: false
    def coding_arg  = coding_bed ? "--coding ${coding_bed}" : ''
    def capture_arg = interval_bed ? "--capture ${interval_bed}" : ''    // optional for WGS
    // callable territory: tumour depth >= filter_min_tumor_dp and normal depth >= tmb_min_normal_dp
    def min_t_dp    = (params.filter_min_tumor_dp != null) ? params.filter_min_tumor_dp : 20
    def min_n_dp    = (params.tmb_min_normal_dp   != null) ? params.tmb_min_normal_dp   : 10
    def depth_args  = (tumor_depth && normal_depth)
        ? "--tumor-depth ${tumor_depth} --normal-depth ${normal_depth} --min-tumor-depth ${min_t_dp} --min-normal-depth ${min_n_dp}" : ''
    """
    compute_tmb.sh --sample ${meta.sample_id} --vcf ${annotated_vcf} --assay ${meta.assay ?: 'wes'} \\
        ${capture_arg} ${coding_arg} ${depth_args} --cutoff ${cutoff} --include-splice ${incl_splice}
    """
}


// ============================================================
//  STAGE 7 classification: ANNOTATE_SV - SV and CNV annotation - AnnotSV
//
//  PHASE 1: SV/CNV ANNOTATION - AnnotSV
//
//  Annotates a germline SV/CNV VCF (Manta diploidSV or gCNV segments)
//  with gene overlap, ACMG/ClinGen dosage, and known-pathogenic overlap.
//  Turns raw breakpoints/segments into interpretable, reportable calls.
//  Also can be plugged in somatic framework
//  'source_tag' distinguishes the two callers feeding this same process
//  (manta | gcnv) so outputs never collide and the report can label them.
//
//  Requires: an AnnotSV container (params.container_annotsv) and the
//  AnnotSV annotation-DB directory (params.annotsv_dir), mounted read-only.
//  -SVinputInfo 1 preserves the original VCF INFO in the TSV.
// ============================================================

process ANNOTATE_SV {
    tag { "${meta.sample_id}:${source_tag}" }
    label 'process_medium'

    publishDir { "${pubDir('annotation', true)}/sv" }, mode: 'copy'

    maxRetries 2

    input:
    tuple val(meta), val(source_tag), path(vcf), path(tbi)
    path  annotsv_dir

    output:
    tuple val(meta), val(source_tag),
          path("${meta.sample_id}_${source_tag}_annotsv.tsv"), emit: annotated_sv

    script:
    // WES: restrict to 'split' (per-gene) lines is NOT forced here - full+split
    // are both emitted so whole-gene and per-exon overlaps are visible. Genome
    // build is taken from the pipeline-wide GRCh38 assumption used by VEP.
    """
    set -euo pipefail

    AnnotSV \\
        -SVinputFile ${vcf} \\
        -annotationsDir ${annotsv_dir} \\
        -genomeBuild GRCh38 \\
        -SVinputInfo 1 \\
        -outputFile ${meta.sample_id}_${source_tag}_annotsv.tsv \\
        -outputDir .

    # AnnotSV writes <outputFile> directly; guarantee the declared path exists
    # even when zero SVs pass (empty TSV with header) so the channel never stalls.
    if [ ! -s ${meta.sample_id}_${source_tag}_annotsv.tsv ]; then
        echo -e "AnnotSV_ID\\tSV_chrom\\tSV_start\\tSV_end\\tSV_type\\tGene_name\\tACMG_class" \\
            > ${meta.sample_id}_${source_tag}_annotsv.tsv
    fi
    """
}


// ============================================================
//  STAGE 7 classification: CLASSIFY_VARIANTS - variant pre-classification (not a clinical tier)
//
//  Aggregates automatable evidence - gene role/tier (COSMIC Cancer Gene
//  Census), VEP consequence/impact, population frequency, ClinVar (if
//  annotated) - into a computational pre-classification per variant.
//  This is decision SUPPORT for a molecular pathologist, NOT a final
//  AMP/ASCO/CAP tier (which requires a curated drug-variant-disease
//  knowledge base and expert review).
// ============================================================

process CLASSIFY_VARIANTS {
    tag { "${meta.sample_id}" }
    label 'process_low'

    publishDir { pubDir('classification', true) }, mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf), path(annotated_tbi)
    path  cgc_tsv
    path  af_thresholds, stageAs: 'af_thr/*'    // params.germline_af_thresholds (GENE BA1 BS1) or []

    output:
    tuple val(meta), path("${meta.sample_id}_preclassification.tsv"), emit: classification

    script:
    def af_field          = params.classify_af_field    ?: 'MAX_AF'
    def af_cutoff         = params.classify_common_af    ?: 0.01
    def cosmic_min        = params.cosmic_recurrence_min ?: 10
    def cosmic_strong     = params.cosmic_strong_min     ?: 50
    def clinvar_min_stars = (params.clinvar_min_stars != null) ? params.clinvar_min_stars : 1   // 0 = gating off (elvis would turn 0 into 1)
    // germline rules for normals (germline arm, or matched normals with germline_on_normals)
    def cls_mode          = (params.analysis == 'germline' || meta.type == 'DN') ? 'germline' : 'somatic'
    def ba1_af            = (params.germline_ba1_af != null) ? params.germline_ba1_af : 0.05
    def germ_low_vaf      = (params.germline_low_vaf != null) ? params.germline_low_vaf : 0.30
    def mosaic_min_dp     = (params.germline_mosaic_min_dp != null) ? params.germline_mosaic_min_dp : 20
    """
    classify_variants.sh --sample ${meta.sample_id} --vcf ${annotated_vcf} --mode ${cls_mode} \\
        ${cgc_tsv ? "--cgc ${cgc_tsv}" : ''} ${af_thresholds ? "--af-thresholds ${af_thresholds}" : ''} \\
        --af-field ${af_field} --af-cutoff ${af_cutoff} \\
        --cosmic-min ${cosmic_min} --cosmic-strong ${cosmic_strong} --min-stars ${clinvar_min_stars} \\
        --ba1 ${ba1_af} --germ-low-vaf ${germ_low_vaf} --mosaic-min-dp ${mosaic_min_dp}
    """
}


// ============================================================
//  STAGE 7 classification: VCF2MAF - VCF to MAF
// ============================================================

process VCF2MAF {
    tag { "${meta.sample_id}" }
    label 'process_tiny'

    // -- 07 classification (analysis-separated)----------------
    publishDir { "${pubDir('classification', true)}/maf" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta), path("${meta.sample_id}.maf"), emit: maf

   script:
    def is_germline = (params.analysis == 'germline' || meta.type == 'DN')
    def tumor_id  = "${meta.sample_id}"
    def normal_id = is_germline ? "${meta.sample_id}" : meta.normal_id

    """
    set -euo pipefail

    zcat ${vcf} > input.vcf

    # Align the MAF transcript with the VEP PICK transcript that classification and TMB use
    # (vcf2maf otherwise applies its own canonical-transcript choice). No PICK -> unchanged.
    awk -F'\\t' '
      /^##INFO=<ID=CSQ/ { s = \$0; sub(/.*Format: /, "", s); sub(/">.*/, "", s); n = split(s, f, "|")
                          for (i = 1; i <= n; i++) { if (f[i] == "PICK") ip = i; if (f[i] == "Feature") ift = i }; next }
      /^#/ { next }
      ip && ift { m = split(\$8, kv, ";"); for (k = 1; k <= m; k++) if (kv[k] ~ /^CSQ=/) {
                    nc = split(substr(kv[k], 5), cs, ","); for (j = 1; j <= nc; j++) { split(cs[j], x, "|")
                        if (x[ip] == "1" && x[ift] != "") print x[ift] } } }' input.vcf | sort -u > pick_enst.txt
    enst_flag=""
    if [ -s pick_enst.txt ]; then enst_flag="--custom-enst pick_enst.txt"; fi

    vcf2maf.pl \\
        --input-vcf      input.vcf \\
        --output-maf     full.maf \\
        --ref-fasta      ${ref_fasta} \\
        --tumor-id       ${tumor_id} \\
        --normal-id      ${normal_id} \\
        --vcf-tumor-id   ${tumor_id} \\
        --vcf-normal-id  ${normal_id} \\
        --ncbi-build     GRCh38 \\
        --inhibit-vep \$enst_flag \\
        --retain-fmt     GT,AD,AF,DP \\
        --retain-ann     ClinVar_CLNSIG,ClinVar_CLNREVSTAT,MAX_AF,gnomAD_fafmax_faf95_max_joint


    PY=\$(command -v python3 || echo /opt/conda/bin/python3)
    "\$PY" "\$(command -v trim_maf.py)" ${meta.sample_id}.maf ${is_germline ? 'germline' : 'somatic'}
    """
}


// ============================================================
//  STAGE 0 setup: PREFLIGHT_REFERENCES - reference checks before any job is submitted
//
//  PREFLIGHT_REFERENCES - fail fast on genome-build / contig-naming mismatches
//  between the reference .dict and the resource VCFs / capture BED. Runs once,
//  in parallel with mapping; 'terminate' stops the run immediately on failure
//  instead of letting hours of alignment finish first.
// ============================================================

process PREFLIGHT_REFERENCES {
    tag 'references'
    label 'process_tiny'
    errorStrategy 'terminate'

    publishDir { "${pubDir('manifest', true)}" }, mode: 'copy'

    input:
    path ref_dict
    path vcfs, stageAs: 'vcf/*'        // resource VCFs that are set (dbsnp, gnomAD, PON, ...)
    path beds, stageAs: 'bed/*'        // interval_bed (+ call_regions_gz)

    output:
    path 'preflight_references.txt', emit: report

    script:
    """
    preflight_references.sh ${ref_dict} vcf bed
    """
}


// ============================================================
//  STAGE 5 variant_calling: EVIDENCE_BAM - read slices around given regions
//
//  EVIDENCE_BAM: reads overlapping the evidence regions (BED) for a specific project, indexed.
//  BAM or CRAM input (CRAM decoded with the reference).
// ============================================================

process EVIDENCE_BAM {
    tag { "${meta.sample_id}" }
    label 'process_low'

    publishDir { "${pubDir('variant_calling', true)}/evidence_bam" }, mode: 'copy'

    input:
    tuple val(meta), path(aln), path(idx)
    path  regions_bed
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta), path("${meta.sample_id}.evidence.bam"), path("${meta.sample_id}.evidence.bam.bai"), emit: bam

    script:
    """
    set -euo pipefail

    # -M: indexed multi-region iterator (without it samtools reads the whole BAM)
    samtools view -@ ${task.cpus} -T ${ref_fasta} -M -L ${regions_bed} -b -o ${meta.sample_id}.evidence.bam ${aln}
    samtools index -@ ${task.cpus} ${meta.sample_id}.evidence.bam
    """
}
