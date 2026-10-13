#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Shared JVM-options helper.
include { gatkJavaOpts; pubDir; umiLeadingSkip} from './utils'
// ============================================================
//  components.nf  -  All Process Definitions
// ============================================================

// ============================================================
//  STAGE 0A: Run manifest (per-sample audit record)
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
//  STAGE 0B: Write manifest (per-sample audit record)
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
//   Multisatge: MultiQC Aggregation
// ============================================================

process MULTIQC {
    tag { "${label}" }
    label 'process_low'

    // MultiQC aggregates ACROSS all samples of a run - it is not a per-stage
    // -- 0x multiqc  -----------------------------------------
    publishDir { "${params.outdir}/0x_multiqc/${params.assay == 'rna' ? 'rna_fusion' : params.analysis}" }, mode: 'copy'

    input:
    tuple val(label), path(files)

    output:
    path "${label}_multiqc_report.html", emit: report

    script:
    """
    set -euo pipefail
    multiqc ${files.join(' ')} --filename ${label}_multiqc_report.html --force
    """
}


// ============================================================
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
//  STAGE 1A ALIGN DNA: TRIM + ALIGN + FIXMATE  (streaming single pass)

//  Architecture: fastp - bwa-mem2 - samtools fixmate
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
//  STAGE 1B ALIGN UMI: EXTRACT UMI + ALIGN  (dedup_mode = umi_*)
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
    def fastq_heap   = 8  
    def zipper_heap  = 4  
    
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
        --read-structures ${params.umi_read_structure ?: '3M2S+T'} ${params.umi_read_structure_r2 ?: (params.umi_read_structure ?: '3M2S+T')} \\
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
    | samtools sort -@ 2 -m 4G -T "\${tmp_dir}"/sam_sort_tmp -o ${meta.uid}.umi.mapped.bam -

    # Clean up intermediate unmapped data
    rm -f ./chunk_unmapped.bam
    rm -rf "\${tmp_dir}"
    """
}


// ============================================================
//  STAGE 2A MARKDUP: MERGE + SORT + MARK DUPLICATES
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
//  STAGE 2B CONSENSUS UMI: MERGE + UMI CONSENSUS  (dedup_mode = umi_consensus)
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

    // Initial standalone sort buffer (14 threads x 5G = 70GB footprint)
    def phase1_sort_mem  = "5G"

    // Streamed Java heaps balanced to fit inside 150GB
    def group_heap       = 24
    def call_heap        = 16
    def filter_heap      = 16
    def zipper_heap      = 8
    def phase4_sort_mem  = "1500M" // ~21GB final coordinate sort pool

    // Consensus parameters
    def min_reads    = params.consensus_min_reads ?: 1
    def err_pre      = params.consensus_error_rate_pre  ?: 45
    def err_post     = params.consensus_error_rate_post ?: 40
    def max_read_err = params.consensus_max_read_error ?: 0.05
    def min_base_q   = params.consensus_min_base_quality ?: 20

    """
    set -euo pipefail

    # Setup temporary directory structure based on your format
    tmp_dir="./mer_con_tmp"
    mkdir -p "\$tmp_dir"
    
    # Explicitly force the JVM runtime environment to use your workspace folder
    export TMPDIR="\$tmp_dir"
    
    # Safety trap to clear out all intermediate assets on exit
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

    # Wipe the raw template file immediately to free up workspace disk space
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

    # Clear out the intermediate grouped file
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
//  STAGE 2C MARKDUP UMI: MERGE + UMI-AWARE MARKDUP  (dedup_mode = umi_markdup)
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
    def system_buffer = Math.max((task.memory.toGiga() * 0.15).toInteger(), 10)
    def mem_gb        = (task.memory.toGiga() - system_buffer).toInteger()
    // def java_opts     = "-Xmx${mem_gb}g -XX:ParallelGCThreads=${parallel_gc_threads} -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_read_samtools=true"
    def java_opts = "-Xmx${mem_gb}g -XX:ParallelGCThreads=${parallel_gc_threads}"
    def bam_list      = chunk_bams instanceof List ? chunk_bams : [chunk_bams]
    def input_list    = bam_list.collect { "-I $it" }.join(' ')
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
//  STAGE 3A MOSEDEPTH QC: COVERAGE QC - mosdepth
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

    script:
    // WES: per-target regions (--by capture BED) enable --thresholds + --quantize.
    // WGS: no capture regions. mosdepth --thresholds REQUIRES --by, and --fast-mode
    //      is incompatible with --quantize - so WGS uses windowed --by instead, which
    //      keeps thresholds/quantize working and gives genome coverage in windows.
    def is_wgs = (meta.assay == 'wgs')
    def by_flag = is_wgs
        ? "--by ${params.mosdepth_wgs_window ?: 1000}"
        : ((meta.assay == 'wes' && interval_bed) ? "--by ${interval_bed}" : '')
    // --thresholds is a TARGET-capture metric (% of target bases >= X). Meaningful
    // for WES only. On WGS windows it measures "% of arbitrary windows crossing X",
    // which is not comparable to the WES gate - so WGS coverage QC uses summary +
    // global.dist instead. --thresholds also REQUIRES --by, hence the WES+BED guard.
    def thresholds_flag = (!is_wgs && by_flag) ? "--thresholds 1,10,20,30,50,100" : ''
    // Per-gene coverage of reportable genes (WES + WGS): a negative result is only
    // reportable where the gene is covered. MAPQ-filtered (usable) depth.
    def gene_mapq  = (params.coverage_report_min_mapq != null) ? params.coverage_report_min_mapq : 20
    def gene_block = !coverage_report_bed ? '' : """
    gene_md() {
        mosdepth \\
            --threads ${task.cpus} --no-per-base --mapq ${gene_mapq} \\
            --by "\$1" --thresholds 10,20,30,50,100 \\
            ${meta.sample_id}.genes ${bam}
    }
    gene_md ${coverage_report_bed}
    if [ -z "\$(zcat ${meta.sample_id}.genes.regions.bed.gz | head -1)" ]; then
        # no region matched the BAM: most often chr1-vs-1 naming. Retry once with toggled names.
        echo "WARN: coverage_report_bed matched no BAM contig - retrying with toggled chr naming" >&2
        awk -F'\\t' 'BEGIN { OFS = "\\t" } /^(#|track|browser)/ { next }
            { c = \$1; d = (c ~ /^chr/) ? substr(c, 4) : "chr" c
              if (d == "chrMT") d = "chrM"; else if (d == "M") d = "MT"; \$1 = d; print }' ${coverage_report_bed} > report.toggled.bed
        gene_md report.toggled.bed
    fi
    if [ -z "\$(zcat ${meta.sample_id}.genes.regions.bed.gz | head -1)" ]; then
        echo "WARN: coverage_report_bed shares no contig with the BAM - gene coverage NOT computed" >&2
        printf '# ERROR: coverage_report_bed shares no contig with the BAM (check chr naming / build)\\n' > ${meta.sample_id}.gene_coverage.tsv
    else

    # length-weighted per-gene summary (col4 of the report BED = gene)
    zcat ${meta.sample_id}.genes.regions.bed.gz \\
        | awk -F'\\t' '{ L = \$3 - \$2; bp[\$4] += L; dsum[\$4] += L * \$5 }
                   END { for (g in bp) printf "%s\\t%d\\t%.2f\\n", g, bp[g], dsum[g] / bp[g] }' \\
        | LC_ALL=C sort -k1,1 > depth.tmp
    zcat ${meta.sample_id}.genes.thresholds.bed.gz \\
        | awk -F'\\t' '!/^#/ { g = \$4; for (i = 5; i <= 9; i++) t[g, i] += \$i }
                   END { for (k in t) { split(k, a, SUBSEP); seen[a[1]] = 1 }
                         for (g in seen) printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", g, t[g,5], t[g,6], t[g,7], t[g,8], t[g,9] }' \\
        | LC_ALL=C sort -k1,1 > thr.tmp
    {   printf 'gene\\ttarget_bp\\tmean_depth\\tpct_10x\\tpct_20x\\tpct_30x\\tpct_50x\\tpct_100x\\n'
        LC_ALL=C join -t \$'\\t' depth.tmp thr.tmp \\
            | awk -F'\\t' 'BEGIN { OFS = "\\t" }
                { L = \$2; printf "%s\\t%d\\t%s", \$1, L, \$3
                  for (i = 4; i <= 8; i++) printf "\\t%.4f", (L > 0 ? \$i / L : 0); printf "\\n" }'
    } > ${meta.sample_id}.gene_coverage.tsv
    rm -f depth.tmp thr.tmp
    fi
    """
    """
    set -euo pipefail
    mosdepth \\
        --threads ${task.cpus} \\
        --quantize 0:1:10:20:30:50:100: \\
        ${thresholds_flag} \\
        ${by_flag} \\
        ${meta.sample_id} ${bam}
    ${gene_block}
    """

}



// ============================================================
//  STAGE 3B: INFER SEX - mosdepth
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
    """
    set -euo pipefail

    # WES (--by capture): use the *_region (on-target) rows - whole-chromosome means are
    # dominated by off-target reads. chr prefix optional (GRCh38 chr / b37 style).
    if grep -q '_region' ${mosdepth_summary}; then suf="_region"; else suf=""; fi
    auto=\$(awk -v s="\$suf" '{ c = \$1; if (s != "") { if (c !~ /_region\$/) next; sub(/_region\$/, "", c) } else if (c ~ /_region\$/) next }
                  c ~ /^(chr)?([1-9]|1[0-9]|2[0-2])\$/ { sum += \$4; n++ }
                  END { if (n > 0) printf "%.4f", sum / n; else print "NA" }' ${mosdepth_summary})
    xdep=\$(awk -v s="\$suf" '\$1 == ("chrX" s) || \$1 == ("X" s) { print \$4; exit }' ${mosdepth_summary})
    xdep="\${xdep:-NA}"

    ratio="NA"; inferred="UNKNOWN"
    if [ "\$auto" != "NA" ] && [ "\$xdep" != "NA" ] && awk "BEGIN{exit !(\$auto>0)}"; then
        ratio=\$(awk -v x="\$xdep" -v a="\$auto" 'BEGIN{printf "%.3f", x/a}')
        inferred=\$(awk -v r="\$ratio" -v fmin="${female_min}" -v mmax="${male_max}" \
            'BEGIN{ if(r+0>=fmin+0) print "female"; else if(r+0<=mmax+0) print "male"; else print "AMBIGUOUS" }')
    fi

    case "${meta.sex}" in
        0) expected="female" ;;
        1) expected="male" ;;
        *) expected="unknown" ;;
    esac

    if [ "\$inferred" = "UNKNOWN" ] || [ "\$inferred" = "AMBIGUOUS" ] || [ "\$expected" = "unknown" ]; then
        flag="SKIP"
    elif [ "\$inferred" = "\$expected" ]; then
        flag="PASS"
    elif [ "${meta.type}" = "DT" ]; then
        flag="REVIEW_TUMOR"    # chrX loss / loss of Y in tumours can explain it - judge on the normal
    else
        flag="MISMATCH"
    fi

    {   echo -e "sample_id\\t${meta.sample_id}"
        echo -e "chrX_autosome_ratio\\t\$ratio"
        echo -e "inferred_sex\\t\$inferred"
        echo -e "expected_sex\\t\$expected"
        echo -e "sex_concordance\\t\$flag"
    } > ${meta.sample_id}_sex.tsv
    """
}

// ============================================================
//  STAGE 3C SAMTOOLS STATS QC: ALIGNMENT STATS
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
//  STAGE 3D METRICS QC: METRICS - CollectMultipleMetrics + CollectHsMetrics
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

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_metrics_tmp"
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
//  STAGE 3E VERIFYBAMID2 QC: CONTAMINATION QC
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
//  STAGE 3F: QC SAMPLE-SWAP DETECTION - somalier  
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
//  STAGE 3G: QC SAMPLE-SWAP DETECTION - somalier  
//============================================================

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
//  STAGE 3H : QC SAMPLE-SWAP DETECTION - somalier  
//============================================================
//         END { print "# flagged_pairs="flagged"\\t(threshold="min")" > "/dev/stderr" }
//     ' ${pairs_tsv} > sample_swap_report.tsv

//     flagged=\$(grep -c "FLAG_LOW_RELATEDNESS" sample_swap_report.tsv || true)
//     echo "CHECK_SAMPLE_SWAPS: \${flagged} pair(s) below relatedness ${min_rel}" >&2
//     """
// }


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
    // Somatic/tumor-only: each patient's samples MUST match, cross-patient must not.
    // Germline: samples are usually 1-per-patient and relatives may legitimately
    // exist, so we only surface UNEXPECTED identity-level matches (possible
    // duplicates/swaps) and never flag mere low relatedness between distinct people.
    def is_germline = (analysis == 'germline')
    """
    set -euo pipefail

    hdr=\$(head -1 ${pairs_tsv} | sed 's/^#//')
    col() { echo "\$hdr" | tr '\\t' '\\n' | grep -nx "\$1" | cut -d: -f1; }
    ca=\$(col sample_a); cb=\$(col sample_b); cr=\$(col relatedness)
    ci0=\$(col ibs0); cn=\$(col n)     # optional: absent -> relatedness-only rule
    if [ -z "\$ca" ] || [ -z "\$cb" ] || [ -z "\$cr" ]; then
        echo "ERROR: could not find sample_a/sample_b/relatedness columns in ${pairs_tsv}" >&2
        head -1 ${pairs_tsv} >&2; exit 1
    fi

    awk -F'\\t' -v ca="\$ca" -v cb="\$cb" -v cr="\$cr" -v ci0="\${ci0:-0}" -v cn="\${cn:-0}" \\
        -v minrel=${min_rel} -v maxunrel=${max_unrel} -v maxibs0=${max_ibs0} -v germ="${is_germline}" '
        function patient(s,   n){ n=s; sub(/_[^_]+\$/, "", n); return n }
        BEGIN { OFS="\\t"; print "sample_a","sample_b","relatedness","same_patient","verdict"
                swap=0; contam=0; dup=0; review=0 }
        NR>1 && \$ca!="" {
            a=\$ca; b=\$cb; rel=\$cr+0
            same=(patient(a)==patient(b)) ? "yes" : "no"
            verdict="PASS"

            if (germ=="true") {
                # Germline: no T/N expectation. Only flag identity-level matches
                # between what should be different samples (duplicate/swap). A
                # relatedness ~1.0 between two distinct sample IDs is the signal;
                # family relationships (~0.5/0.25) are expected and NOT flagged.
                if (rel > 0.80) { verdict="FLAG_POSSIBLE_DUPLICATE"; dup++;
                                  print a,b,\$cr,same,verdict }
            } else {
                # Somatic / tumor-only: same-patient pairs must be related;
                # different-patient pairs must not be.
                if (same=="yes" && rel < minrel) {
                    if (ci0 > 0 && cn > 0 && \$cn + 0 > 0 && (\$ci0 / \$cn) <= maxibs0) {
                        verdict="REVIEW_LOW_RELATEDNESS_IBS0_OK"; review++      # LOH-compatible, same person
                    } else { verdict="FLAG_SWAP_LOW_RELATEDNESS"; swap++ }
                }
                if (same=="no"  && rel > maxunrel) { verdict="FLAG_UNEXPECTED_RELATEDNESS"; contam++ }
                if (same=="yes" || verdict!="PASS") print a,b,\$cr,same,verdict
            }
        }
        END {
            if (germ=="true")
                print "# possible_duplicates="dup"  (relatedness>0.80 between distinct sample IDs)" > "/dev/stderr"
            else
                print "# REVIEW_LOW_RELATEDNESS_IBS0_OK = IBS0~0: tumour LOH in one person, OR a parent-child pair (also IBS0~0) - confirm identity"
                print "# same_patient_swaps="swap"  loh_compatible_reviews="review"  unexpected_relatedness_pairs="contam \\
                      "  (min_rel="minrel", max_unrel="maxunrel")" > "/dev/stderr"
        }
    ' ${pairs_tsv} > sample_swap_report.tsv

    n=\$(( \$(wc -l < sample_swap_report.tsv) - 1 ))
    echo "CHECK_SAMPLE_SWAPS (${analysis}): \${n} flagged/relevant pair(s)" >&2
    """

}


// ============================================================
//  STAGE 3I: Soft QC gate (report-only - NEVER halts)
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
          path(samtools_stats, stageAs: 'samtools.stats')

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
    def min_ontarget  = params.qc_min_on_target       ?: 0.50
    def min_mapping   = params.qc_min_mapping_rate    ?: 0.95
    def max_fold80    = params.qc_max_fold80          ?: 3.0
    """
    set -euo pipefail

    cov=\$(awk -F'\\t' '\$1=="total_region"{print \$4}' ${mosdepth_summary})
    if [ -z "\${cov:-}" ]; then cov=\$(awk -F'\\t' '\$1=="total"{print \$4}' ${mosdepth_summary}); fi
    cov="\${cov:-NA}"

    contam="NA"
    if [ -s "${selfsm}" ]; then
        contam=\$(awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="FREEMIX") c=i} NR==2 && c{print \$c}' ${selfsm})
        contam="\${contam:-NA}"
    fi

    # Picard HsMetrics (WES only): breadth, on-target, fold-80 uniformity.
    breadth="NA"; ontarget="NA"; fold80="NA"
    if [ -s "hs_metrics.txt" ]; then
        breadth=\$(awk -F'\\t'  '/^BAIT_SET/{for(i=1;i<=NF;i++) if(\$i=="PCT_TARGET_BASES_30X")c=i; getline; if(c)print \$c}' hs_metrics.txt)
        ontarget=\$(awk -F'\\t' '/^BAIT_SET/{for(i=1;i<=NF;i++) if(\$i=="PCT_SELECTED_BASES")c=i; getline; if(c)print \$c}' hs_metrics.txt)
        fold80=\$(awk -F'\\t'   '/^BAIT_SET/{for(i=1;i<=NF;i++) if(\$i=="FOLD_80_BASE_PENALTY")c=i; getline; if(c)print \$c}' hs_metrics.txt)
        breadth="\${breadth:-NA}"; ontarget="\${ontarget:-NA}"; fold80="\${fold80:-NA}"
    fi

    # samtools stats: mapping rate = reads mapped / raw total (SN section).
    mapping="NA"
    if [ -s "samtools.stats" ]; then
        mapped=\$(awk -F'\\t' '/^SN\\treads mapped:/{print \$3; exit}' samtools.stats)
        total=\$(awk -F'\\t'  '/^SN\\traw total sequences:/{print \$3; exit}' samtools.stats)
        if [ -n "\${mapped:-}" ] && [ -n "\${total:-}" ] && [ "\${total:-0}" != "0" ]; then
            mapping=\$(awk -v m="\$mapped" -v t="\$total" 'BEGIN{printf "%.4f", m/t}')
        fi
    fi

    cov_flag=\$(awk -v v="\$cov" -v t="${min_cov}" 'BEGIN{print (v=="NA")?"UNKNOWN":((v+0>=t+0)?"PASS":"WARN")}')
    con_flag=\$(awk -v v="\$contam" -v t="${max_contam}" 'BEGIN{print (v=="NA")?"UNKNOWN":((v+0<=t+0)?"PASS":"WARN")}')
    brd_flag=\$(awk -v v="\$breadth" -v t="${min_breadth}" 'BEGIN{print (v=="NA")?"SKIP":((v+0>=t+0)?"PASS":"WARN")}')
    ont_flag=\$(awk -v v="\$ontarget" -v t="${min_ontarget}" 'BEGIN{print (v=="NA")?"SKIP":((v+0>=t+0)?"PASS":"WARN")}')
    map_flag=\$(awk -v v="\$mapping" -v t="${min_mapping}" 'BEGIN{print (v=="NA")?"UNKNOWN":((v+0>=t+0)?"PASS":"WARN")}')
    f80_flag=\$(awk -v v="\$fold80" -v t="${max_fold80}" 'BEGIN{print (v=="NA")?"SKIP":((v+0<=t+0)?"PASS":"WARN")}')

    {   echo -e "sample_id\\t${meta.sample_id}"
        echo -e "metric\\tvalue\\tthreshold\\tflag"
        echo -e "mean_coverage\\t\$cov\\t>=${min_cov}\\t\$cov_flag"
        echo -e "contamination_freemix\\t\$contam\\t<=${max_contam}\\t\$con_flag"
        echo -e "pct_target_30x\\t\$breadth\\t>=${min_breadth}\\t\$brd_flag"
        echo -e "pct_on_target\\t\$ontarget\\t>=${min_ontarget}\\t\$ont_flag"
        echo -e "mapping_rate\\t\$mapping\\t>=${min_mapping}\\t\$map_flag"
        echo -e "fold_80_penalty\\t\$fold80\\t<=${max_fold80}\\t\$f80_flag"
        echo -e "note\\tT/N_concordance_assessed_in_somalier.pairs.tsv\\t-\\t-"
    } > ${meta.sample_id}_qc_gate.tsv

    exit 0
    """
}


// ============================================================
//  STAGE 3J QC SUMMARY: Soft QC gate (report-only - NEVER halts)
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
    set -euo pipefail

    metrics="mean_coverage contamination_freemix pct_target_30x pct_on_target mapping_rate fold_80_penalty"

    { printf "SAMPLE"; for m in \$metrics; do printf "\\t%s" "\$m"; done; printf "\\tOVERALL\\n"; } > cohort_qc_summary.tsv

    for f in ${gate_tsvs}; do
        [ -s "\$f" ] || continue
        sid=\$(awk -F'\\t' '\$1=="sample_id"{print \$2; exit}' "\$f")
        overall="PASS"
        line="\$sid"
        for m in \$metrics; do
            val=\$(awk -F'\\t' -v M="\$m" '\$1==M{print \$2; exit}' "\$f")
            flag=\$(awk -F'\\t' -v M="\$m" '\$1==M{print \$4; exit}' "\$f")
            val="\${val:-NA}"; flag="\${flag:-NA}"
            line="\$line\t\$val"
            if [ "\$flag" = "WARN" ]; then overall="WARN"; fi
        done
        printf "%b\\t%s\\n" "\$line" "\$overall" >> cohort_qc_summary.tsv
    done

    head -1 cohort_qc_summary.tsv > .hdr
    tail -n +2 cohort_qc_summary.tsv | sort -t\$'\\t' -k1,1 > .body
    cat .hdr .body > cohort_qc_summary.tsv
    rm -f .hdr .body
    """
}

// ============================================================
//  STAGE 4A: BQSR - BaseRecalibrator
// ============================================================
process BASE_RECALIBRATOR {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    maxRetries    2

    input:
    tuple val(meta), path(bam), path(bai)
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
    tuple val(meta), path("${meta.sample_id}_recal.table"), emit: recal_table

    script:
    def interval_flag = (meta.assay == 'wes' && interval_bed) ? "-L ${interval_bed}" : ''
    def java_opts = gatkJavaOpts(task)

    """
    set -euo pipefail

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_bqsr_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" BaseRecalibrator \\
        -R ${ref_fasta} -I ${bam} \\
        --known-sites ${dbsnp} \\
        --known-sites ${mills} \\
        --known-sites ${known_indels} \\
        ${interval_flag} \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}_recal.table
    """
}


// ============================================================
//  STAGE 4B: BQSR - ApplyBQSR
// ============================================================
process APPLY_BQSR {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    maxRetries    2

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

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_applybqsr_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" ApplyBQSR \\
        -R ${ref_fasta} -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        --create-output-bam-index true \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}_bqsr.bam
    """
}




// // ============================================================
// //  STAGE 5D: SOMATIC ENSEMBLE - Mutect2 + Strelka2
// //  (bcftools norm fix applied - see CHANGELOG v5.1)
// // ============================================================



// ============================================================
//  STAGE 6A-LOH: FACETS - allele-specific CN + LOH + purity/ploidy
//  (Phase 4, tumor-normal)
//
//  FACETS infers allele-specific copy number from the B-allele
//  frequencies at heterozygous germline SNPs, comparing tumor vs
//  matched normal. This is what CNVkit cannot give you:
//    - LOH (copy-neutral and copy-loss)
//    - allele-specific (major/minor) copy number per segment
//    - tumor purity + ploidy
//  Matches the FACETS/Sequenza approach used by BostonGene.
//
//  Two steps in one process:
//    1. snp-pileup  - pile up reads at polymorphic SNP sites (dbSNP)
//    2. FACETS R    - two-pass (purity run seeds high-sensitivity run)
//
//  facets_snps = a polymorphic-SNP VCF (use the SAME dbSNP you
//  reconciled to d1.vd1). cval = segmentation stringency.
// ============================================================


// ============================================================
//  STAGE 5L-CONTAM-1: GetPileupSummaries (tumor and normal)
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

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_${label}_pile_tmp"
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
//  STAGE 5M-CONTAM-2: CalculateContamination (tumor vs matched normal)
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
//  STAGE 5N-ORIENT: LearnReadOrientationModel (per sample)
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
    tuple val(meta), path(f1r2_tarballs)

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
//  STAGE 5A: Mutational signatures (COSMIC fitting)
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
    def cosmic_v = params.sig_cosmic_version  ?: 3.4
    // WGS must use genome-wide COSMIC signatures regardless of params.sig_exome.
    def exome    = (meta.assay == 'wgs') ? 'False'
                 : ((params.sig_exome == null ? true : params.sig_exome) ? 'True' : 'False')
    def min_snv  = params.sig_min_snv ?: 50
    """
    set -euo pipefail
    mkdir -p sig_input ${meta.sample_id}_signatures

    # PASS-only (upstream filters run in 'tag' mode). bcftools reads .vcf and .vcf.gz.
    bcftools view -f PASS,. "${vcf}" -O v -o sig_input/${meta.sample_id}.vcf

    # Guard: COSMIC signature assignment is not statistically meaningful below a
    # minimum SNV count. Skip cleanly (an empty/low-count VCF returns noise or errors).
    n_snv=\$(bcftools view -H -v snps sig_input/${meta.sample_id}.vcf | wc -l)
    echo "signatures: \$n_snv SNVs (floor ${min_snv})"
    if [ "\$n_snv" -lt ${min_snv} ]; then
        printf 'SKIPPED: %s SNVs (< %s) - signature fit not meaningful.\\nn_snv\\t%s\\n' \\
            "\$n_snv" "${min_snv}" "\$n_snv" > ${meta.sample_id}_signatures/INSUFFICIENT_MUTATIONS.txt
        exit 0
    fi

    python3 - <<'PYEOF'
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(
    samples="sig_input",
    output="${meta.sample_id}_signatures",
    input_type="vcf",
    context_type="96",
    genome_build="${genome}",
    cosmic_version=${cosmic_v},
    exome=${exome},
    make_plots=True,
    verbose=False,
)
PYEOF
    """
}

// ============================================================
//  STAGE 5B CNV: CNVkit (tumor-normal)
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
          path("${meta.sample_id}.cnr"),
          path("${meta.sample_id}-call.cns"), emit: cnv

    script:
    def paired  = normal_bam ? 'true' : 'false'
    def has_pon = cnvkit_pon ? 'true' : 'false'
    def is_wgs  = (meta.assay == 'wgs') ? 'true' : 'false'

    // cnvkit rejects --reference together with --targets/--fasta, so PoN presence
    // determines the ENTIRE argument set (not two independent flags):
    //   - PoN given      - --reference ONLY (it already encodes bins + genome)
    //   - no PoN         - build a reference on the fly from the normal:
    //                        --normal (if paired) + binning (WES --targets / WGS
    //                        --method wgs [+ --access]) + --fasta
    def batch_args
    if (cnvkit_pon) {
        batch_args = "--reference ${cnvkit_pon}"
    } else {
        def normal_arg  = normal_bam ? "--normal normal.bam " : ""
        def binning_arg = (meta.assay == 'wgs')
            ? ("--method wgs" + (interval_bed ? " --access ${interval_bed}" : ""))
            : "--targets ${interval_bed}"
        batch_args = "${normal_arg}${binning_arg} --fasta ${ref_fasta}"
    }

    // Only symlink the normal when one is present.
    def link_normal = normal_bam
        ? "ln -sf ${normal_bam} normal.bam && ln -sf ${normal_bai} normal.bam.bai"
        : "true"
    """
    set -euo pipefail

    # ---- guards: fail loud rather than emit silently-wrong copy-number ----
    if [ "${paired}" != "true" ] && [ "${has_pon}" != "true" ]; then
        echo "ERROR: CNVkit tumor-only requires a panel-of-normals (cnvkit_pon) as --reference." >&2
        exit 1
    fi
    if [ "${has_pon}" != "true" ] && [ "${is_wgs}" != "true" ] && [ ! -s "${interval_bed}" ]; then
        echo "ERROR: assay=wes without a PoN requires a target/interval BED (--targets)." >&2
        exit 1
    fi
    if [ "${is_wgs}" = "true" ] && [ "${has_pon}" = "true" ]; then
        # An exome-style reference has 'Antitarget' bins in the gene column of its
        # DATA rows; a WGS reference does not. Read the .cnn directly.
        reader=cat
        case "${cnvkit_pon}" in *.gz) reader=zcat ;; esac
        if "\$reader" "${cnvkit_pon}" | awk 'NR>1 && \$4 ~ /[Aa]ntitarget/ {f=1; exit} END{exit !f}'; then
            echo "ERROR: cnvkit_pon has Antitarget bins (exome-style) but assay=wgs." >&2
            echo "       Build a WGS PoN (cnvkit.py batch --method wgs) or run reference-free." >&2
            exit 1
        fi
    fi

    # ---- fixed input names so output naming is independent of dedup_mode
    #      (_markdup / _consensus / _bqsr all become tumor.bam) ----
    ln -sf ${tumor_bam} tumor.bam
    ln -sf ${tumor_bai} tumor.bam.bai
    ${link_normal}

    cnvkit.py batch tumor.bam \\
        ${batch_args} \\
        --processes ${task.cpus} \\
        --output-dir .

    # CNVkit names outputs after the input basename (tumor.*); canonicalize.
    mv tumor.cns ${meta.sample_id}.cns
    mv tumor.cnr ${meta.sample_id}.cnr
    cnvkit.py call ${meta.sample_id}.cns -o ${meta.sample_id}-call.cns
    """
}

// ============================================================
//  STAGE 5C: MSI - MSIsensor-pro
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
    if (normal_bam) {
        // paired - msi submodule, msisensor_ref = scan list
        """
        set -euo pipefail
        msisensor-pro msi -d ${msisensor_list} \\
            -t ${tumour_bam} -n ${normal_bam} ${region_flag} \\
            -o ${meta.sample_id}_msi -b ${task.cpus} -c 20
        """
    } else {
        // tumor-only - pro submodule, msisensor_ref = pre-built baseline
        """
        set -euo pipefail
        msisensor-pro pro -d ${msisensor_list} \\
            -t ${tumour_bam} ${region_flag} \\
            -o ${meta.sample_id}_msi -b ${task.cpus} -c 20
        """
    }
}


// ============================================================
//  Gene mapping of MSIsensor-pro unstable microsatellite sites.
//  Reads the paired-mode _somatic file; base R (runs in an R
//  container, NOT the MSIsensor image). Keyed on Ensembl gene ID.
// ============================================================
process MSISENSOR_GENE_MAP {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling (analysis-separated) -----------------------
    publishDir { "${pubDir('variant_calling', true)}/msi" }, mode: 'copy'

    input:
    tuple val(meta), path(somatic)     // ${sample}_msi_somatic, extracted at the call site
    path  gene_bed                     // ensembl115_genes.bed (chr,start,end,symbol,ENSG)

    output:
    tuple val(meta), path("${meta.sample_id}_msi_genes.tsv"), emit: msi_genes

    script:
    """
    set -euo pipefail
    cat > msi_gene_map.R <<'RSCRIPT'
    args <- commandArgs(trailingOnly = TRUE)
    sites <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    req <- c("chromosome","location","repeat_times","repeat_unit_bases")
    if (!all(req %in% colnames(sites)))
        stop(sprintf("_somatic missing columns: %s", paste(setdiff(req, colnames(sites)), collapse=",")), call.=FALSE)
    bed <- read.delim(args[2], sep="\\t", header=FALSE, stringsAsFactors=FALSE)
    if (ncol(bed) < 5) stop("gene BED needs 5 columns: chr,start,end,symbol,ENSG", call.=FALSE)
    bed <- bed[, 1:5]; colnames(bed) <- c("chr","start","end","gene","gene_id")
    norm_chr <- function(x) { x <- sub("^chr", "", as.character(x)); x[x %in% c("M","MT")] <- "MT"; x }
    sites\$chromosome <- norm_chr(sites\$chromosome); bed\$chr <- norm_chr(bed\$chr)
    bed\$start      <- bed\$start + 1L
    sites\$location <- suppressWarnings(as.integer(sites\$location))
    have_diff <- "difference" %in% colnames(sites); have_fdr <- "FDR" %in% colnames(sites)
    bbc <- split(seq_len(nrow(bed)), bed\$chr)
    gsym <- character(0); gid <- character(0); si <- integer(0)
    for (i in seq_len(nrow(sites))) {
        idx <- bbc[[ sites\$chromosome[i] ]]; if (is.null(idx)) next
        hit <- which(bed\$start[idx] <= sites\$location[i] & bed\$end[idx] >= sites\$location[i])
        for (j in idx[hit]) { gsym <- c(gsym, bed\$gene[j]); gid <- c(gid, bed\$gene_id[j]); si <- c(si, i) }
    }
    out <- data.frame(
        Gene_ID=gid, Hugo_Symbol=gsym, Chromosome=sites\$chromosome[si], Position=sites\$location[si],
        Repeat_Unit=sites\$repeat_unit_bases[si], Ref_Repeat_Times=sites\$repeat_times[si],
        Repeat_Diff = if (have_diff) sites\$difference[si] else NA,
        FDR         = if (have_fdr)  sites\$FDR[si]        else NA,
        Status="SOMATIC_UNSTABLE", stringsAsFactors=FALSE)
    write.table(out, args[3], sep="\\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("Mapped %d unstable site(s) to %d gene(s)\\n", nrow(out), length(unique(gid))))
    RSCRIPT
    Rscript --vanilla msi_gene_map.R ${somatic} ${gene_bed} ${meta.sample_id}_msi_genes.tsv
    rm -f msi_gene_map.R
    """
}

// ============================================================
//  STAGE 5D MSI Status: MSI status classification
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
    def cutoff    = params.msi_cutoff ?: 20
    def min_sites = (params.msi_min_sites != null) ? params.msi_min_sites : 50
    """
    set -euo pipefail

    main="${meta.sample_id}_msi"
    if [ ! -s "\$main" ]; then
        echo "ERROR: MSIsensor-pro main output '\$main' not found among staged files" >&2
        ls -l >&2
        exit 1
    fi

    score=\$(awk 'NR==2 && \$3!=""{print \$3}' "\$main")
    if [ -z "\${score:-}" ]; then
        score="NA"
        status="NA"
    else
        status=\$(awk -v s="\$score" -v c="${cutoff}" 'BEGIN{ if(s+0>=c+0) print "MSI-H"; else print "MSS" }')
    fi

    total=\$(awk 'NR==2{print \$1}' "\$main")
    somatic=\$(awk 'NR==2{print \$2}' "\$main")
    # Too few evaluable microsatellites -> the score is not interpretable.
    if [ -n "\${total:-}" ] && [ "\$status" != "NA" ] && [ "\${total%.*}" -lt ${min_sites} ]; then
        status="NA_LOW_SITES"
    fi

    {   echo -e "sample_id\\t${meta.sample_id}"
        echo -e "total_sites\\t\${total:-NA}"
        echo -e "somatic_unstable_sites\\t\${somatic:-NA}"
        echo -e "msi_score_pct\\t\$score"
        echo -e "msi_cutoff_pct\\t${cutoff}"
        echo -e "msi_status\\t\$status"
    } > ${meta.sample_id}_msi_status.txt
    """
}


// ============================================================
//  STAGE 6A: ANNOTATION - VEP
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
    def clinvar_flag = clinvar ? "--custom ${clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN" : ''
    // COSMIC is a somatic-mutation catalog - apply only in somatic contexts, not germline.
    def cosmic_flag = (cosmic && params.analysis in ['somatic', 'tumor_only'])
        ? "--custom ${cosmic},COSMIC,vcf,exact,0,GENOME_SCREEN_SAMPLE_COUNT,TIER" : ''
    """
    set -euo pipefail
    syn_flag=""
    if [ -s "${syn_path}" ]; then syn_flag="--synonyms ${syn_path}"; fi

    vep \\
        --input_file  ${vcf} \\
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
        ${clinvar_flag} \\
        ${cosmic_flag} \\
        --stats_file  ${meta.sample_id}_vep_summary.html \\
        --buffer_size 5000 --no_escape

    tabix -p vcf ${meta.sample_id}_vep.vcf.gz
    """
}


// ============================================================
//  STAGE 6B: Tumor Mutational Burden (TMB)
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
    tuple val(meta), path(annotated_vcf), path(annotated_tbi)
    path  interval_bed
    path  coding_bed, stageAs: 'coding_bed/*'   // params.tmb_coding_bed (CDS BED) or [] -> capture-only territory

    output:
    tuple val(meta), path("${meta.sample_id}_tmb.txt"), emit: tmb

    script:
    def cutoff      = params.tmb_cutoff         ?: 10
    def incl_splice = params.tmb_include_splice ?: false
    def terms = [
        'missense_variant','stop_gained','stop_lost','start_lost',
        'frameshift_variant','inframe_insertion','inframe_deletion',
        'protein_altering_variant'
    ]
    if (incl_splice) { terms += ['splice_acceptor_variant','splice_donor_variant'] }
    def term_re = terms.join('|')
    // WGS: interval_bed is the CNVkit access BED, not a coding territory - a coding
    // numerator over it is not a TMB. Emit NA (same schema) until a WGS
    // callable-territory denominator is implemented.
    if (meta.assay == 'wgs') {
    """
    set -euo pipefail
    echo "COMPUTE_TMB: assay=wgs - TMB not computed (no valid WGS denominator yet)" >&2
    {   echo -e "sample_id\\t${meta.sample_id}"
        echo -e "nonsynonymous_coding_variants\\tNA"
        echo -e "covered_megabases\\tNA"
        echo -e "tmb_mut_per_mb\\tNA"
        echo -e "tmb_cutoff\\t${cutoff}"
        echo -e "tmb_status\\tNA"
        echo -e "include_splice\\t${incl_splice}"
        echo -e "territory\\tNA"
        echo -e "consequence_source\\tNA"
    } > ${meta.sample_id}_tmb.txt
    """
    } else {
    """
    set -euo pipefail

    # Merge a BED (plain or .gz) into sorted, non-overlapping BED3 on stdout.
    merge_bed() {
        if echo "\$1" | grep -q '[.]gz\$'; then CAT=zcat; else CAT=cat; fi
        \$CAT "\$1" | grep -v -E '^(#|track|browser)' \\
            | awk -F'\\t' 'NF >= 3 { print \$1 "\\t" \$2 "\\t" \$3 }' \\
            | LC_ALL=C sort -k1,1 -k2,2n \\
            | awk -F'\\t' 'BEGIN { OFS = "\\t" }
                { if (\$1 == c && \$2 <= e) { if (\$3 > e) e = \$3 }
                  else { if (c != "") print c, s, e; c = \$1; s = \$2; e = \$3 } }
                END { if (c != "") print c, s, e }'
    }

    # ---- TERRITORY (denominator AND numerator region): capture, intersected with
    #      the coding BED when provided (sweep over +1/-1 boundary events). ----
    # Contig naming (chr1 vs 1, chrM vs MT): bring every BED to the VCF's naming. A mismatch
    # would otherwise give an empty territory and a silent NA.
    bcftools view -h ${annotated_vcf} | sed -n 's/^##contig=<ID=\\([^,>]*\\).*/\\1/p' > vcf.contigs
    harmonize() {
        if [ ! -s vcf.contigs ]; then cat "\$1"; return; fi
        awk -F'\\t' 'BEGIN { OFS = "\\t" } NR == FNR { ok[\$1] = 1; next }
            { c = \$1
              if (!(c in ok)) { d = (c ~ /^chr/) ? substr(c, 4) : "chr" c
                                if (d == "chrMT") d = "chrM"; else if (d == "M") d = "MT"
                                if (d in ok) c = d; else next }
              \$1 = c; print }' vcf.contigs "\$1" | LC_ALL=C sort -k1,1 -k2,2n
    }
    merge_bed "${interval_bed}" > capture.raw.bed
    harmonize capture.raw.bed > capture.merged.bed
    if [ -n "${coding_bed}" ]; then
        merge_bed "${coding_bed}" > coding.raw.bed
        harmonize coding.raw.bed > coding.merged.bed
        {   awk -F'\\t' '{ print \$1 "\\t" \$2 "\\t1\\tA"; print \$1 "\\t" \$3 "\\t-1\\tA" }' capture.merged.bed
            awk -F'\\t' '{ print \$1 "\\t" \$2 "\\t1\\tB"; print \$1 "\\t" \$3 "\\t-1\\tB" }' coding.merged.bed
        } | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \\
          | awk -F'\\t' 'BEGIN { OFS = "\\t" }
              \$1 != c { c = \$1; a = 0; b = 0; open = 0 }
              { if (\$4 == "A") a += \$3; else b += \$3
                if (!open && a > 0 && b > 0) { st = \$2; open = 1 }
                else if (open && !(a > 0 && b > 0)) { if (\$2 > st) print c, st, \$2; open = 0 } }' \\
          > territory.bed
        territory="capture_x_coding_bed"
    else
        cp capture.merged.bed territory.bed
        territory="capture_bed_only(includes_non-coding_target;set_params.tmb_coding_bed)"
    fi
    if [ ! -s territory.bed ]; then
        echo "ERROR: TMB territory is empty - interval_bed / tmb_coding_bed share no contig (or no overlap) with the VCF." >&2
        echo "       VCF contigs: \$(head -3 vcf.contigs | tr '\\n' ' ')  capture: \$(cut -f1 capture.raw.bed | sort -u | head -3 | tr '\\n' ' ')" >&2
        exit 1
    fi
    covered_mb=\$(awk -F'\\t' '{ t += \$3 - \$2 } END { printf "%.6f", t / 1e6 }' territory.bed)

    # ---- NUMERATOR: PASS nonsynonymous variants inside the territory. One consequence
    #      per variant: VEP PICK transcript (MANE-first pick_order) when annotated,
    #      else worst-across-transcripts (VCFs annotated before PICK was enabled). ----
    csq_fields=\$(bcftools +split-vep -l ${annotated_vcf} | awk '{ print \$2 }')
    if grep -qx PICK <<< "\$csq_fields"; then
        csq_src="VEP_PICK"
        # PICK is per allele+gene: count each VARIANT once if any picked consequence qualifies.
        nonsyn=\$(bcftools view -f PASS,. -T territory.bed -Ou ${annotated_vcf} \\
            | bcftools +split-vep -d -f '%PICK\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%Consequence\\n' - \\
            | awk -F'\\t' -v re='(${term_re})' '\$1 == "1" && \$6 ~ re { print \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 }' \\
            | sort -u | wc -l)
    else
        csq_src="worst_across_transcripts"
        nonsyn=\$(bcftools view -f PASS,. -T territory.bed -Ou ${annotated_vcf} \\
            | bcftools +split-vep -s worst -f '%Consequence\\n' - \\
            | grep -E -c '(${term_re})' || true)
    fi

    # ---- TMB = numerator / denominator ----
    tmb=\$(awk -v n="\$nonsyn" -v mb="\$covered_mb" 'BEGIN{ if(mb>0) printf "%.4f", n/mb; else print "NA" }')
    status=\$(awk -v t="\$tmb" -v c="${cutoff}" 'BEGIN{ if(t=="NA"){print "NA"} else if(t+0>=c+0){print "TMB-High"} else {print "TMB-Low"} }')

    {   echo -e "sample_id\\t${meta.sample_id}"
        echo -e "nonsynonymous_coding_variants\\t\$nonsyn"
        echo -e "covered_megabases\\t\$covered_mb"
        echo -e "tmb_mut_per_mb\\t\$tmb"
        echo -e "tmb_cutoff\\t${cutoff}"
        echo -e "tmb_status\\t\$status"
        echo -e "include_splice\\t${incl_splice}"
        echo -e "territory\\t\$territory"
        echo -e "consequence_source\\t\$csq_src"
    } > ${meta.sample_id}_tmb.txt
    """
    }
}



// ============================================================
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
//  STAGE 7: Variant pre-classification (NOT a final clinical tier)
//
//  Aggregates automatable evidence - gene role/tier (COSMIC Cancer Gene
//  Census), VEP consequence/impact, population frequency, ClinVar (if
//  annotated) - into a computational pre-classification per variant.
//  This is decision SUPPORT for a molecular pathologist, NOT a final
//  AMP/ASCO/CAP tier (which requires a curated drug-variant-disease
//  knowledge base and expert review).
// ============================================================
//         }
//         {
//             chrom=\$1; pos=\$2; ref=\$3; alt=\$4; sym=\$5; csq=\$6; impact=\$7; af=\$8; clin=\$9; cosmic_cnt=\$10
//             in_cgc = (sym in role) ? "yes" : "no"
//             cgc_role = (sym in role) ? role[sym] : "NA"
//             cgc_tier = (sym in tier) ? tier[sym] : "NA"
//             common = (af!="" && af!="." && af+0>afc+0) ? "yes" : "no"
//             recurrent = (cosmic_cnt!="" && cosmic_cnt!="." && cosmic_cnt+0>=cosmic_min+0) ? "yes" : "no"
//             cls="VUS"
//             if(common=="yes"){ cls="LIKELY_BENIGN" }
//             else if(clin ~ /[Pp]athogenic/ && in_cgc=="yes"){ cls="ONCOGENIC_STRONG" }
//             else if(clin ~ /[Bb]enign/){ cls="LIKELY_BENIGN" }
//             else if(in_cgc=="yes"){
//                 if(impact=="HIGH") cls="ONCOGENIC_STRONG"
//                 else if(impact=="MODERATE") cls="ONCOGENIC_MODERATE"
//                 else cls="VUS"
//             } else {
//                 if(impact=="HIGH"||impact=="MODERATE") cls="VUS"
//                 else cls="LIKELY_BENIGN"
//             }
//             if(recurrent=="yes" && in_cgc=="yes" && cls=="VUS"){ cls="ONCOGENIC_MODERATE" }
//             print chrom,pos,ref,alt,sym,csq,impact,(af==""?".":af),(clin==""?".":clin),(cosmic_cnt==""?".":cosmic_cnt),cgc_role,cgc_tier,cls
//         }
//     ' variants_raw.tsv > body.tsv

//     {   echo "# COMPUTATIONAL PRE-CLASSIFICATION FOR PATHOLOGIST REVIEW - NOT A FINAL AMP/ASCO/CAP CLINICAL TIER"
//         echo "# Final tiering requires a curated knowledge base (OncoKB/CIViC) and expert review."
//         echo -e "CHROM\\tPOS\\tREF\\tALT\\tGENE\\tCONSEQUENCE\\tIMPACT\\t${af_field}\\tClinVar\\tCOSMIC_RECURRENCE\\tCGC_ROLE\\tCGC_TIER\\tPRE_CLASSIFICATION"
//         cat body.tsv
//     } > ${meta.sample_id}_preclassification.tsv
//     """
// }


process CLASSIFY_VARIANTS {
    tag { "${meta.sample_id}" }
    label 'process_low'

    publishDir { pubDir('classification', true) }, mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf), path(annotated_tbi)
    path  cgc_tsv

    output:
    tuple val(meta), path("${meta.sample_id}_preclassification.tsv"), emit: classification

    script:
    def af_field          = params.classify_af_field    ?: 'MAX_AF'
    def af_cutoff         = params.classify_common_af    ?: 0.01
    def cosmic_min        = params.cosmic_recurrence_min ?: 10
    def clinvar_min_stars = (params.clinvar_min_stars != null) ? params.clinvar_min_stars : 1   // 0 = gating off (elvis would turn 0 into 1)
    def cls_mode          = (params.analysis == 'germline') ? 'germline' : 'somatic'
    def ba1_af            = (params.germline_ba1_af != null) ? params.germline_ba1_af : 0.05
    """
    set -euo pipefail

    # CGC gene -> role/tier, columns located by HEADER NAME (layout differs between
    # CGC releases; fixed indices fail silently). Missing CGC -> empty lookup.
    if [ -n "${cgc_tsv}" ]; then
        if echo "${cgc_tsv}" | grep -q '[.]gz\$'; then CATC=zcat; else CATC=cat; fi
        "\$CATC" "${cgc_tsv}" | awk -F'\\t' '
            NR==1 { for (i = 1; i <= NF; i++) { h = toupper(\$i); gsub(/[ ]+/, "_", h)
                        if (h == "GENE_SYMBOL")    g = i
                        if (h == "ROLE_IN_CANCER") r = i
                        if (h == "TIER")           t = i }
                    if (!g || !r || !t) { print "ERROR: CGC header lacks GENE_SYMBOL/ROLE_IN_CANCER/TIER" > "/dev/stderr"; bad = 1; exit 1 }
                    next }
            \$g != "" { print \$g "\\t" \$r "\\t" \$t }
            END { if (bad) exit 1 }' > cgc_lookup.tsv
    else
        : > cgc_lookup.tsv
    fi

    fields="\$(bcftools +split-vep -l "${annotated_vcf}" 2>/dev/null | awk '{print \$2}')"

    if grep -qx "COSMIC_GENOME_SCREEN_SAMPLE_COUNT" <<< "\$fields"; then
        COSMIC_FMT='\\t%COSMIC_GENOME_SCREEN_SAMPLE_COUNT'
    else
        COSMIC_FMT=''
    fi

    if grep -qx "ClinVar_CLNSIG" <<< "\$fields"; then
        CLIN_FMT='%ClinVar_CLNSIG'
    else
        echo "WARN: ClinVar_CLNSIG not in VCF; using VEP cache CLIN_SIG" >&2
        CLIN_FMT='%CLIN_SIG'
    fi

    if grep -qx "ClinVar_CLNREVSTAT" <<< "\$fields"; then
        REVSTAT_FMT='%ClinVar_CLNREVSTAT'
    else
        REVSTAT_FMT='.'
    fi

    # PASS-only: consensus/filter run in 'tag' mode, so failed calls are still in
    # the VCF; they must not be pre-classified or reported.
    # One consequence per variant: the VEP PICK transcript (MANE-first pick_order) when
    # present, else worst-across-transcripts (VCFs annotated before PICK was enabled).
    CSQ_FMT="%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%IMPACT\\t%${af_field}\\t\${CLIN_FMT}\\t\${REVSTAT_FMT}\${COSMIC_FMT}\\n"
    if grep -qx "PICK" <<< "\$fields"; then
        bcftools view -f PASS,. -Ou "${annotated_vcf}" \\
        | bcftools +split-vep -d -f "%PICK\\t\${CSQ_FMT}" - \\
        | awk -F'\\t' '\$1 == "1"' | cut -f2- | sort -u > variants_raw.tsv
    else
        echo "WARN: no VEP PICK field; using worst consequence across transcripts" >&2
        bcftools view -f PASS,. -Ou "${annotated_vcf}" \\
        | bcftools +split-vep -s worst -d -f "\${CSQ_FMT}" - \\
        | sort -u > variants_raw.tsv
    fi

    cat > classify.awk <<'AWK'
BEGIN {
    OFS = "\t"
    while ((getline line < "cgc_lookup.tsv") > 0) {
        split(line, a, "\t"); role[a[1]] = a[2]; tier[a[1]] = a[3]
    }
}
function clinvar_class(s,   t, i, n, tok, hp, hb) {
    if (s == "" || s == ".") return "other"
    s = tolower(s)
    if (s ~ /conflicting/) return "conflicting"
    gsub(/[/,&|]+/, " ", s)
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
function clinvar_stars(s) {
    if (s == "" || s == ".")            return 0
    s = tolower(s)
    if (s ~ /no_assertion/)             return 0
    if (s ~ /no_classification/)        return 0
    if (s ~ /practice_guideline/)       return 4
    if (s ~ /reviewed_by_expert_panel/) return 3
    if (s ~ /multiple_submitters/)      return 2
    if (s ~ /criteria_provided/)        return 1
    return 0
}
{
    chrom=\$1; pos=\$2; ref=\$3; alt=\$4; gene=\$5; csq=\$6; impact=\$7
    af=\$8; clin=\$9; revstat=\$10; cosmic=\$11
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
    lof = (csq ~ /(stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_lost)/)
    tsg = (tolower(cgc_role) ~ /tsg/)
    onc = (tolower(cgc_role) ~ /oncogene/)
    if (MODE == "germline") {
        # ACMG-style vocabulary. ClinVar (star-gated) is a GERMLINE pathogenicity source.
        if      (cv == "pathogenic")                cls = "PATHOGENIC_OR_LIKELY"
        else if (cv == "benign")                    cls = "BENIGN_OR_LIKELY"
        else if (af != "" && af != "." && af+0 > BA1+0) cls = "BENIGN_OR_LIKELY"          # BA1-like
        else if (is_common)                         cls = "LIKELY_BENIGN"
        else if (lof && in_cgc)                     cls = "LOF_CANCER_GENE_REVIEW"
        else if (high || moderate)                  cls = "VUS"
        else                                        cls = "LIKELY_BENIGN"
    } else {
        # Somatic ONCOGENICITY is gene-role aware; ClinVar (germline assertions) does not drive it.
        #   TSG + loss-of-function         -> oncogenic (strong if COSMIC-recurrent)
        #   oncogene + recurrent non-LoF   -> oncogenic (hotspot-like)
        #   ClinVar P/LP (star-gated) in a TSG -> supporting (loss-of-function mechanism transfers)
        #   truncating an oncogene alone   -> NOT oncogenic evidence (VUS_CANCER_GENE)
        if      (cv == "benign" || is_common)       cls = "LIKELY_BENIGN"                  # population / germline leakage
        else if (tsg && lof)                        cls = is_recurrent ? "ONCOGENIC_STRONG" : "ONCOGENIC_MODERATE"
        else if (tsg && cv == "pathogenic")         cls = "ONCOGENIC_MODERATE"   # germline-P/LP in a TSG = loss-of-function evidence
        else if (onc && !lof && moderate && is_recurrent) cls = "ONCOGENIC_STRONG"
        else if (in_cgc && is_recurrent && (high || moderate)) cls = "ONCOGENIC_MODERATE"
        else if (in_cgc && (high || moderate))      cls = "VUS_CANCER_GENE"
        else if (high || moderate)                  cls = "VUS"
        else                                        cls = "LIKELY_BENIGN"
    }
    print chrom, pos, ref, alt, gene, csq, impact, (af==""?".":af), (clin==""?".":clin), stars, (cosmic==""?".":cosmic), cgc_role, cgc_tier, cls
}
AWK

    awk -F'\\t' -v AF_CUTOFF="${af_cutoff}" -v COSMIC_MIN="${cosmic_min}" -v MIN_STARS="${clinvar_min_stars}" \\
        -v MODE="${cls_mode}" -v BA1="${ba1_af}" \\
        -f classify.awk variants_raw.tsv > body.tsv

    {
        echo "# COMPUTATIONAL PRE-CLASSIFICATION FOR PATHOLOGIST REVIEW - NOT A FINAL AMP/ASCO/CAP CLINICAL TIER"
        echo "# Final tiering requires a curated knowledge base (OncoKB/CIViC) and expert review."
        echo "# ClinVar source: ClinVar_CLNSIG when present, else VEP cache CLIN_SIG."
        echo "# Rules (${cls_mode}): somatic = CGC role-aware (TSG loss-of-function; recurrent oncogene change); germline = ClinVar P/LP (star-gated) + BA1 AF > ${ba1_af}."
        printf 'CHROM\\tPOS\\tREF\\tALT\\tGENE\\tCONSEQUENCE\\tIMPACT\\t%s\\tClinVar\\tCLINVAR_STARS\\tCOSMIC_RECURRENCE\\tCGC_ROLE\\tCGC_TIER\\tPRE_CLASSIFICATION\\n' "${af_field}"
        cat body.tsv
    } > ${meta.sample_id}_preclassification.tsv
    """
}

// ============================================================
//  STAGE 7: VCF to MAF generation 
// ============================================================
// process VCF2MAF {
//     tag { "${meta.sample_id}" }
//     label 'process_tiny'

//     // -- 07 classification (analysis-separated)----------------
//     publishDir { "${pubDir('classification', true)}/maf" }, mode: 'copy'

//     input:
//     tuple val(meta), path(vcf), path(tbi)
//     path  ref_fasta
//     path  ref_fai

//     output:
//     tuple val(meta), path("${meta.sample_id}.maf"), emit: maf

//    script:
//     def tumor_id  = "${meta.sample_id}"
//     def normal_id = "${meta.sample_id}".replaceAll(/_T$/, '_N')
//     """
//     set -euo pipefail

//     zcat ${vcf} > input.vcf

//     vcf2maf.pl \\
//         --input-vcf      input.vcf \\
//         --output-maf     full.maf \\
//         --ref-fasta      ${ref_fasta} \\
//         --tumor-id       ${tumor_id} \\
//         --normal-id      ${normal_id} \\
//         --vcf-tumor-id   ${tumor_id} \\
//         --vcf-normal-id  ${normal_id} \\
//         --ncbi-build     GRCh38 \\
//         --inhibit-vep \\
//         --retain-fmt     GT,AD,AF,DP

//     cat > trim_maf.py <<'PYEOF'
// import csv, sys
// want = ["Hugo_Symbol","Gene","dbSNP_RS","NCBI_Build","Chromosome",
//         "Start_Position","End_Position","Strand","Variant_Classification",
//         "Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
//         "HGVSc","HGVSp","HGVSp_Short","Transcript_ID","Exon_Number",
//         "t_GT","n_GT","t_depth","t_ref_count","t_alt_count",
//         "n_depth","n_ref_count","n_alt_count","Allele","Protein_position",
//         "Amino_acids","Codons","EXON","INTRON","Tumor_VAF","Normal_VAF",
//         "all_effects","Feature"]
// computed = {"Tumor_VAF", "Normal_VAF"}
// with open("full.maf") as f:
//     lines = f.readlines()
// comments = [l for l in lines if l.startswith("#")]
// data     = [l for l in lines if not l.startswith("#")]
// reader = csv.DictReader(data, delimiter="\\t")
// have = set(reader.fieldnames or [])
// miss = [c for c in want if c not in computed and c not in have]
// if miss:
//     sys.stderr.write("WARN missing columns: " + ",".join(miss) + "\\n")
// def vaf(a, d):
//     try:
//         a, d = float(a), float(d)
//         return round(a/d, 4) if d > 0 else ""
//     except (ValueError, TypeError):
//         return ""
// with open(sys.argv[1], "w", newline="") as out:
//     for c in comments:
//         out.write(c)
//     w = csv.DictWriter(out, fieldnames=want, delimiter="\\t", extrasaction="ignore")
//     w.writeheader()
//     for row in reader:
//         r = {c: row.get(c, "") for c in want if c not in computed}
//         r["Tumor_VAF"]  = vaf(row.get("t_alt_count"), row.get("t_depth"))
//         r["Normal_VAF"] = vaf(row.get("n_alt_count"), row.get("n_depth"))
//         w.writerow(r)
// PYEOF

//     PY=\$(command -v python3 || echo /opt/conda/bin/python3)
//     "\$PY" trim_maf.py ${meta.sample_id}.maf
//     """
// }

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
    def is_germline = (params.analysis == 'germline')
    def tumor_id  = "${meta.sample_id}"
    def normal_id = is_germline ? "${meta.sample_id}" : "${meta.sample_id}".replaceAll(/_T$/, '_N')

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
        --retain-fmt     GT,AD,AF,DP

    cat > trim_maf.py <<'PYEOF'

import csv, sys
germline = ${is_germline ? 'True' : 'False'}
if germline:
    want = ["Hugo_Symbol","Entrez_Gene_Id","Gene","dbSNP_RS","NCBI_Build","Chromosome",
            "Start_Position","End_Position","Strand","Variant_Classification","Variant_Type",
            "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
            "Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","HGVSc","HGVSp","HGVSp_Short",
            "Transcript_ID","Exon_Number","t_GT","n_GT","n_depth","n_ref_count","n_alt_count",
            "Allele","Protein_position","Amino_acids","Codons","EXON","INTRON","Normal_VAF",
            "all_effects","Feature","FILTER"]
    computed = {"Normal_VAF"}
else:
    want = ["Hugo_Symbol","Gene","dbSNP_RS","NCBI_Build","Chromosome",
            "Start_Position","End_Position","Strand","Variant_Classification",
            "Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
            "HGVSc","HGVSp","HGVSp_Short","Transcript_ID","Exon_Number",
            "t_GT","n_GT","t_depth","t_ref_count","t_alt_count",
            "n_depth","n_ref_count","n_alt_count","Allele","Protein_position",
            "Amino_acids","Codons","EXON","INTRON","Tumor_VAF","Normal_VAF",
            "all_effects","Feature","FILTER"]
    computed = {"Tumor_VAF", "Normal_VAF"}

with open("full.maf") as f:
    lines = f.readlines()
comments = [l for l in lines if l.startswith("#")]
data     = [l for l in lines if not l.startswith("#")]
reader = csv.DictReader(data, delimiter="\\t")
have = set(reader.fieldnames or [])
miss = [c for c in want if c not in computed and c not in have]
if miss:
    sys.stderr.write("WARN missing columns: " + ",".join(miss) + "\\n")
def vaf(a, d):
    try:
        a, d = float(a), float(d)
        return round(a/d, 4) if d > 0 else ""
    except (ValueError, TypeError):
        return ""
with open(sys.argv[1], "w", newline="") as out:
    for c in comments:
        out.write(c)
    w = csv.DictWriter(out, fieldnames=want, delimiter="\t", extrasaction="ignore")
    w.writeheader()

    for row in reader:
        r = {c: row.get(c, "") for c in want if c not in computed}
        if not germline:
            r["Tumor_VAF"] = vaf(row.get("t_alt_count"), row.get("t_depth"))
        r["Normal_VAF"] = vaf(row.get("n_alt_count"), row.get("n_depth"))
        w.writerow(r)
PYEOF

    PY=\$(command -v python3 || echo /opt/conda/bin/python3)
    "\$PY" trim_maf.py ${meta.sample_id}.maf
    """
}

// ============================================================
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
    set -uo pipefail
    awk -F'\\t' '\$1 == "@SQ" { n = ""; l = ""
        for (i = 2; i <= NF; i++) { if (\$i ~ /^SN:/) n = substr(\$i, 4); if (\$i ~ /^LN:/) l = substr(\$i, 4) }
        print n "\\t" l }' ${ref_dict} > dict.tsv
    fail=0; : > preflight_references.txt
    note() { echo "\$1" | tee -a preflight_references.txt >&2; }
    # exact contig lookup in the .dict (reads the whole file: no early pipe close -> no SIGPIPE)
    in_dict() { awk -F'\\t' -v c="\$1" '\$1 == c { f = 1 } END { exit !f }' dict.tsv; }
    # NOTE: Nextflow runs task scripts with 'bash -e'; together with pipefail, any
    # 'producer | head' on a large file exits 141 (SIGPIPE) - hence 'set +o pipefail' subshells.

    for v in vcf/*; do
        [ -e "\$v" ] || continue
        name=\$(basename "\$v")
        hdr=\$(bcftools view -h "\$v" 2>/dev/null) || { note "FAIL  \$name: not readable by bcftools"; fail=1; continue; }
        # build check: PRIMARY-chromosome lengths in the header must equal the reference .dict
        # (alt/decoy/unplaced contigs are ignored to avoid false alarms)
        bad=\$(printf '%s\\n' "\$hdr" | sed -n 's/^##contig=<ID=\\([^,>]*\\),length=\\([0-9]*\\).*/\\1\\t\\2/p' \\
              | awk -F'\\t' 'NR == FNR { L[\$1] = \$2; next }
                  \$1 ~ /^(chr)?([0-9]+|X|Y)\$/ && (\$1 in L) && L[\$1] != \$2 { printf "%s(%s vs ref %s) ", \$1, \$2, L[\$1] }' dict.tsv -)
        first=\$(set +o pipefail; bcftools view -H "\$v" 2>/dev/null | head -1 | cut -f1)
        if [ -n "\$bad" ]; then note "FAIL  \$name: contig lengths differ from the reference (different build?): \$bad"; fail=1
        elif [ -n "\$first" ] && ! in_dict "\$first"; then
            case "\$name" in
              *clinvar*|*ClinVar*|*cosmic*|*COSMIC*) note "WARN  \$name: contig '\$first' not in reference (OK only if VEP --synonyms maps it)";;
              *) note "FAIL  \$name: contig '\$first' not in the reference .dict (chr naming / build mismatch)"; fail=1;;
            esac
        else note "OK    \$name"; fi
    done

    for b in bed/*; do
        [ -e "\$b" ] || continue
        name=\$(basename "\$b")
        if echo "\$name" | grep -q '[.]gz\$'; then CAT=zcat; else CAT=cat; fi
        first=\$(set +o pipefail; \$CAT "\$b" | grep -v -E '^(#|track|browser)' | head -1 | cut -f1)
        if [ -z "\$first" ]; then note "FAIL  \$name: empty"; fail=1
        elif ! in_dict "\$first"; then note "FAIL  \$name: contig '\$first' not in the reference .dict"; fail=1
        else note "OK    \$name"; fi
    done

    [ "\$fail" -eq 0 ] || { echo "PREFLIGHT_REFERENCES failed - see messages above" >&2; exit 1; }
    """
}
