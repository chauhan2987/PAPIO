#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Shared JVM-options helper (single source of truth in utils), same import
// idiom as components_common / _somatic / _germline.
include { gatkJavaOpts; pubDir; rnaStrandedness } from './utils'

// ============================================================
//  components_rna.nf  -  RNA fusion process definitions
//
//  Stage map (same --start_at/--stop_at vocabulary as the DNA arm):
//    mapping         - STAR_ALIGN, SALMON_QUANT
//    metrics         - RNASEQ_METRICS, samtools stats
//    variant_calling - ARRIBA, STARFUSION, FUSIONCATCHER -> FUSION_CONSENSUS
//                      -> FILTER_FUSION_CANDIDATES -> FUSIONINSPECTOR;
//                      RNA_IMMUNE_PREP -> TRUST4, TELESCOPE
//    annotation      - ANNOTATE_FUSIONS
//    classification  - CONSOLIDATE_FUSION_REPORT
//
//  One STAR run feeds two callers: its coordinate-sorted BAM goes to Arriba and
//  its Chimeric.out.junction to STAR-Fusion ("kickstart" mode, no second STAR
//  run; params.rna_shared_star = false gives STAR-Fusion its own STAR run).
//  FusionCatcher aligns the raw FASTQs itself. The STAR index must use the same
//  genome and chromosome naming as the CTAT library; its GTF may be a newer
//  GENCODE release, because FUSION_CONSENSUS matches callers on Ensembl gene IDs.
//  Fusion callers read the raw STAR BAM / FASTQs; no RNA duplicate marking.
// ============================================================


// ============================================================
//  STAGE 1 mapping: STAR_ALIGN - STAR alignment (Arriba BAM and STAR-Fusion junctions)
//
//  Multi-lane/chunk safe: STAR accepts comma-separated R1/R2 file lists, so
//  no upstream merge step is needed (unlike the DNA markdup merge). Emits:
//    - <sid>.Aligned.sortedByCoord.out.bam (+ .bai)  - Arriba, QC, CRAM
//    - <sid>.Chimeric.out.junction                   - STAR-Fusion kickstart
//    - <sid>.SJ.out.tab, <sid>.Log.final.out         - QC/MultiQC
//
//  Chimeric parameters follow Arriba's published recommendations AND add
//  `Junctions` + `--chimOutJunctionFormat 1` so the same run feeds STAR-Fusion.
// ============================================================

process STAR_ALIGN {
    tag { "${meta.sample_id}" }
    label 'process_star'

    maxRetries 2

    // -- 01 mapping -----------------------------------------------------
    publishDir { "${pubDir('mapping', true)}/star" }, mode: 'copy',
        pattern: '*{Log.final.out,SJ.out.tab,Chimeric.out.junction,ReadsPerGene.out.tab}'

    input:
    tuple val(meta), path(reads1), path(reads2)
    path  star_index        // built with the GTF (splice junctions + GeneCounts come from the index)

    output:
    tuple val(meta),
          path("${meta.sample_id}.Aligned.sortedByCoord.out.bam"),
          path("${meta.sample_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sample_id}.Chimeric.out.junction"),  emit: chimeric
    path  "${meta.sample_id}.Log.final.out",                           emit: log
    path  "${meta.sample_id}.SJ.out.tab",                              emit: sj
    // Gene counts (cols 2/3/4 = unstranded / fwd / rev strand) from the GTF STAR loaded.
    tuple val(meta), path("${meta.sample_id}.ReadsPerGene.out.tab"),  emit: gene_counts
    // Transcriptome alignments for Salmon (salmon_transcripts set); work directory only, not published.
    tuple val(meta), path("${meta.sample_id}.Aligned.toTranscriptome.out.bam"), optional: true, emit: tx_bam


    script:
    // Reserve 2 threads for the concurrent samtools sort/index; STAR gets the rest.
    def star_threads = Math.max(task.cpus - 2, 4)
    def r1_list      = (reads1 instanceof List ? reads1 : [reads1]).join(',')
    def r2_list      = (reads2 && (reads2 instanceof List ? reads2.size() : 1))
                       ? (reads2 instanceof List ? reads2 : [reads2]).join(',') : ''
    def read_files   = r2_list ? "${r1_list} ${r2_list}" : "${r1_list}"
    // gzip FASTQ -> zcat; plain FASTQ -> read directly (zcat on plain text fails).
    def r1_first     = (reads1 instanceof List ? reads1[0] : reads1).toString()
    def rf_cmd       = r1_first.endsWith('.gz') ? '--readFilesCommand zcat' : ''
    def twopass      = (params.star_twopass == false) ? 'None' : 'Basic'
    def rg_line      = "ID:${meta.sample_id} SM:${meta.sample_id} LB:${meta.batch} PL:ILLUMINA"
    // Salmon input: transcriptome BAM keeping alignments with indels and soft-clips (soft-clips
    // extended); only pairs aligned as pairs (BanSingleEnd). STAR >= 2.7.11a option.
    def quant_mode   = params.salmon_transcripts
                       ? '--quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeSAMoutput BanSingleEnd_ExtendSoftclip'
                       : '--quantMode GeneCounts'

    // -- Chimeric-detection strategy (VERIFIED against STAR issue #1288 + Arriba wiki) --
    //   STAR's NEW chimeric detection (--chimMultimapNmax > 0) ONLY works with
    //   --chimOutType Junctions and is INCOMPATIBLE with WithinBAM. So:
    //   - shared mode  - --chimOutType WithinBAM HardClip Junctions with multimap
    //     OMITTED (old detection). One run feeds Arriba (WithinBAM) AND STAR-Fusion
    //     kickstart (Junctions) - the STAR-maintainer-endorsed way to get both.
    //   - separate mode (params.rna_shared_star=false) - WithinBAM HardClip only,
    //     with Arriba's recommended --chimMultimapNmax 50 for max multimap
    //     sensitivity; STAR-Fusion then runs its own bundled STAR (full mode).
    def shared       = (params.rna_shared_star != false)
    def chim_seg_min = params.star_chim_segment_min ?: 10
    def chim_out     = shared ? 'WithinBAM HardClip Junctions' : 'WithinBAM HardClip'
    // Mate-overlap merging (--peOverlapNbasesMin, Arriba's recommendation) is accepted by STAR
    // only when chimeric output is WithinBAM alone; with Junctions also written (shared mode)
    // STAR stops with a PARAMETERS error, so shared mode runs without it.
    def pe_overlap   = shared ? '' : '--peOverlapNbasesMin 10'
    def chim_mm      = (!shared && (params.star_chim_multimap_nmax ?: 0) as int > 0)
                       ? "--chimMultimapNmax ${params.star_chim_multimap_nmax}" : ''
    """
    set -euo pipefail

    tmp_dir="\$PWD/${meta.sample_id}_star_tmp"
    # STAR requires --outTmpDir to NOT pre-exist; it creates and manages it.
    rm -rf "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    STAR \\
        --runThreadN ${star_threads} \\
        --genomeDir ${star_index} --genomeLoad NoSharedMemory \\
        --readFilesIn ${read_files} ${rf_cmd} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within --outBAMcompression 6 \\
        --outSAMattrRGline ${rg_line} \\
        --twopassMode ${twopass} \\
        ${quant_mode} \\
        --outFilterMultimapNmax 50 ${pe_overlap} \\
        --alignSplicedMateMapLminOverLmate 0.5 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentMin ${chim_seg_min} \\
        --chimOutType ${chim_out} \\
        --chimOutJunctionFormat 1 \\
        --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \\
        --chimSegmentReadGapMax 3 ${chim_mm} \\
        --outFileNamePrefix ${meta.sample_id}. \\
        --outTmpDir "\$tmp_dir"

    # STAR names the junction 'Chimeric.out.junction' unless multimap>0, in which
    # case it is written next to the prefix; normalise the name defensively.
    if [ ! -f "${meta.sample_id}.Chimeric.out.junction" ]; then
        if [ -f "Chimeric.out.junction" ]; then
            mv Chimeric.out.junction "${meta.sample_id}.Chimeric.out.junction"
        else
            # produce an empty (header-only) junction so downstream stays graceful
            printf '# empty chimeric junction (no chimeric reads)\\n' > "${meta.sample_id}.Chimeric.out.junction"
        fi
    fi

    samtools index -@ ${task.cpus} ${meta.sample_id}.Aligned.sortedByCoord.out.bam
    """
}


// ============================================================
//  STAGE 1 mapping: SALMON_QUANT - transcript quantification - Salmon (expression only)
//
//  Alignment mode on STAR's Aligned.toTranscriptome.out.bam (as nf-core/rnaseq
//  star_salmon): no second pass over the reads and no Salmon index (a few GB of
//  memory instead of ~20 GB). The transcript FASTA (salmon_transcripts) must be built
//  from the same gtf + genome as the STAR index (bin/gtf_to_transcripts.sh): Salmon
//  stops on a transcript of the BAM missing from it. Library type inferred (-l A).
//  Fragment-GC bias needs pairs.
// ============================================================

process SALMON_QUANT {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // -- 01 mapping -----------------------------------------------------
    publishDir { "${pubDir('mapping', true)}/salmon" }, mode: 'copy'

    input:
    tuple val(meta), path(tx_bam)
    path  transcripts

    output:
    tuple val(meta), path("${meta.sample_id}"), emit: results

    script:
    def gc_flag = meta.single_end ? '' : '--gcBias'
    """
    set -euo pipefail

    salmon --no-version-check quant \\
        -t ${transcripts} \\
        -l A \\
        -a ${tx_bam} \\
        --seqBias ${gc_flag} \\
        -p ${task.cpus} \\
        -o ${meta.sample_id}
    """
}


// ============================================================
//  STAGE 3 metrics: RNASEQ_METRICS - RNA-seq metrics and strandedness
//
//  Captures the RNA-specific QC MultiQC understands: %coding / %UTR /
//  %intronic / %intergenic bases, median 5'-3' bias, %rRNA. Needs a
//  refFlat (params.refflat) and, ideally, an rRNA interval_list
//  (params.rrna_interval_list) - the latter is optional (skipped if null).
// ============================================================

process RNASEQ_METRICS {
    tag { "${meta.sample_id}" }
    label 'process_qc'

    scratch false

    memory { 16.GB * task.attempt }
    time   { 60.min * task.attempt }

    // -- 03 metrics / QC ------------------------------------------------
    publishDir { "${pubDir('metrics', true)}/rnaseq_qc" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(gene_counts)   // STAR ReadsPerGene.out.tab or [] (re-entry)
    path  refflat
    path  rrna_interval_list

    output:
    tuple val(meta), path("${meta.sample_id}.rna_metrics"),      emit: metrics, optional: true   // needs refflat
    tuple val(meta), path("${meta.sample_id}.strandedness.txt"), emit: strandedness

    script:
    def java_opts = gatkJavaOpts(task)
    // STRAND per library prep. 'auto' = infer from STAR gene counts (col3 = read-1 sense,
    // col4 = read-1 antisense); without counts (re-entry) 'auto' falls back to NONE.
    def strand = rnaStrandedness()
    // same vocabulary as ARRIBA (yes = forward-stranded, as Arriba's -s yes)
    def strand_spec = strand in ['reverse','rf','fr-firststrand'] ? 'SECOND_READ_TRANSCRIPTION_STRAND'
                    : strand in ['forward','yes','fr','fr-secondstrand'] ? 'FIRST_READ_TRANSCRIPTION_STRAND'
                    : 'NONE'
    def strand_arg = (strand == 'auto' && gene_counts) ? "auto --gene-counts ${gene_counts}" : strand_spec
    // Paired-end: keep properly paired reads (Picard stops on interchromosomal discordant
    // pairs). Single-end reads never carry the proper-pair flag, so they are used as is.
    def qc_filter = meta.single_end ? '-F 0x900' : '-f 2'
    def refflat_arg = refflat ? "--refflat ${refflat}" : ''
    def rrna_arg    = rrna_interval_list ? "--rrna ${rrna_interval_list}" : ''
    """
    rnaseq_metrics.sh --sample ${meta.sample_id} --bam ${bam} \\
        --strand ${strand_arg} --basis "rna_strandedness=${strand}" \\
        ${refflat_arg} ${rrna_arg} --qc-filter "${qc_filter}" \\
        --threads ${task.cpus} --java-opts "${java_opts}"
    """
}


// ============================================================
//  STAGE 5 variant_calling: ARRIBA - fusion caller - Arriba
//
//  Reads the shared STAR BAM directly. Arriba's own blacklist +
//  known-fusions + protein-domains resources ship with the release and are
//  passed as params; blacklist is effectively required for precision.
// ============================================================

process ARRIBA {
    tag { "${meta.sample_id}" }
    label 'process_arriba'

    maxRetries 2

    // -- 05 variant_calling / fusion ------------------------------------
    publishDir { "${pubDir('variant_calling', true)}/arriba" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  gtf
    path  arriba_blacklist
    path  arriba_known_fusions
    path  arriba_protein_domains

    output:
    tuple val(meta), path("${meta.sample_id}_arriba_fusions.tsv"),           emit: fusions
    path  "${meta.sample_id}_arriba_fusions.discarded.tsv",                  emit: discarded

    script:
    def bl_flag  = arriba_blacklist       ? "-b ${arriba_blacklist}"       : '-f blacklist'   // -f blacklist disables only if truly absent
    def kf_flag  = arriba_known_fusions   ? "-k ${arriba_known_fusions}"    : ''
    def pd_flag  = arriba_protein_domains ? "-p ${arriba_protein_domains}"  : ''
    // Arriba's -s accepts ONLY: auto | yes | no | reverse. Translate common
    // library-prep terms; default to 'auto' (robust) for anything unrecognised.
    def strand   = rnaStrandedness()
    def strand_map = [ 'auto':'auto', 'unstranded':'no', 'no':'no',
                       'forward':'yes', 'yes':'yes', 'fr':'yes', 'fr-secondstrand':'yes',
                       'reverse':'reverse', 'rf':'reverse', 'fr-firststrand':'reverse' ]
    def strand_flag = "-s ${strand_map.getOrDefault(strand, 'auto')}"
    """
    set -euo pipefail

    arriba \\
        -x ${bam} \\
        -o ${meta.sample_id}_arriba_fusions.tsv \\
        -O ${meta.sample_id}_arriba_fusions.discarded.tsv \\
        -a ${ref_fasta} \\
        -g ${gtf} \\
        ${bl_flag} ${kf_flag} ${pd_flag} ${strand_flag}
    """
}


// ============================================================
//  STAGE 5 variant_calling: STARFUSION - fusion caller - STAR-Fusion
//
//  Two modes:
//    kickstart (default, params.rna_shared_star=true): consume the shared
//      STAR Chimeric.out.junction - NO second STAR run.
//    full: run STAR-Fusion's own bundled STAR on the raw FASTQs (fidelity
//      fallback). Selected when a junction is not supplied (empty channel).
//  Requires the CTAT genome library (params.ctat_lib).
// ============================================================

process STARFUSION {
    tag { "${meta.sample_id}" }
    label 'process_star'

    maxRetries 2

    // -- 05 variant_calling / fusion ------------------------------------
    publishDir { "${pubDir('variant_calling', true)}/starfusion" }, mode: 'copy'

    input:
    tuple val(meta), path(junction), path(reads1), path(reads2)
    path  ctat_lib

    output:
    tuple val(meta), path("${meta.sample_id}_starfusion.tsv"), emit: fusions

    script:
    // Optional inputs arrive as a real file OR as [] (empty list = falsy),
    // the same idiom the DNA arm uses for pon_flag/clinvar_flag. junction is a
    // real file in kickstart mode and [] in full (FASTQ) mode.
    def has_junction = junction as boolean
    def r1_list = reads1 ? (reads1 instanceof List ? reads1 : [reads1]).join(',') : ''
    def r2_list = reads2 ? (reads2 instanceof List ? reads2 : [reads2]).join(',') : ''
    def fq_args = r2_list ? "--left_fq ${r1_list} --right_fq ${r2_list}"
                : r1_list ? "--left_fq ${r1_list}" : ''
    def mode_args = has_junction ? "-J ${junction}" : fq_args
    """
    set -euo pipefail

    tmp_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/${meta.sample_id}_sf_tmp.XXXXXX")
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    STAR-Fusion \\
        --genome_lib_dir ${ctat_lib} \\
        ${mode_args} \\
        --CPU ${task.cpus} \\
        --output_dir starfusion_out \\
        --tmpdir "\$tmp_dir"

    # Prefer the abridged prediction table; fall back to the full one; guarantee
    # a (possibly header-only) output so consensus never sees a missing file.
    if   [ -s starfusion_out/star-fusion.fusion_predictions.abridged.tsv ]; then
        cp starfusion_out/star-fusion.fusion_predictions.abridged.tsv ${meta.sample_id}_starfusion.tsv
    elif [ -s starfusion_out/star-fusion.fusion_predictions.tsv ]; then
        cp starfusion_out/star-fusion.fusion_predictions.tsv ${meta.sample_id}_starfusion.tsv
    else
        printf '#FusionName\\tJunctionReadCount\\tSpanningFragCount\\tLeftGene\\tLeftBreakpoint\\tRightGene\\tRightBreakpoint\\n' \\
            > ${meta.sample_id}_starfusion.tsv
    fi
    """
}


// ============================================================
//  STAGE 5 variant_calling: FUSIONCATCHER - fusion caller - FusionCatcher
//
//  Runs its own aligners on the raw FASTQs with a curated DB
//  (params.fusioncatcher_db). Slow (hours) - gate with
//  params.use_fusioncatcher. Its orthogonal alignment strategy is what makes
//  the ensemble robust, so keep it on when turnaround allows.
// ============================================================

process FUSIONCATCHER {
    tag { "${meta.sample_id}" }
    label 'process_fusioncatcher'

    scratch false

    maxRetries 1

    // -- 05 variant_calling / fusion ------------------------------------
    publishDir { "${pubDir('variant_calling', true)}/fusioncatcher" }, mode: 'copy'

    input:
    tuple val(meta), path(reads1), path(reads2)
    path  fusioncatcher_db

    output:
    tuple val(meta), path("${meta.sample_id}_fusioncatcher.txt"), emit: fusions

    script:
    def r1_files = (reads1 instanceof List ? reads1 : [reads1]).join(' ')
    def r2_files = reads2 ? (reads2 instanceof List ? reads2 : [reads2]).join(' ') : ''
    def fq_ext   = (reads1 instanceof List ? reads1[0] : reads1).toString().endsWith('.gz') ? 'fastq.gz' : 'fastq'

    """
    set -euo pipefail

    tmp_dir="\$PWD/fc_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir" fc_input fc_out' EXIT

    # lanes concatenated per mate; the extension follows the input (gzip or plain FASTQ)
    mkdir -p fc_input
    cat ${r1_files} > fc_input/${meta.sample_id}_R1.${fq_ext}
    ${ reads2 ? "cat ${r2_files} > fc_input/${meta.sample_id}_R2.${fq_ext}" : "" }

    fusioncatcher.py \\
        -d ${fusioncatcher_db} \\
        -i fc_input \\
        -o fc_out \\
        -p ${task.cpus} \\
        --tmp "\$tmp_dir" \\
        --skip-blat

    if [ -s fc_out/final-list_candidate-fusion-genes.txt ]; then
        cp fc_out/final-list_candidate-fusion-genes.txt ${meta.sample_id}_fusioncatcher.txt
    else
        printf 'Gene_1_symbol(5end_fusion_partner)\\tGene_2_symbol(3end_fusion_partner)\\tFusion_point_for_gene_1(5end_fusion_partner)\\tFusion_point_for_gene_2(3end_fusion_partner)\\tSpanning_pairs\\tSpanning_unique_reads\\n' \\
            > ${meta.sample_id}_fusioncatcher.txt
    fi

    rm -rf fc_input fc_out "\$tmp_dir"
    """
}


// ============================================================
//  STAGE 5 variant_calling: FUSION_CONSENSUS - fusion consensus across callers
//
//  Normalises the three callers' heterogeneous outputs to a canonical
//  (gene5'--gene3') key + chrom:pos breakpoints, tallies caller support, and
//  emits a consensus verdict:
//    - PASS_CONSENSUS  : called by >= params.fusion_min_callers callers, OR
//    - PASS_KNOWN      : a known oncogenic pair (params.known_fusions_tsv)
//                        supported by >=1 caller (single-caller rescue), OR
//    - REVIEW          : otherwise (kept, flagged).
//  Read-through / adjacent-gene calls (Arriba 'type', or same-chrom within
//  params.fusion_readthrough_kb) are flagged, mirroring the COSMIC-recurrence
//  rescue logic in CLASSIFY_VARIANTS.
// ============================================================

process FUSION_CONSENSUS {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 05 variant_calling / fusion ------------------------------------
    publishDir { "${pubDir('variant_calling', true)}/consensus" }, mode: 'copy'

    input:
    tuple val(meta),
          path(arriba,        stageAs: 'arriba.tsv'),
          path(starfusion,    stageAs: 'starfusion.tsv'),
          path(fusioncatcher, stageAs: 'fusioncatcher.txt')
    path  known_fusions_tsv, stageAs: 'known_ref/*'   // staged in a subfolder: no clash with the working copy's name
    path  gtf,               stageAs: 'gene_ref/*'    // RNA arm GTF: reported gene symbols
    path  ctat_lib,          stageAs: 'ctat_ref/*'    // CTAT library ([] when not set): its GTF names genes for FusionInspector

    output:
    tuple val(meta), path("${meta.sample_id}_fusion_consensus.tsv"),  emit: consensus
    tuple val(meta), path("${meta.sample_id}_fusion_candidates.txt"), emit: candidates
    tuple val(meta), path("${meta.sample_id}_fusion_names.tsv"),      emit: names
    tuple val(meta), path("${meta.sample_id}_fusion_consensus_stats.txt"), emit: stats

    script:
    def min_callers  = params.fusion_min_callers      ?: 2
    def readthru_kb  = params.fusion_readthrough_kb   ?: 200
    def known_arg    = known_fusions_tsv ? "${known_fusions_tsv}" : 'NONE'
    def ctat_gtf     = ctat_lib ? "${ctat_lib}/ref_annot.gtf" : ''
    """
    fusion_consensus.sh ${meta.sample_id} ${min_callers} ${readthru_kb} ${known_arg} ${gtf} ${ctat_gtf}
    """
}


// ============================================================
//  STAGE 5 variant_calling: FILTER_FUSION_CANDIDATES - FusionInspector candidates present in the CTAT annotation
//
//  Filter unknown genes to avoid FusionInspector crash
// ============================================================

process FILTER_FUSION_CANDIDATES {
    tag "${meta.sample_id}"
    label 'process_low'

    input:
    tuple val(meta), path(candidates)
    path  ctat_lib

    output:
    tuple val(meta), path("${meta.sample_id}_filtered_candidates.txt"), emit: filtered

    script:
    """
    filter_fusion_candidates.sh ${candidates} ${ctat_lib} ${meta.sample_id}_filtered_candidates.txt
    """
}


// ============================================================
//  STAGE 5 variant_calling: FUSIONINSPECTOR - in-silico fusion validation - FusionInspector
//
//  Re-aligns reads to the consensus candidate fusion contigs and reports
//  per-fusion supporting evidence (junction/spanning reads, FFPM), the
//  strongest orthogonal confirmation short of RT-PCR. Only the PASS_*
//  candidate list is validated (bounded runtime). Skips gracefully when the
//  candidate list is empty.
// ============================================================

process FUSIONINSPECTOR {
    tag { "${meta.sample_id}" }
    label 'process_fusion_validate'

    maxRetries 2

    // -- 05 variant_calling / fusion ------------------------------------
    publishDir { "${pubDir('variant_calling', true)}/fusioninspector" }, mode: 'copy'

    input:
    tuple val(meta), path(candidates), path(reads1), path(reads2)
    path  ctat_lib

    output:
    tuple val(meta), path("${meta.sample_id}_fusioninspector.tsv"), emit: validated

    script:
    def r1_list = (reads1 instanceof List ? reads1 : [reads1]).join(',')
    def r2_list = reads2 ? (reads2 instanceof List ? reads2 : [reads2]).join(',') : ''
    def fq_args = r2_list ? "--left_fq ${r1_list} --right_fq ${r2_list}" : "--left_fq ${r1_list}"

    """
    set -euo pipefail

    export TMPDIR="\$PWD/fi_tmp"
    mkdir -p "\$TMPDIR"
    trap 'rm -rf "\$TMPDIR"' EXIT

    # Empty candidate list - emit a header-only result and exit 0
    if [ ! -s ${candidates} ]; then
        printf '#FusionName\\tJunctionReadCount\\tSpanningFragCount\\tFFPM\\tannots\\n' \\
            > ${meta.sample_id}_fusioninspector.tsv
        exit 0
    fi

    FusionInspector \\
        --fusions ${candidates} \\
        --genome_lib_dir ${ctat_lib} \\
        ${fq_args} \\
        --out_prefix FI \\
        --vis \\
        --CPU ${task.cpus} \\
        -O fi_out

    # FusionInspector names its table <out_prefix>.FusionInspector.fusions.abridged.tsv
    # (here FI.*). The first match wins; the older names are kept as fallbacks.
    fi_tsv=""
    for cand in fi_out/FI.FusionInspector.fusions.abridged.tsv \\
                fi_out/FI.fusion_predictions.final.abridged.tsv \\
                fi_out/finspector.FusionInspector.fusions.abridged.tsv; do
        if [ -s "\$cand" ]; then fi_tsv="\$cand"; break; fi
    done
    if [ -n "\$fi_tsv" ]; then
        cp "\$fi_tsv" ${meta.sample_id}_fusioninspector.tsv
    else
        echo "WARN: no FusionInspector abridged table found in fi_out/ - writing header-only result" >&2
        ls fi_out >&2 || true
        printf '#FusionName\\tJunctionReadCount\\tSpanningFragCount\\tFFPM\\tannots\\n' \\
            > ${meta.sample_id}_fusioninspector.tsv
    fi

    rm -rf "\$TMPDIR"
    """
}


// ============================================================
//  STAGE 5 variant_calling: RNA_IMMUNE_PREP - BAM preparation for TRUST4 and Telescope
//
//  STAR BAM or archived CRAM in; coordinate BAM (TRUST4) and collated BAM (Telescope) out.
// ============================================================

process RNA_IMMUNE_PREP {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    input:
    tuple val(meta), path(aln), path(idx)
    path  ref_fasta
    path  ref_fai

    output:
    tuple val(meta), path("${meta.sample_id}.coord.bam"),    emit: coord
    tuple val(meta), path("${meta.sample_id}.collated.bam"), emit: collated

    script:
    def is_cram = aln.name.endsWith('.cram')
    """
    set -euo pipefail

    ${ is_cram
        ? "samtools view -@ ${task.cpus} -T ${ref_fasta} -b -o ${meta.sample_id}.coord.bam ${aln}"
        : "ln -s ${aln} ${meta.sample_id}.coord.bam" }

    samtools collate -@ ${task.cpus} --reference ${ref_fasta} \\
        -o ${meta.sample_id}.collated.bam ${aln} \$PWD/collate_tmp
    """
}


// ============================================================
//  STAGE 5 variant_calling: TRUST4 - TCR/BCR repertoire - TRUST4
// ============================================================

process TRUST4 {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    memory { 16.GB * task.attempt }
    time   { 9.h   * task.attempt }

    // -- 05 variant_calling ---------------------------------------------
    publishDir { "${pubDir('variant_calling', true)}/trust4" }, mode: 'copy',
        pattern: '*{_report.tsv,_airr.tsv}'

    input:
    tuple val(meta), path(bam)
    path  trust4_ref_dir

    output:
    tuple val(meta), path("${meta.sample_id}_report.tsv"), emit: report
    path  "${meta.sample_id}_airr.tsv", optional: true,     emit: airr

    script:
    """
    set -euo pipefail

    run-trust4 \\
        -b ${bam} \\
        -f ${trust4_ref_dir}/hg38_bcrtcr.fa \\
        --ref ${trust4_ref_dir}/human_IMGT+C.fa \\
        --abnormalUnmapFlag \\
        -t ${task.cpus} \\
        -o ${meta.sample_id}
    """
}


// ============================================================
//  STAGE 5 variant_calling: TELESCOPE - HERV expression - Telescope
// ============================================================

process TELESCOPE {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // -- 05 variant_calling ---------------------------------------------
    publishDir { "${pubDir('variant_calling', true)}/telescope" }, mode: 'copy',
        pattern: '*-telescope_report.tsv'

    input:
    tuple val(meta), path(collated_bam)
    path  herv_gtf

    output:
    tuple val(meta), path("${meta.sample_id}-telescope_report.tsv"), emit: report

    script:
    """
    set -euo pipefail

    telescope assign \\
        --outdir . \\
        --exp_tag ${meta.sample_id} \\
        --theta_prior 200000 --max_iter 200 \\
        ${collated_bam} ${herv_gtf}
    """
}


// ============================================================
//  STAGE 6 annotation: ANNOTATE_FUSIONS - fusion annotation and tiers
//
//  Joins consensus fusions (+ FusionInspector confirmation) to a known-fusion
//  catalog (params.known_fusions_tsv) and the COSMIC Cancer Gene Census
//  (params.cosmic_cgc_tsv, reused from the DNA arm) to assign a computational
//  tier - the fusion analog of CLASSIFY_VARIANTS. NOT a final clinical call.
// ============================================================

process ANNOTATE_FUSIONS {
    tag { "${meta.sample_id}" }
    label 'process_low'

    // -- 06 annotation --------------------------------------------------
    publishDir { "${pubDir('annotation', true)}" }, mode: 'copy'


    input:
    tuple val(meta),
          path(consensus,       stageAs: 'consensus.tsv'),
          path(fusioninspector, stageAs: 'fi.tsv'),
          path(names,           stageAs: 'names.tsv')
    path  cgc_tsv,           stageAs: 'cgc_ref/*'
    path  known_fusions_tsv, stageAs: 'known_ref/*'

    output:
    tuple val(meta), path("${meta.sample_id}_fusions_annotated.tsv"), emit: annotated

    script:
    def cgc_arg    = cgc_tsv ? "cgc.tsv" : 'NONE'
    def known_arg  = known_fusions_tsv ? "${known_fusions_tsv}" : 'NONE'
    // Resolve .gz-ness in Groovy (avoids fragile shell filename tests).
    def cgc_prep   = !cgc_tsv ? '' :
                     (cgc_tsv.name.endsWith('.gz') ? "zcat ${cgc_tsv} > cgc.tsv" : "cp ${cgc_tsv} cgc.tsv")
    def min_known  = (params.fusion_known_min_reads != null) ? params.fusion_known_min_reads : 3
    """
    set -euo pipefail
    ${cgc_prep}


    PY=\$(command -v python3 || echo /opt/conda/bin/python3)
    "\$PY" "\$(command -v annotate_fusions.py)" ${cgc_arg} ${known_arg} ${meta.sample_id} ${min_known}
    """
}


// ============================================================
//  STAGE 7 classification: CONSOLIDATE_FUSION_REPORT - consolidated fusion report per patient
//
//  (analog of CONSOLIDATE_PATIENT_REPORT)
// ============================================================

process CONSOLIDATE_FUSION_REPORT {
    tag { patient_id }
    label 'process_low'

    // -- 07 classification ----------------------------------------------
    publishDir { "${pubDir('classification', true)}/reports" }, mode: 'copy'

    input:
    tuple val(patient_id), val(sample_id),
          path(annotated, stageAs: 'annotated.tsv'),
          path(stats,     stageAs: 'consensus_stats.txt')

    output:
    tuple val(patient_id), path("${patient_id}_fusion_report.tsv"), emit: report

    script:
    """
    fusion_report.sh ${patient_id} ${sample_id} consensus_stats.txt annotated.tsv
    """
}
