#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
//  workflows_rna.nf  -  RNA Fusion Sub-Workflow + RNA channel factory
//
//  Stage gating uses the SAME mechanism as the DNA workflow: a stageRuns('stage')
//  predicate built on the shared stageRankMap() (single source of truth in
//  utils_v3), so no gate hard-codes a rank literal.
//
//  Rank map (utils_v3.validStages order):
//    1 mapping | 2 markdup | 3 metrics | 4 bqsr | 5 variant_calling
//    6 annotation | 7 classification
//  The RNA arm exercises 1,3,5,6,7 (markdup=2 is QC-only inside metrics;
//  bqsr=4 is DNA-only and never run here).
//
//  Flow (gold-standard triad, matching nf-core/rnafusion):
//    STAR (1 shared run) -> Arriba + STAR-Fusion(kickstart) + FusionCatcher
//    -> consensus (N-of-M + known-fusion rescue) -> FusionInspector validation
//    -> annotation (known-fusion DB + COSMIC CGC) -> per-patient report.
// ============================================================

//  ------------------- stage 0 -------------------------------
//  ------------------------ utils ----------------------------
include { GENERATE_RNA_META_INFO      } from './utils'
include { stageRankMap; stageRuns     } from './utils'

//  ------------------------ utils ----------------------------
include { toFileChannel; toOptionalFileChannel; stripChunkMeta } from './workflows_dna'

//  ---------------- components_common ------------------------
include { MANIFEST_ROW                } from './components_common'
include { WRITE_MANIFEST              } from './components_common'

//  --------------- multistage processes ----------------------
// Reused, assay-agnostic processes from the DNA arm (single source of truth).
include { MULTIQC                     } from './components_common'
include { CRAM_CONVERT                } from './components_common'
include { SAMTOOLS_STATS              } from './components_common'


// RNA-specific processes.
include {
    STAR_ALIGN; RNASEQ_METRICS;
    ARRIBA; STARFUSION; FUSIONCATCHER;
    FUSION_CONSENSUS; FILTER_FUSION_CANDIDATES; FUSIONINSPECTOR;
    ANNOTATE_FUSIONS; CONSOLIDATE_FUSION_REPORT
} from './components_rna'



// ============================================================
//  RNA_WORKFLOW  -  entry sub-workflow (called from main_v3.nf)
// ============================================================
workflow RNA_WORKFLOW {

    take:
    samplesheet

    main:
    // -- stage gating: integer ranks, same mechanism as the DNA workflow --
    // int start_rank = stageRankMap()[params.start_at]
    // int stop_rank  = stageRankMap()[params.stop_at]

    // -- reference channels (real file OR [] for absent optionals) --------
    // One reference for the whole RNA arm: the CTAT genome + annotation the STAR index was
    // built from (STAR, Arriba, CRAM, stats, gene counts). rna_use_ctat_reference=false
    // restores params.ref_fasta / ref_fai / gtf.
    def use_ctat = params.rna_use_ctat_reference
    ch_ref_fasta = toFileChannel(use_ctat ? "${params.ctat_lib}/ref_genome.fa"     : params.ref_fasta)
    ch_ref_fai   = toFileChannel(use_ctat ? "${params.ctat_lib}/ref_genome.fa.fai" : params.ref_fai)
    ch_gtf       = toFileChannel(use_ctat ? "${params.ctat_lib}/ref_annot.gtf"     : params.gtf)


    // In shared mode the STAR index MUST be the CTAT lib's index (kickstart
    // consistency); default to it unless params.star_index is set explicitly.
    ch_ctat_lib    = toFileChannel(params.ctat_lib)
    ch_star_index  = params.star_index ? toFileChannel(params.star_index): toFileChannel("${params.ctat_lib}/ref_genome.fa.star.idx")

    ch_arriba_bl = params.arriba_blacklist        ? toFileChannel(params.arriba_blacklist)        : channel.value([])
    ch_arriba_kf = params.arriba_known_fusions    ? toFileChannel(params.arriba_known_fusions)    : channel.value([])
    ch_arriba_pd = params.arriba_protein_domains  ? toFileChannel(params.arriba_protein_domains)  : channel.value([])

    ch_fc_db     = toOptionalFileChannel(params.fusioncatcher_db)
    ch_known_fus = toOptionalFileChannel(params.known_fusions_tsv)
    ch_cgc_tsv   = toOptionalFileChannel(params.cosmic_cgc_tsv)
    ch_refflat   = toOptionalFileChannel(params.refflat)
    ch_rrna      = toOptionalFileChannel(params.rrna_interval_list)

    // -- read ingestion + per-sample grouping -----------------------------
    GENERATE_RNA_META_INFO(samplesheet)
    ch_manifest = GENERATE_RNA_META_INFO.out.reads

    if (params.write_manifest) {
            // One row per sample (input file is the first read / BAM), then consolidate.
            ch_manifest_input = ch_manifest.map { meta, reads ->
                // 'reads' may be a list (paired FASTQ) or a single file; take the first.
                def f = (reads instanceof List) ? reads[0] : reads
                tuple(meta, f)
            }
            MANIFEST_ROW(ch_manifest_input)
            WRITE_MANIFEST(MANIFEST_ROW.out.row.collect())
    }
        

    // -- Entry-point: Per-sample reads (unconditional) ----------------------
    // Grouped from the chunk-level manifest. Needed by BOTH the align path and
    // the re-entry path (which reads meta to locate staged BAMs), so it is NOT
    // gated behind the mapping stage - gating it left it unbound on re-entry.
    ch_reads = ch_manifest.map { meta, reads -> tuple(meta.sample_id, meta, reads) }
        .groupTuple()
        .map { sid, metas, reads_list ->
            def m0          = metas[0]
            def sample_meta = m0.findAll { k, v -> !(k in ['uid', 'rg']) }
            def r1 = reads_list.collect { it[0] }
            def r2 = sample_meta.single_end ? [] : reads_list.collect { it[1] }
            tuple(sample_meta, r1, r2)
        }

    // --------------------------------------------------------------------
    //  STAGE 1 (rank 1) - STAR alignment, OR stage prior BAM/CRAM for re-entry
    //  (same stageRuns('mapping') gate + cram-preferred idiom as the DNA workflow)
    // --------------------------------------------------------------------
    if (stageRuns('mapping')) {
        STAR_ALIGN(ch_reads, ch_star_index, ch_gtf)
        ch_bam      = STAR_ALIGN.out.bam           // tuple(meta, bam, bai)
        ch_chimeric = STAR_ALIGN.out.chimeric      // tuple(meta, junction)
        ch_star_log = STAR_ALIGN.out.log
    }
    else {
        // re-entry (start_rank >= 3): stage STAR output from input_bam_dir,
        // preferring the CRAM archive if present (as the DNA workflow does).
        ch_bam = ch_reads.map { meta, _r1, _r2 ->
            def cram = file("${params.input_bam_dir}/${meta.sample_id}.cram")
            def bam  = file("${params.input_bam_dir}/${meta.sample_id}.Aligned.sortedByCoord.out.bam")
            def active_bam = cram.exists() ? cram : bam
            def active_idx = cram.exists()
                ? file("${params.input_bam_dir}/${meta.sample_id}.cram.crai", checkIfExists: true)
                : file("${params.input_bam_dir}/${meta.sample_id}.Aligned.sortedByCoord.out.bam.bai", checkIfExists: true)
            return tuple(meta, active_bam, active_idx)
        }
        ch_chimeric = ch_reads.map { meta, _r1, _r2 ->
            def j = file("${params.input_bam_dir}/${meta.sample_id}.Chimeric.out.junction")
            tuple(meta, j.exists() ? j : [])
        }
        ch_star_log = channel.empty()
    }

    // disk: archive the (large) STAR BAM as CRAM - only when we actually aligned
    if (stageRuns('mapping') && params.rna_cram_convert) {
        CRAM_CONVERT('rna', ch_bam, ch_ref_fasta, ch_ref_fai)
    }

    // --------------------------------------------------------------------
    //  STAGE 3 (rank 3) - RNA QC  (MarkDuplicates QC-only, CollectRnaSeqMetrics)
    // --------------------------------------------------------------------
    if (stageRuns('metrics')) {
        // STAR gene counts (this run only) drive strandedness inference; [] on re-entry.
        ch_counts = stageRuns('mapping')
            ? STAR_ALIGN.out.gene_counts.map { m, f -> tuple(m.sample_id, f) }
            : channel.empty()
        ch_rnaqc_in = ch_bam
            .map { m, b, i -> tuple(m.sample_id, m, b, i) }
            .join(ch_counts, by: 0, remainder: true)
            .filter { row -> row[1] != null }
            .map { row -> tuple(row[1], row[2], row[3], (row.size() > 4 && row[4]) ? row[4] : []) }
        RNASEQ_METRICS(ch_rnaqc_in, ch_refflat, ch_rrna)
        SAMTOOLS_STATS(ch_bam, ch_ref_fasta, ch_ref_fai)

        ch_qc = ch_star_log
            .mix(RNASEQ_METRICS.out.metrics.map { meta, f -> f })
            .mix(SAMTOOLS_STATS.out.stats.map    { meta, f -> f })
            .collect()
            .map { files -> tuple('rna_fusion', files) }
        MULTIQC(ch_qc)
    }

    // --------------------------------------------------------------------
    //  STAGE 5 (rank 5) - fusion calling -> consensus -> in-silico validation
    // --------------------------------------------------------------------
    if (stageRuns('variant_calling')) {

        // ---- caller 1: Arriba (reads the shared BAM) ----
        if (params.use_arriba) {
            ARRIBA(ch_bam, ch_ref_fasta, ch_gtf, ch_arriba_bl, ch_arriba_kf, ch_arriba_pd)
            ch_arriba = ARRIBA.out.fusions
        } else {
            ch_arriba = ch_bam.map { meta, bam, bai -> tuple(meta, []) }
        }

        // ---- caller 2: STAR-Fusion (kickstart from junction, or full mode) ----
        if (params.use_starfusion) {
            if (params.rna_shared_star != false) {
                ch_sf_in = ch_chimeric.map { meta, j -> tuple(meta, j, [], []) }
            } else {
                ch_sf_in = ch_reads.map { meta, r1, r2 -> tuple(meta, [], r1, r2) }
            }
            STARFUSION(ch_sf_in, ch_ctat_lib)
            ch_starf = STARFUSION.out.fusions
        } else {
            ch_starf = ch_bam.map { meta, bam, bai -> tuple(meta, []) }
        }

        // ---- caller 3: FusionCatcher (PE only; SE<130bp unsupported upstream) ----
        if (params.use_fusioncatcher) {
            ch_fc_branch = ch_reads.branch {
                run:  !it[0].single_end || params.fusioncatcher_allow_se
                skip: true
            }
            FUSIONCATCHER(ch_fc_branch.run.map { meta, r1, r2 -> tuple(meta, r1, r2) }, ch_fc_db)
            ch_fc = FUSIONCATCHER.out.fusions
                .mix(ch_fc_branch.skip.map { meta, r1, r2 -> tuple(meta, []) })
        } else {
            ch_fc = ch_bam.map { meta, bam, bai -> tuple(meta, []) }
        }

        // ---- consensus: join the three callers on sample_id, carry meta from BAM ----
        ch_consensus_in = ch_bam.map { meta, bam, bai -> tuple(meta.sample_id, meta) }
            .join(ch_arriba.map { meta, f -> tuple(meta.sample_id, f) })
            .join(ch_starf.map  { meta, f -> tuple(meta.sample_id, f) })
            .join(ch_fc.map     { meta, f -> tuple(meta.sample_id, f) })
            .map { sid, meta, a, s, f -> tuple(meta, a, s, f) }

        FUSION_CONSENSUS(ch_consensus_in, ch_known_fus)

        // ---- in-silico validation: FusionInspector on the PASS_* candidate list ----
        // if (params.use_fusioninspector) {
        //     ch_fi_in = FUSION_CONSENSUS.out.candidates
        //         .map { meta, c -> tuple(meta.sample_id, meta, c) }
        //         .join(ch_reads.map { meta, r1, r2 -> tuple(meta.sample_id, r1, r2) })
        //         .map { sid, meta, c, r1, r2 -> tuple(meta, c, r1, r2) }
        //     FUSIONINSPECTOR(ch_fi_in, ch_ctat_lib)
        //     ch_validated = FUSIONINSPECTOR.out.validated
        // } else {
        //     ch_validated = FUSION_CONSENSUS.out.consensus.map { meta, c -> tuple(meta, []) }
        // }

        if (params.use_fusioninspector) {
            // 1. Filter consensus candidates to drop unannotated GTF genes
            FILTER_FUSION_CANDIDATES(
                FUSION_CONSENSUS.out.candidates,
                ch_ctat_lib
            )

            // 2. Join filtered candidates back with ch_reads on sample_id
            ch_fi_in = FILTER_FUSION_CANDIDATES.out.filtered
                .map { meta, c -> tuple(meta.sample_id, meta, c) }
                .join(ch_reads.map { meta, r1, r2 -> tuple(meta.sample_id, r1, r2) })
                .map { sid, meta, c, r1, r2 -> tuple(meta, c, r1, r2) }

            // 3. Run FusionInspector on clean, filtered candidates
            FUSIONINSPECTOR(ch_fi_in, ch_ctat_lib)
            ch_validated = FUSIONINSPECTOR.out.validated
        } else {
            ch_validated = FUSION_CONSENSUS.out.consensus.map { meta, c -> tuple(meta, []) }
        }

    }

    // --------------------------------------------------------------------
    //  STAGE 6 (rank 6) - annotation (known-fusion DB + COSMIC CGC tiering)
    // --------------------------------------------------------------------
    if (stageRuns('annotation')) {
        ch_annot_in = FUSION_CONSENSUS.out.consensus
            .map { meta, c -> tuple(meta.sample_id, meta, c) }
            .join(ch_validated.map { meta, v -> tuple(meta.sample_id, v) })
            .map { sid, meta, c, v -> tuple(meta, c, v) }
        ANNOTATE_FUSIONS(ch_annot_in, ch_cgc_tsv, ch_known_fus)
    }

    // --------------------------------------------------------------------
    //  STAGE 7 (rank 7) - per-patient consolidated fusion report
    // --------------------------------------------------------------------
    if (stageRuns('classification')) {
        ch_report_in = ANNOTATE_FUSIONS.out.annotated
            .map { meta, a -> tuple(meta.sample_id, meta, a) }
            .join(FUSION_CONSENSUS.out.stats.map { meta, s -> tuple(meta.sample_id, s) })
            .map { sid, meta, a, s -> tuple(meta.patient_id, meta.sample_id, a, s) }
        CONSOLIDATE_FUSION_REPORT(ch_report_in)
    }
}
