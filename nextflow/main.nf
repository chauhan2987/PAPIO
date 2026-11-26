#!/usr/bin/env nextflow

// ============================================================
//  main.nf  -  Pipeline Entry Point   (STRICT DSL2 )
//
//  Strict-syntax rules honoured here:
//    - NO executable statements at top level. Only includes, function
//      declarations, the entry workflow, and completion handlers.
//    - All parameter DEFAULTS live in nextflow.config (params {} block),
//      NOT as top-level `params.x = ...` assignments.
//    - Constants and validation logic are FUNCTIONS, called from workflow{}.
//
//  Responsibilities: validate params, log a run summary,
//  route to the DNA workflow.
// ============================================================

include { DNA_WORKFLOW } from './workflows_dna'
include { RNA_WORKFLOW } from './workflows_rna'   // RNA fusion arm

// Shared constants live in utils (single source of truth) so validation
// here and execution in DNA_WORKFLOW can never disagree about stage order.
include { stageRankMap; validStages; validAssays;
          validAnalyses; validCallers; validDedupModes; rnaStrandedness } from './utils'


// ============================================================
//  CONSTANTS - thin local alias over the shared stage-rank map.
//  (Every existing call site uses stageRank()[...]; keeping the alias
//  means none of validateParams() below needs to change.)
// ============================================================
def stageRank() { stageRankMap() }


// ============================================================
//  HELPERS - proper function declarations (allowed at top level)
// ============================================================

// does the run touch a given stage? (start_rank <= rank <= stop_rank)
def runsStage(String stageName, int startRank, int stopRank) {
    def r = stageRank()[stageName]
    return (startRank <= r && r <= stopRank)
}

// assert a list of params are set, with context
def requireParams(List names, String reason) {
    names.each { p ->
        if (!params.containsKey(p) || !params[p]) {
            error "ERROR: params.${p} is required ${reason}."
        }
    }
}

// warn (not fail) for soft requirements
def warnParams(List names, String reason) {
    names.each { p ->
        if (!params[p]) {
            log.warn "params.${p} not set - ${reason}"
        }
    }
}


// ============================================================
//  VALIDATION - STAGE-AWARE. Returns [startRank, stopRank].
//  Called once at the top of workflow{} (strict: no top-level run).
// ============================================================
def validateParams() {
    ['germline_low_vaf', 'germline_ba1_af', 'classify_common_af'].each { k ->
        def v = params[k]
        if (v != null && !(v instanceof Number && v > 0 && v < 1))
            error "params.${k} must be a number between 0 and 1 (got '${v}')"
    }
    if (params.germline_mosaic_min_dp != null && !(params.germline_mosaic_min_dp instanceof Integer && params.germline_mosaic_min_dp >= 1))
        error "params.germline_mosaic_min_dp must be a whole number >= 1 (got '${params.germline_mosaic_min_dp}')"
    if (params.germline_ba1_af != null && params.classify_common_af != null && params.germline_ba1_af < params.classify_common_af)
        error "params.germline_ba1_af (${params.germline_ba1_af}) must be >= classify_common_af (${params.classify_common_af}) (BA1 >= BS1)"
    if (params.classify_variants && !params.gnomad_faf_vcf)
        log.warn "gnomad_faf_vcf not set: classification uses VEP MAX_AF (includes founder/small populations); see README Appendix A, step 7.6"
    // -- Enum validation ---------------------------------------
    if (!validDedupModes().contains(params.dedup_mode))
        error "ERROR: Invalid params.dedup_mode '${params.dedup_mode}'. One of: ${validDedupModes().join(', ')}"
    if (!validAssays().contains(params.assay))
        error "ERROR: Invalid params.assay '${params.assay}'. One of: ${validAssays().join(', ')}"
    if (!validAnalyses().contains(params.analysis))
        error "ERROR: Invalid params.analysis '${params.analysis}'. One of: ${validAnalyses().join(', ')}"
    if (!validStages().contains(params.start_at))
        error "ERROR: Invalid params.start_at '${params.start_at}'. One of: ${validStages().join(', ')}"
    if (!validStages().contains(params.stop_at))
        error "ERROR: Invalid params.stop_at '${params.stop_at}'. One of: ${validStages().join(', ')}"

    // Mode switches: a misspelt value must stop the run, not fall back silently (a misspelt
    // consensus_filter_mode would drop sub-threshold calls instead of tagging them).
    def enumParams = [
        somatic_filter_mode  : ['tag', 'remove'],
        consensus_filter_mode: ['tag', 'remove'],
        consensus_backbone   : ['mutect2', 'union'],
        joint_genotyper      : ['classic', 'gnarly'],
        nv_tensor_type       : ['reference', 'read_tensor'],
        manifest_checksum    : ['none', 'header', 'full'],
        rna_strandedness     : ['auto', 'yes', 'no', 'unstranded', 'forward', 'reverse', 'fr', 'rf',
                                'fr-firststrand', 'fr-secondstrand'],
    ]
    enumParams.each { k, allowed ->
        def v = (k == 'rna_strandedness') ? rnaStrandedness() : params[k]
        if (v != null && !(v.toString().toLowerCase() in allowed))
            error "ERROR: Invalid params.${k} '${v}'. One of: ${allowed.join(', ')}"
    }

    def startRank = stageRank()[params.start_at]
    def stopRank  = stageRank()[params.stop_at]
    if (startRank > stopRank)
        error "ERROR: params.start_at '${params.start_at}' cannot come after params.stop_at '${params.stop_at}'"

    // bind_paths: non-empty list; no entry may be a filesystem root
    if (!params.bind_paths || !(params.bind_paths instanceof List) || params.bind_paths.isEmpty())
        error "ERROR: params.bind_paths must be a non-empty LIST of directories the containers can " +
              "read/write (e.g. ['/home/groups/g/u', '/scratch/groups/g/u'])."
    params.bind_paths.each { bp ->
        if (java.nio.file.Paths.get(bp.toString()).normalize().toString() in ['/', '/home', '/scratch', '/data'])
            error "ERROR: bind_paths entry '${bp}' is a filesystem root. Bind your OWNED directory " +
                  "(e.g. /home/groups/<grp>/<user>), not a top-level mount."
    }

    def sheet = params.samplesheet ?:
        (params.assay == 'rna' ? params.rna_samplesheet : params.dna_samplesheet)
    if (!sheet)
        error "No samplesheet for assay='${params.assay}'"

    // -- Guards against SILENT misconfiguration (each of these has bitten a real run) --
    // (1) '--flag false' on the command line arrives as the STRING "false" (truthy).
    def boolText = params.findAll { k, v -> v instanceof String && v.toLowerCase() in ['true', 'false'] }
    if (boolText)
        error "ERROR: parameter(s) ${boolText.keySet().join(', ')} hold the TEXT 'true'/'false' " +
              "(command-line booleans arrive as strings, and 'false' then counts as TRUE). Set them in the params YAML."
    // (2) params-file keys the pipeline does not declare are silently ignored (e.g. 'clinvar' vs 'clinvar_vcf').
    def pfMatch = ((workflow.commandLine ?: '') =~ /-params-file[= ]+(\S+)/)
    def cfgFile = file("${projectDir}/nextflow.config")
    if (pfMatch.find() && cfgFile.exists()) {
        def pfArg = pfMatch.group(1)
        def pf    = file(pfArg.startsWith('/') ? pfArg : "${workflow.launchDir}/${pfArg}")
        if (pf.exists() && pf.name ==~ /.*\.ya?ml$/) {
            def declared = [] as Set
            def inParams = false
            cfgFile.readLines().each { l ->
                if (l ==~ /^params\s*\{.*/) { inParams = true }
                else if (inParams && l ==~ /^\}.*/) { inParams = false }
                else if (inParams) { def m = (l =~ /^\s+([A-Za-z_]\w*)\s*=/); if (m.find()) declared << m.group(1) }
            }
            def given   = pf.readLines().collect { l -> def m = (l =~ /^([A-Za-z_]\w*)\s*:/); m.find() ? m.group(1) : null }.findAll { k -> k }
            def unknown = given.findAll { k -> !(k in declared) }
            if (declared && unknown)
                log.warn "Params file '${pf.name}': key(s) NOT declared by the pipeline and IGNORED: ${unknown.join(', ')} " +
                         "(typo or renamed parameter?)"
        }
    }
    // (3) params.samplesheet silently overrides dna_/rna_samplesheet.
    if (params.samplesheet && (params.dna_samplesheet || params.rna_samplesheet))
        log.warn "params.samplesheet='${params.samplesheet}' OVERRIDES dna_samplesheet/rna_samplesheet - using it."
    // (4) Values nextflow.config reads at parse time must come from the params file, not a later -c config.
    if (workflow.workDir?.toString()?.contains('null/nf_work'))
        error "ERROR: workDir is '${workflow.workDir}' - params.work_root was unset when nextflow.config was read. " +
              "Set work_root in the params file (a later -c config is read too late)."
    def procCfg = workflow.session?.config?.navigate('process')
    def nullCtr = (procCfg instanceof Map) ? procCfg.findAll { k, v ->
        v instanceof Map && v.container != null && v.container.toString().startsWith('null/') }.keySet() : []
    if (nullCtr)
        error "ERROR: container paths start with 'null/' (${nullCtr.take(3).join(', ')}...) - params.container_dir was unset " +
              "when nextflow.config was read. Set container_dir in the params file (a later -c config is read too late)."
    def runOpts = workflow.session?.config?.navigate('singularity.runOptions')?.toString() ?: ''
    def unbound = (params.bind_paths instanceof List ? params.bind_paths : []).findAll { bp -> !runOpts.contains("--bind ${java.nio.file.Paths.get(bp.toString()).normalize()}") }
    if (unbound)
        error "ERROR: bind_paths ${unbound} are missing from singularity.runOptions ('${runOpts}'). " +
              "Set bind_paths in the params file (a later -c config is read too late)."

    // -- RNA fusion arm: dedicated validation; skip all DNA-specific checks -
    if (params.assay == 'rna')
        return validateRnaParams(startRank, stopRank, sheet)

    // -- Mid-pipeline entry: require the staged-input directory -
    if (params.start_at == 'annotation') {
        if (!params.input_vcf_dir)
            error "ERROR: params.start_at='annotation' requires params.input_vcf_dir (dir of <sample>.vcf.gz or <sample>_filtered.vcf.gz)."
    } else if (params.start_at == 'classification') {
        if (!params.input_annot_dir)
            error "ERROR: params.start_at='classification' requires params.input_annot_dir (dir of <sample>_vep.vcf.gz)."
    } else if (startRank > stageRank()['mapping']) {
        if (!params.input_bam_dir)
            error "ERROR: params.start_at='${params.start_at}' requires params.input_bam_dir " +
                  "(dir of staged BAM/CRAM inputs, e.g. fastquorum consensus BAMs named <sample_id>_markdup.bam)."
    }
    if (params.dedup_mode == 'external' && startRank <= stageRank()['markdup'])
        error "ERROR: dedup_mode='external' means BAMs are already deduplicated/consensus - " +
              "set start_at to 'metrics' (or later) and point input_bam_dir at them."

    // -- Stage 1 (markdup path): make the inline-UMI trim visible. A umi_read_structure
    //    left set for a library WITHOUT inline UMIs would cut its leading bases.
    if (params.dedup_mode == 'markdup' && params.umi_read_structure && runsStage('mapping', startRank, stopRank))
        log.warn "dedup_mode='markdup' with umi_read_structure='${params.umi_read_structure}'" +
                 (params.umi_read_structure_r2 ? " / R2 '${params.umi_read_structure_r2}'" : "") +
                 ": fastp hard-trims the leading UMI+linker bases of every read. " +
                 "For a library WITHOUT inline UMIs set umi_read_structure: null in the params file."

    // -- Stage 2 (dedup) sanity for UMI modes -------------------
    if (params.dedup_mode in ['umi_consensus', 'umi_markdup']) {
        if (!params.umi_read_structure)
            error "ERROR: dedup_mode='${params.dedup_mode}' requires params.umi_read_structure (e.g. '3M2S+T')."
        if (params.start_at == 'markdup')
            log.warn "dedup_mode='${params.dedup_mode}' with start_at='markdup': input chunk BAMs MUST already " +
                     "carry UMIs in the RX tag (produced by ALIGN_DNA_UMI). Raw fixmate BAMs will fail at consensus."
        if (params.dedup_mode == 'umi_consensus' && params.bqsr && runsStage('bqsr', startRank, stopRank))
            error "ERROR: bqsr=true with dedup_mode='umi_consensus'. Consensus base qualities are model " +
                  "outputs (fgbio), not instrument qualities - recalibrating them is not valid. Set bqsr=false."
    }

    // -- TIER 1 - ESSENTIAL (always) ----------------------------
    requireParams(
        ['pipeline_root', 'work_root', 'outdir', 'container_dir',
         'ref_fasta', 'ref_fai', 'ref_dict'],
        "for any run (core reference + I/O)"
    )

    // -- STAGE 1 mapping --------------------------------------------
    if (runsStage('mapping', startRank, stopRank))
        requireParams(['ref_index_dir', 'fastq_root'], "for the mapping stage (BWA-MEM2 alignment)")

    // -- STAGE 3 metrics --------------------------------------------
    if (runsStage('metrics', startRank, stopRank)) {
        if (params.assay == 'wes')
            requireParams(['interval_bed'], "for WES QC (mosdepth target coverage)")
        if (params.bait_interval_list && !params.target_interval_list)
            log.warn "bait_interval_list set but target_interval_list missing - HsMetrics needs BOTH; it will be skipped."
        if (params.target_interval_list && !params.bait_interval_list)
            log.warn "target_interval_list set but bait_interval_list missing - HsMetrics needs BOTH; it will be skipped."
        if (!params.bait_interval_list && !params.target_interval_list)
            log.warn "bait/target interval_lists unset - CollectHsMetrics will be skipped."
        warnParams(['verifybamid_db'], "contamination QC (VerifyBamID2) will be skipped")
    }

    // -- STAGE 4 bqsr (bqsr: true) ----------------------------------
    if (runsStage('bqsr', startRank, stopRank) && params.bqsr) {
        requireParams(['dbsnp', 'mills', 'known_indels'],
            "when params.bqsr = true (BQSR known-sites: dbSNP + Mills + 1000G indels)")
        if (params.assay == 'wes')
            requireParams(['interval_bed'], "to restrict BQSR to WES capture regions")
    }

    // -- STAGE 5 variant_calling ------------------------------------
    if (runsStage('variant_calling', startRank, stopRank)) {
        requireParams(['scatter_intervals'], "for the variant_calling stage (interval scatter for HC/Mutect2)")

        if (!(params.germline_filter in ['nv', 'none', 'vqsr']))
            error "ERROR: Invalid params.germline_filter '${params.germline_filter}'. One of: nv, none, vqsr"

        if (params.analysis == 'germline') {
            if (!validCallers().contains(params.germline_caller))
                error "ERROR: Invalid params.germline_caller '${params.germline_caller}'. One of: ${validCallers().join(', ')}"
            if (params.germline_caller == 'haplotypecaller') {
                requireParams(['dbsnp'], "for HaplotypeCaller (dbSNP annotation + GenotypeGVCFs)")
                // FilterVariantTranches (nv) and VQSR both pass each of these as a resource;
                // a missing one renders an empty '--resource' and fails hours into the run.
                if (params.germline_filter in ['nv', 'vqsr'])
                    requireParams(['hapmap', 'omni', 'thousandg_snps', 'mills'],
                        "for germline_filter='${params.germline_filter}' (FilterVariantTranches / VQSR resources)")
            }
            if (params.run_gcnv)
                requireParams(['preprocessed_intervals', 'contig_ploidy_priors', 'gcnv_scatter_intervals'],
                    "when params.run_gcnv = true (germline CNV cohort model)")
        }

        if (params.analysis == 'somatic' && params.germline_on_normals) {
            if (!validCallers().contains(params.germline_caller))
                error "ERROR: Invalid params.germline_caller '${params.germline_caller}'. One of: ${validCallers().join(', ')}"
            if (params.germline_caller == 'haplotypecaller')
                requireParams(['dbsnp'] + (params.germline_filter in ['nv', 'vqsr'] ? ['hapmap', 'omni', 'thousandg_snps', 'mills'] : []),
                    "when germline_on_normals = true with HaplotypeCaller (dbSNP + filtering resources)")
            if (params.run_gcnv && params.gcnv_in_somatic)
                requireParams(['preprocessed_intervals', 'contig_ploidy_priors', 'gcnv_scatter_intervals'],
                    "when gcnv_in_somatic = true (germline CNV cohort model on normals)")
        }

        if (params.analysis == 'somatic') {
            requireParams(['germline_resource'], "for somatic calling (Mutect2 gnomAD germline resource)")
            if (params.use_muse_in_consensus)
                requireParams(['dbsnp'], "for MuSE sump (-D dbSNP) when use_muse_in_consensus = true")
            if (!(params.consensus_backbone in ['mutect2', 'union']))
                error "params.consensus_backbone must be 'mutect2' or 'union' (got '${params.consensus_backbone}')"
            warnParams(['panel_of_normals'], "Mutect2 will run without a panel-of-normals (more false positives)")
            warnParams(['cnvkit_pon'], "CNVkit will build its reference from the matched normal (no panel-of-normals)")
            if (params.containsKey('msi_cutoff'))
                error "ERROR: params.msi_cutoff was replaced by msi_cutoff_paired_low / msi_cutoff_paired_high (paired) and msi_cutoff_tumor_only"
            warnParams(['msisensor_list'], "MSI (MSIsensor-pro) will be skipped - needs a pre-scanned microsatellite list")
            warnParams(['pileup_variants'],
                "contamination filtering (GetPileupSummaries/CalculateContamination) will be skipped")
            if (params.orientation_bias)
                log.info "Orientation-bias filtering ENABLED (Mutect2 f1r2 - ob-priors) - recommended for FFPE."
        }

        if (params.analysis == 'tumor_only')
            error "Tumor-only workflow not yet implemented."
    }

    // -- STAGE 6 annotation -----------------------------------------
    if (runsStage('annotation', startRank, stopRank)) {
        requireParams(['vep_cache_dir'], "for the annotation stage (offline VEP cache)")
        // TMB numerator = coding consequences, so the denominator must be coding territory:
        // capture-only (incl. UTR/intronic baits) deflates TMB by the non-coding fraction.
        if (params.analysis == 'somatic' && params.compute_tmb && params.assay == 'wes')
            requireParams(['tmb_coding_bed'],
                "for WES TMB (CDS BED of the genome build; territory = capture x CDS)")
        warnParams(['clinvar_vcf'], "VEP ClinVar custom annotation will be skipped (no clinvar VCF)")
        warnParams(['vep_plugins_dir'], "VEP plugins custom annotation will be skipped (no plugins)")
    }

    return [startRank, stopRank, sheet]
}


// ============================================================
//  VALIDATION - RNA FUSION ARM (called from validateParams when assay=='rna')
//  Returns [startRank, stopRank]. Uses the shared stage vocabulary; the RNA
//  path exercises mapping - metrics - variant_calling - annotation -
//  classification (markdup/bqsr are DNA-only and never gated here).
// ============================================================
def validateRnaParams(int startRank, int stopRank, String sheet) {

    // TIER 1 - always required (core reference + I/O for the fusion arm).
    requireParams(
        ['pipeline_root', 'work_root', 'outdir', 'container_dir'],
        "for the RNA fusion workflow (I/O)"
    )
    // The CTAT library is needed for STAR-Fusion / FusionInspector, for its genome + GTF
    // (rna_use_ctat_reference) and for its STAR index when star_index is not given; an
    // Arriba-only run on its own reference and STAR index does not need it.
    def needs_ctat = params.rna_use_ctat_reference || params.use_starfusion || params.use_fusioninspector || !params.star_index
    if (needs_ctat)
        requireParams(['ctat_lib'], "for the RNA arm (CTAT library: STAR-Fusion/FusionInspector, rna_use_ctat_reference, or its STAR index)")
    if (params.rna_use_ctat_reference && params.star_index)
        log.warn "star_index='${params.star_index}' is set while rna_use_ctat_reference=true: STAR aligns to that index " +
                 "but Arriba/CRAM/counts use the CTAT genome+GTF. Use star_index: null (CTAT index) or rna_use_ctat_reference: false."
    if (params.rna_use_ctat_reference) {
        ['ref_genome.fa', 'ref_genome.fa.fai', 'ref_annot.gtf'].each { f ->
            if (!file("${params.ctat_lib}/${f}").exists())
                error "ERROR: rna_use_ctat_reference=true but ${params.ctat_lib}/${f} is missing" +
                      (f.endsWith('.fai') ? " (run: samtools faidx ${params.ctat_lib}/ref_genome.fa)" : "")
        }
    } else {
        requireParams(['ref_fasta', 'ref_fai', 'gtf'], "for the RNA arm when rna_use_ctat_reference = false")
    }

    // Annotation/classification consume FUSION_CONSENSUS outputs from the SAME run;
    // there is no staged-input path for them, so re-entry there cannot work.
    if (startRank > stageRank()['variant_calling'])
        error "ERROR: assay='rna' supports start_at up to 'variant_calling' (annotation/classification " +
              "need this run's fusion consensus). Got start_at='${params.start_at}'."

    // Mid-pipeline re-entry needs staged STAR outputs.
    if (startRank > stageRank()['mapping'] && !params.input_bam_dir)
        error "ERROR: start_at='${params.start_at}' requires params.input_bam_dir " +
              "(dir of <sample_id>.Aligned.sortedByCoord.out.bam [+ .Chimeric.out.junction])."

    // -- STAGE 1 mapping (RNA: STAR) --------------------------------
    if (runsStage('mapping', startRank, stopRank)) {
        requireParams(['fastq_root'], "for STAR alignment (FASTQ input root)")
        if (!params.star_index && !params.ctat_lib)
            error "ERROR: mapping needs params.star_index, or params.ctat_lib (its STAR index is used in shared mode)."
    }

    // -- STAGE 3 metrics (RNA) --------------------------------------
    // Without REF_FLAT, Picard CollectRnaSeqMetrics is skipped (strandedness is still reported).
    if (runsStage('metrics', startRank, stopRank))
        warnParams(['refflat'], "RNA-seq metrics (Picard CollectRnaSeqMetrics) will be skipped")

    // -- STAGE 5 variant_calling (RNA: fusion calling) --------------
    if (runsStage('variant_calling', startRank, stopRank)) {
        if (!(params.use_arriba || params.use_starfusion || params.use_fusioncatcher))
            error "ERROR: at least one fusion caller must be enabled (use_arriba / use_starfusion / use_fusioncatcher)."
        if (params.use_arriba)
            warnParams(['arriba_blacklist'], "Arriba precision drops sharply without its blacklist")
        if (params.use_starfusion || params.use_fusioninspector)
            requireParams(['ctat_lib'], "for STAR-Fusion / FusionInspector (CTAT genome library)")
        if (params.use_fusioncatcher)
            requireParams(['fusioncatcher_db'], "for FusionCatcher (or set use_fusioncatcher=false)")
        // STAR writes the chimeric junction format STAR-Fusion >= 1.9 reads (header + multimapping
        // columns) only with --chimMultimapNmax > 0, which STAR refuses together with the WithinBAM
        // output Arriba needs; one shared STAR run therefore cannot feed both callers.
        if (params.rna_shared_star == true && params.use_starfusion)
            error "ERROR: rna_shared_star: true cannot feed STAR-Fusion (STAR cannot write STAR-Fusion's chimeric junction format together with Arriba's WithinBAM output). Set rna_shared_star: false (STAR-Fusion then runs its own STAR), or use_starfusion: false."
        warnParams(['known_fusions_tsv'], "single-caller known-fusion RESCUE and known-fusion tiering will be limited")
    }

    // -- STAGE 6 annotation (RNA) -----------------------------------
    if (runsStage('annotation', startRank, stopRank))
        warnParams(['cosmic_cgc_tsv', 'known_fusions_tsv'],
            "fusion tiering will be limited (no COSMIC CGC gene roles / known-fusion catalog)")

    return [startRank, stopRank, sheet]
}


// ============================================================
//  RUN SUMMARY  (function - called from workflow{})
// ============================================================
// pad helper for the summary table (top-level function; strict-safe)
// Resume guard. With the strict (v2) syntax parser, default since Nextflow 26.04, a
// parameter used inside a process script is not part of the task hash: after a parameter
// change, -resume reuses the old task results (verified on 26.04.1: a changed
// consensus_backbone re-used the old consensus). Store the parameters next to the work
// directory and refuse a resumed run whose parameters changed, unless allowed.
// Parameters that cannot change task results are not compared.
def checkResumeParams() {
    // start_at / stop_at only choose which stages run (extending a run with -resume is the
    // normal use); they change no task's command, so they are not treated as a change.
    def skip = ['outdir', 'work_root', 'detailed_reports', 'bind_paths', 'write_manifest',
                'allow_param_change_on_resume', 'start_at', 'stop_at'] as Set
    def keys = params.keySet().findAll { k -> !(k in skip) && !k.startsWith('lsf_') && !k.contains('-') }.sort()
    def cur  = keys.collectEntries { k -> [(k): (params[k] == null ? 'null' : params[k].toString())] }
    def snap = file("${params.work_root}/.papio_params.tsv")
    if (workflow.resume && snap.exists()) {
        def old = snap.text.readLines().findAll { l -> l }.collectEntries { l ->
            def kv = l.split('\t', 2)
            [(kv[0]): (kv.size() > 1 ? kv[1] : '')]
        }
        def changed = (old.keySet() + cur.keySet()).findAll { k -> !(k in skip) && old[k] != cur[k] }.sort()
        if (changed) {
            def lines = changed.collect { k -> "    ${k}: ${old[k] ?: '(unset)'} -> ${cur[k] ?: '(unset)'}" }.join('\n')
            def msg = "Parameters changed since the run being resumed:\n${lines}\n" +
                      "  Nextflow's strict parser does not re-run cached tasks for a parameter change, so the affected\n" +
                      "  results would be stale. Run without -resume, or with a new work_root, or set\n" +
                      "  allow_param_change_on_resume: true if no cached task uses these parameters."
            if (params.allow_param_change_on_resume)
                log.warn msg
            else
                error "ERROR: ${msg}"
        }
    }
    snap.parent.mkdirs()
    snap.text = cur.collect { k, v -> "${k}\t${v}" }.join('\n') + '\n'
}

def pad(s, n) { (s ?: 'not set').toString().padRight(n) }

// ============================================================
def printSummary() {
    def col = 110
    log.info """
----------------------------------------------------------------
-              CLINICAL EXOME PIPELINE                         -
----------------------------------------------------------------
-  Assay           : ${pad(params.assay, col)}-
-  Analysis        : ${pad(params.analysis, col)}-
-  Germline caller : ${pad(params.germline_caller, col)}-
-  BQSR enabled    : ${pad(params.bqsr.toString(), col)}-
-  Dedup mode      : ${pad(params.dedup_mode, col)}-
-  Start stage     : ${pad(params.start_at, col)}-
-  Stop  stage     : ${pad(params.stop_at, col)}-
----------------------------------------------------------------
-  Pipeline root   : ${pad(params.pipeline_root, col)}-
-  Ref FASTA       : ${pad(params.ref_fasta, col)}-
-  Ref index dir   : ${pad(params.ref_index_dir, col)}-
-  Scatter ints    : ${pad(params.scatter_intervals, col)}-
-  Fastq root      : ${pad(params.fastq_root, col)}-
-  Containers      : ${pad(params.container_dir, col)}-
-  Work root       : ${pad(params.work_root, col)}-
-  Outdir          : ${pad(params.outdir, col)}-
----------------------------------------------------------------
-  Interval BED    : ${pad(params.interval_bed, col)}-
-  VEP cache       : ${pad(params.vep_cache_dir, col)}-
-  CNVkit PON      : ${pad(params.cnvkit_pon, col)}-
-  MSIsensor list  : ${pad(params.msisensor_list, col)}-
----------------------------------------------------------------
""".stripIndent()
}


// ============================================================
//  ENTRY WORKFLOW  (strict 26.04: main / onComplete / onError sections)
//  All executable code lives in sections of the entry workflow.
// ============================================================
workflow {
    main:
    def (startRank, stopRank, sheet) = validateParams()
    checkResumeParams()
    log.info "Samplesheet: ${sheet}"
    printSummary()

    if (params.assay in ['wes', 'wgs']) {
        DNA_WORKFLOW(sheet)
    } else if (params.assay == 'rna') {
        RNA_WORKFLOW(sheet)
    } else {
        error "ERROR: Unhandled assay type '${params.assay}'"
    }

    onComplete:
    log.info """
------------------------------------------------------------
-  Pipeline completed: ${(workflow.success ? 'SUCCESS' : 'FAILED').padRight(68)}-
-  Duration         : ${workflow.duration.toString().padRight(68)}-
-  Exit status      : ${workflow.exitStatus.toString().padRight(68)}-
-  Output directory : ${(params.outdir ?: 'not set').toString().padRight(68)}-
------------------------------------------------------------
""".stripIndent()
    // Failure details live here, not in an onError: section - on Nextflow 26.04 the onError
    // handler is invoked with an argument the section does not accept, which raised a second,
    // misleading "Invalid method invocation doCall ... _closure" error on every failed run.
    if (!workflow.success) {
        log.error "Pipeline failed. Error: ${workflow.errorMessage}"
        log.error "Nextflow session: ${workflow.sessionId}"
        log.error "Resume with: nextflow run main.nf -resume ${workflow.sessionId}"
    }
}
