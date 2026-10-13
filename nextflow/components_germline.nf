#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Shared JVM-options helper (single source of truth in utils).
include { gatkJavaOpts; pubDir } from './utils'

// ============================================================
//  components_germline.nf  -  All Germline Process Definitions
// ============================================================


// ============================================================
//  STAGE 5A: GERMLINE SNP/INDEL - HaplotypeCaller (GVCF mode)
// ============================================================
process HAPLOTYPECALLER {
    tag { "${meta.sample_id}:${interval_name}" }
    label 'process_variant'

    maxRetries    2

    input:
    tuple val(meta), path(bam), path(bai), val(interval_name), path(interval_file)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(meta), val(interval_name),
          path("${meta.sample_id}.${interval_name}.g.vcf.gz"),
          path("${meta.sample_id}.${interval_name}.g.vcf.gz.tbi"), emit: gvcf_chunk

    script:
    def java_opts = gatkJavaOpts(task)

    """
    set -euo pipefail
    gatk --java-options "${java_opts}" HaplotypeCaller \\
        -R ${ref_fasta} -I ${bam} \\
        --dbsnp ${dbsnp} \\
        -L ${interval_file}\\
        -ERC GVCF \\
        -O ${meta.sample_id}.${interval_name}.g.vcf.gz
    """
}



// ============================================================
//  STAGE 5A-GATHER: Merge per-interval GVCFs - one per sample
// ============================================================
process MERGE_GVCFS {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    publishDir { "${pubDir('variant_calling', true)}/gvcf" }, mode: 'copy'

    input:
    tuple val(meta), path(gvcf_chunks), path(tbi_chunks)
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(meta),
          path("${meta.sample_id}.g.vcf.gz"),
          path("${meta.sample_id}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def java_opts = gatkJavaOpts(task)
    def inputs    = gvcf_chunks.collect { v -> "-I ${v}" }.join(' \\\n        ')

    """
    set -euo pipefail
    gatk --java-options "${java_opts}" MergeVcfs \\
        ${inputs} -O ${meta.sample_id}.g.vcf.gz
    """
}


// ============================================================
//  STAGE 5A-GENOTYPE: GenotypeGVCFs  (NEW)
//
//  Converts per-sample GVCF - genotyped VCF.
//  Run per-sample for single-sample genotyping.
//  For joint cohort genotyping, run on a merged GenomicsDB
//  (not implemented here - add GenomicsDBImport for large cohorts).
//
//  This step is REQUIRED before CNN scoring or VQSR.
//  Raw GVCF output is not valid variant calls.
// ============================================================
process GENOTYPE_GVCFS {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    maxRetries    2

    input:
    tuple val(meta), path(gvcf), path(gvcf_tbi)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  dbsnp
    path  dbsnp_tbi
    path  interval_bed

    output:
    tuple val(meta),
          path("${meta.sample_id}_genotyped.vcf.gz"),
          path("${meta.sample_id}_genotyped.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    // WES: pad targets so splice-site/edge variants are genotyped (GATK exome practice).
    def pad = (params.wes_interval_padding != null) ? params.wes_interval_padding : 100
    def interval_flag = (meta.assay == 'wes' && interval_bed) ? "-L ${interval_bed} -ip ${pad}" : ''

    """
    set -euo pipefail

    tmp_dir="\$PWD/${meta.sample_id}_ggt_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V ${gvcf} \\
        --dbsnp ${dbsnp} \\
        ${interval_flag} \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}_genotyped.vcf.gz
    """
}


// ============================================================
//  STAGE 5A-FILTER: NVScoreVariants  (single-sample DL filter)
//
//  Single-sample germline variant filtering (params.germline_filter=='nv').
//  Replaces VQSR for single samples or small cohorts (<30).
//  Uses the 2D read-level CNN model to score each variant.
//
//  Two-pass approach:
//    1. NVScoreVariants:      the model is set by nv_tensor_type
//    2. FilterVariantTranches: apply score thresholds for PASS/FAIL
//
//  Use the joint cohort path (params.germline_filter=='vqsr') for large
//  cohorts (see JOINT COHORT GENOTYPING section below).
// ============================================================
process NV_SCORE_VARIANTS {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    scratch false

    maxRetries 2

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  bam      // Original BAM needed for 2D (read-tensor) model
    path  bai
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(meta),
          path("${meta.sample_id}_nv_annotated.vcf.gz"),
          path("${meta.sample_id}_nv_annotated.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    def tensor  = params.nv_tensor_type ?: 'reference'          // 'reference' = 1D model, 'read_tensor' = 2D model
    def bam_arg = (tensor == 'read_tensor') ? "-I ${bam}" : ''   // GATK rejects -I with the 1D model
    """
    set -euo pipefail

    tmp_dir="\$PWD/${meta.sample_id}_nv_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts} -Djava.io.tmpdir=\$tmp_dir" NVScoreVariants \\
        -R ${ref_fasta} \\
        ${bam_arg} \\
        -V ${vcf} \\
        --tensor-type ${tensor} \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}_nv_annotated.vcf.gz

    # NVScoreVariants streams through Python and does not write a tabix index;
    # create it so the declared .tbi output exists.
    gatk --java-options "${java_opts} -Djava.io.tmpdir=\$tmp_dir" IndexFeatureFile \
        -I ${meta.sample_id}_nv_annotated.vcf.gz
    """
}




// ============================================================
//  STAGE 5A-FILTER-2: FilterVariantTranches  (NEW)
//
//  Applies CNN score tranches to produce final PASS/FAIL flags.
//  Requires truth/training resources: HapMap, Omni, 1000G, Mills.
//  SNP sensitivity: 99.95%  (tune per lab needs)
//  INDEL sensitivity: 99.4% (tune per lab needs)
// ============================================================
process FILTER_VARIANT_TRANCHES {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    publishDir { "${pubDir('variant_calling', true)}/snv" }, mode: 'copy'

    maxRetries    2

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  hapmap          // HapMap 3.3 sites VCF
    path  hapmap_tbi
    path  omni            // 1000G Omni VCF
    path  omni_tbi
    path  thousandg_snps  // 1000G high-confidence SNPs
    path  thousandg_tbi
    path  mills           // Mills gold standard indels
    path  mills_tbi

    output:
    tuple val(meta),
          path("${meta.sample_id}_filtered.vcf.gz"),
          path("${meta.sample_id}_filtered.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    def info_key = ((params.nv_tensor_type ?: 'reference') == 'read_tensor') ? 'CNN_2D' : 'CNN_1D'

    """
    set -euo pipefail

    gatk --java-options "${java_opts}" FilterVariantTranches \\
        -V ${vcf} \\
        --resource ${hapmap} \\
        --resource ${omni} \\
        --resource ${thousandg_snps} \\
        --resource ${mills} \\
        --info-key ${info_key} \\
        --snp-tranche 99.95 \\
        --indel-tranche 99.4 \\
        -O ${meta.sample_id}_filtered.vcf.gz
    """
}


// ============================================================
//  STAGE 5A-ALT: GERMLINE SNP/INDEL - DeepVariant
// ============================================================
process DEEPVARIANT {
    tag { "${meta.sample_id}" }
    label 'process_high_memory'

    publishDir { "${pubDir('variant_calling', true)}/snv/deepvariant" }, mode: 'copy'

    maxRetries    2

    input:
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai
    path  interval_bed

    output:
    tuple val(meta),
          path("${meta.sample_id}.vcf.gz"),
          path("${meta.sample_id}.vcf.gz.tbi"),
          path("${meta.sample_id}.g.vcf.gz"),
          path("${meta.sample_id}.g.vcf.gz.tbi"), emit: vcf_gvcf

    script:
    // WES: padded, contig-clamped, merged regions (DeepVariant has no padding option).
    def pad          = (params.wes_interval_padding != null) ? params.wes_interval_padding : 100
    def use_regions  = (meta.assay == 'wes' && interval_bed)
    def regions_flag = use_regions ? "--regions regions.padded.bed" : ''
    def make_regions = !use_regions ? '' : """
    if echo "${interval_bed}" | grep -q '[.]gz\$'; then CAT=zcat; else CAT=cat; fi
    \$CAT ${interval_bed} \\
        | awk -F'\\t' -v p=${pad} 'BEGIN { OFS = "\\t"; while ((getline l < "${ref_fai}") > 0) { split(l, a, "\\t"); len[a[1]] = a[2] } }
            /^(#|track|browser)/ { next }
            (\$1 in len) { s = \$2 - p; if (s < 0) s = 0; e = \$3 + p; if (e > len[\$1]) e = len[\$1]; print \$1, s, e }' \\
        | LC_ALL=C sort -k1,1 -k2,2n \\
        | awk -F'\\t' 'BEGIN { OFS = "\\t" }
            { if (\$1 == c && \$2 <= e) { if (\$3 > e) e = \$3 } else { if (c != "") print c, s, e; c = \$1; s = \$2; e = \$3 } }
            END { if (c != "") print c, s, e }' > regions.padded.bed
    [ -s regions.padded.bed ] || { echo "ERROR: padded regions BED is empty (contig names of BED vs FASTA?)" >&2; exit 1; }
    """
    def model_type   = (meta.assay == 'wes') ? 'WES' : 'WGS'

    """
    set -euo pipefail
    ${make_regions}
    run_deepvariant \\
        --model_type=${model_type} \\
        --ref=${ref_fasta} \\
        --reads=${bam} \\
        ${regions_flag} \\
        --output_vcf=${meta.sample_id}.vcf.gz \\
        --output_gvcf=${meta.sample_id}.g.vcf.gz \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=\${TMPDIR:-/tmp}/${meta.sample_id}_dv_tmp
    """
}


// ============================================================
//  STAGE 5B-SV: GERMLINE SV - Manta germline mode  (NEW)
//
//  Manta germline mode detects: DEL, DUP, INV, INS, BND.
//  WES SV limits apply: only breakpoints within/near capture
//  regions are detectable.
//
//  Separate from somatic Manta (MANTA process below) because:
//    - germline uses configManta.py without --tumorBam
//    - output files have different names (diploidSV vs somaticSV)
//    - no indel candidate output needed (no Strelka2 in germline)
// ============================================================
process MANTA_GERMLINE {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    publishDir { "${pubDir('variant_calling', true)}/sv/manta" }, mode: 'copy'

    maxRetries    2

    input:
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai
    path  call_regions_gz
    path  call_regions_tbi

    output:
    tuple val(meta),
          path("${meta.sample_id}_manta_diploidSV.vcf.gz"),
          path("${meta.sample_id}_manta_diploidSV.vcf.gz.tbi"), emit: sv_vcf
    tuple val(meta),
          path("${meta.sample_id}_manta_candidateSV.vcf.gz"),
          path("${meta.sample_id}_manta_candidateSV.vcf.gz.tbi"), emit: candidate_sv

    script:
    def call_regions = (meta.assay == 'wes' && call_regions_gz) ? "--callRegions ${call_regions_gz}" : ''
    def exome_flag   = (meta.assay == 'wes') ? '--exome' : ''

    """
    set -euo pipefail

    # Manta requires <bam>.bai; our index is <sample>_markdup.bai / _consensus.bai.
    if [ ! -f "${bam}.bai" ]; then ln -s ${bai} "${bam}.bai"; fi

    configManta.py \\
        --bam            ${bam} \\
        --referenceFasta ${ref_fasta} \\
        ${call_regions} \\
        ${exome_flag} \\
        --runDir manta_germline_run

    manta_germline_run/runWorkflow.py -m local -j ${task.cpus}

    cp manta_germline_run/results/variants/diploidSV.vcf.gz       ${meta.sample_id}_manta_diploidSV.vcf.gz
    cp manta_germline_run/results/variants/diploidSV.vcf.gz.tbi   ${meta.sample_id}_manta_diploidSV.vcf.gz.tbi
    cp manta_germline_run/results/variants/candidateSV.vcf.gz     ${meta.sample_id}_manta_candidateSV.vcf.gz
    cp manta_germline_run/results/variants/candidateSV.vcf.gz.tbi ${meta.sample_id}_manta_candidateSV.vcf.gz.tbi
    """
}


// ============================================================
//  STAGE 6A-GERMLINE: GERMLINE CNV - GATK gCNV scatter  (NEW)
//
//  GATK gCNV pipeline (3-step):
//    Step 1: CollectReadCounts - per-sample read depth per interval
//    Step 2: DetermineGermlineContigPloidy - cohort ploidy model
//    Step 3: GermlineCNVCaller - cohort-level CNV model
//    Step 4: PostprocessGermlineCNVCalls - per-sample VCF output
//
//  Requires a cohort of -30 normal samples to build a reliable
//  model. For small cohorts (<30), use ExomeDepth instead.
//
//  This process: Step 1 - CollectReadCounts (per sample)
// ============================================================
process GATK_COLLECT_READ_COUNTS {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    // publishDir { "${pubDir('variant_calling', true)}/cnv/gcnv/counts" }, mode: 'copy'

    maxRetries    2

    input:
    tuple val(meta), path(bam), path(bai)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  preprocessed_intervals   // Output of PreprocessIntervals (run once per cohort)

    output:
    tuple val(meta), path("${meta.sample_id}.counts.hdf5"), emit: counts

    script:
    def java_opts = gatkJavaOpts(task)

    """
    set -euo pipefail

    gatk --java-options "${java_opts}" CollectReadCounts \\
        -I ${bam} \\
        -R ${ref_fasta} \\
        -L ${preprocessed_intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O ${meta.sample_id}.counts.hdf5
    """
}



// ============================================================
//  STAGE 6A-GERMLINE-2: GATK gCNV - DetermineGermlineContigPloidy
//
//  Cohort-level ploidy estimation.
//  Requires ALL sample count HDF5 files simultaneously.
//  Run once per cohort batch.
// ============================================================
process GATK_PLOIDY_MODEL {
    tag   "cohort_ploidy"
    label 'process_ploidy'

    scratch false

    // publishDir { "${pubDir('variant_calling', true)}/cnv/gcnv/ploidy_model" }, mode: 'copy'

    maxRetries    1

    input:
    path  count_files             // List of all sample .counts.hdf5 files
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  preprocessed_intervals
    path  contig_ploidy_priors    // TSV: CONTIG_NAME, PLOIDY_PRIOR_0, PLOIDY_PRIOR_1, ...

    output:
    path "ploidy_model",          emit: ploidy_model
    path "ploidy_calls",          emit: ploidy_calls

    script:
    def inputs    = count_files.collect { v -> "-I ${v}" }.join(' \\\n        ')
    def java_opts = gatkJavaOpts(task)
    
    """
    set -euo pipefail

    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export PYTENSORRC="\$PWD/.pytensorrc"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" "\$TMPDIR" .pytensor
    { echo "[global]"; echo "base_compiledir=\$PWD/.pytensor"; } > "\$PYTENSORRC"

    gatk --java-options "${java_opts} -Djava.io.tmpdir=\$TMPDIR" DetermineGermlineContigPloidy \\
        ${inputs} \\
        -L ${preprocessed_intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --contig-ploidy-priors ${contig_ploidy_priors} \\
        --output-prefix cohort \\
        -O ploidy_model

    mkdir -p ploidy_calls
    cp -r ploidy_model/cohort-calls/* ploidy_calls/
    """
}


// ============================================================
//  STAGE 6A-GERMLINE-3: GATK gCNV - GermlineCNVCaller (scatter)
//
//  Runs per-shard of intervals (scatter).
//  Each shard is an independent cohort model.
//  Gather in PostprocessGermlineCNVCalls.
// ============================================================
process GATK_GCNV_CALLER {
    tag { "shard_${shard_index}" }
    label 'process_medium'

    scratch false

    maxRetries    1

    input:
    val   shard_index
    path  count_files             // All sample counts
    path  ploidy_calls            // From GATK_PLOIDY_MODEL
    path  preprocessed_intervals
    path  interval_shard          // This shard's interval list
    path  annotated_intervals     // GC-content annotations (GC-bias correction); [] if unused


    output:
    tuple val(shard_index),
          path("gcnv_model_shard_${shard_index}"),
          path("gcnv_calls_shard_${shard_index}"), emit: shard_output

    script:
    def java_opts = gatkJavaOpts(task)
    def inputs    = count_files.collect { v -> "-I ${v}" }.join(' \\\n        ')
    // GC-bias correction is applied only when an annotated-intervals TSV is provided.
    def annot_flag = annotated_intervals ? "--annotated-intervals ${annotated_intervals}" : ''

    """
    set -euo pipefail

    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export PYTENSORRC="\$PWD/.pytensorrc"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" "\$TMPDIR" .pytensor
    { echo "[global]"; echo "base_compiledir=\$PWD/.pytensor"; } > "\$PYTENSORRC"

    gatk --java-options "${java_opts} -Djava.io.tmpdir=\$TMPDIR" GermlineCNVCaller \\
        --run-mode COHORT \\
        ${inputs} \\
        --contig-ploidy-calls ${ploidy_calls} \\
        -L ${interval_shard} \\
         ${annot_flag} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output-prefix shard_${shard_index} \\
        --output gcnv_out_${shard_index}

        mv gcnv_out_${shard_index}/shard_${shard_index}-model gcnv_model_shard_${shard_index}
        mv gcnv_out_${shard_index}/shard_${shard_index}-calls gcnv_calls_shard_${shard_index}

    """
}



// ============================================================
//  STAGE 6A-GERMLINE-4: GATK gCNV - PostprocessGermlineCNVCalls
//
//  Per-sample: gather all shard outputs - one VCF per sample.
//  sample_index must match the order samples were passed to
//  GermlineCNVCaller (0-based).
// ============================================================
process GATK_GCNV_POSTPROCESS {
    tag { "${meta.sample_id}" }
    label 'process_variant'

    scratch false

    publishDir { "${pubDir('variant_calling', true)}/cnv/gcnv" }, mode: 'copy'

    maxRetries    2

    input:
    tuple val(meta), val(sample_index)
    path  calls_shards    // List of gcnv_calls_shard_* directories
    path  model_shards    // List of gcnv_model_shard_* directories
    path  ploidy_calls
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(meta),
          path("${meta.sample_id}_gcnv_intervals.vcf.gz"),
          path("${meta.sample_id}_gcnv_intervals.vcf.gz.tbi"),
          path("${meta.sample_id}_gcnv_segments.vcf.gz"),
          path("${meta.sample_id}_gcnv_segments.vcf.gz.tbi"), emit: cnv_vcf
    tuple val(meta),
          path("${meta.sample_id}_gcnv_denoised_cr.tsv"), emit: denoised_cr
          
    script:
    def java_opts = gatkJavaOpts(task)
    def calls_args   = calls_shards.collect { v -> "--calls-shard-path ${v}" }.join(' \\\n        ')
    def model_args   = model_shards.collect { v -> "--model-shard-path ${v}" }.join(' \\\n        ')

    """
    set -euo pipefail

    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export PYTENSORRC="\$PWD/.pytensorrc"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" "\$TMPDIR" .pytensor
    { echo "[global]"; echo "base_compiledir=\$PWD/.pytensor"; } > "\$PYTENSORRC"

    gatk --java-options "${java_opts}" PostprocessGermlineCNVCalls \\
        ${calls_args} \\
        ${model_args} \\
        --contig-ploidy-calls ${ploidy_calls} \\
        --sample-index ${sample_index} \\
        --output-genotyped-intervals  ${meta.sample_id}_gcnv_intervals.vcf.gz \\
        --output-genotyped-segments   ${meta.sample_id}_gcnv_segments.vcf.gz \\
        --output-denoised-copy-ratios ${meta.sample_id}_gcnv_denoised_cr.tsv \\
        --sequence-dictionary ${ref_dict}
    """
}

// ============================================================
//  JOINT COHORT GENOTYPING + VQSR   (params.germline_filter == 'vqsr')
//
//  Completes the GenomicsDB path the single-sample GENOTYPE_GVCFS block
//  flags as "not implemented here". Standard GATK cohort best practice:
//    per-interval GenomicsDBImport across ALL samples
//      -> joint GenotypeGVCFs (per interval)
//      -> gather to one cohort VCF
//      -> VQSR: VariantRecalibrator + ApplyVQSR (SNP then INDEL)
//      -> split back to per-sample VCFs (keeps the downstream per-sample
//         (meta, vcf, tbi) contract unchanged).
//
//  VQSR needs a reasonably large callset to train (roughly >=30 WGS or a
//  large exome cohort); below that, keep params.germline_filter='nv'.
//  All truth resources (HapMap/Omni/1000G/Mills/dbSNP) are already the
//  ones plumbed into GERMLINE_BLOCK for the single-sample filter.
// ============================================================

//  Consolidate every sample's gVCF into a per-interval GenomicsDB workspace.
//  The gVCF list is broadcast (value channel); the interval shard drives N tasks.
process GENOMICSDB_IMPORT {
    tag { "${interval_name}" }
    label 'process_variant'

    scratch false

    maxRetries    2

    input:
    tuple val(interval_name), path(interval_file)
    path  gvcfs              // value: ALL sample gVCFs (broadcast to every shard)
    path  gvcf_tbis
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(interval_name), path(interval_file), path("gendb_${interval_name}"), emit: gendb

    script:
    def java_opts = gatkJavaOpts(task)
    // GenomicsDBImport is native/off-heap heavy - cap JVM heap low, leave the
    // rest of the request for native allocation, or LSF kills the job.

    def v_args    = gvcfs.collect { g -> "-V ${g}" }.join(' \\\n        ')

    """
    set -euo pipefail

    tmp_dir="\$PWD/gdb_${interval_name}"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" GenomicsDBImport \\
        ${v_args} \\
        -L ${interval_file}\\
        --genomicsdb-workspace-path gendb_${interval_name} \\
        --batch-size 50 \\
        --merge-input-intervals \\
        --tmp-dir "\$tmp_dir"
    """
}


// Standard Broad prep; collapses ref blocks, drops low-qual, caps allele bloat
process REBLOCK_GVCF {
    tag { "${meta.sample_id}" }
    label 'process_low'

    input:
    tuple val(meta), path(gvcf), path(tbi)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(meta), path("${meta.sample_id}_reblocked.g.vcf.gz"),
                     path("${meta.sample_id}_reblocked.g.vcf.gz.tbi"), emit: gvcf

    script:
    def java_opts = gatkJavaOpts(task)
    // GQ bands mirror the gnomAD "Biggest Practices" aggressive two-band
    // scheme (pass >=20 / fail <20). -do-qual-approx adds the INFO
    // annotation GnarlyGenotyper needs to approximate QUAL without
    // iterating genotypes. --floor-blocks keeps block boundaries clean.
    def gq_bands = params.reblock_gq_bands ?: '20'
    """
    set -euo pipefail

    tmp_dir="\$PWD/reblock_${meta.sample_id}"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" ReblockGVCF \\
        -R ${ref_fasta} \\
        -V ${gvcf} \\
        -do-qual-approx \\
        --floor-blocks \\
        -GQB ${gq_bands} \\
        --tmp-dir "\$tmp_dir" \\
        -O ${meta.sample_id}_reblocked.g.vcf.gz
    """
}

//  Joint-genotype the cohort GenomicsDB for one interval shard.
process GENOTYPE_GVCFS_JOINT {
    tag { "${interval_name}" }
    label 'process_long'

    maxRetries    2

    input:
    tuple val(interval_name), path(interval_file), path(gendb)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(interval_name),
          path("joint.${interval_name}.vcf.gz"),
          path("joint.${interval_name}.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task) 

    """
    set -euo pipefail

    tmp_dir="\$PWD/ggvcf_${interval_name}"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts} -DGENOMICSDB_PRINT_TIMERS=0" GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V gendb://${gendb} \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --max-alternate-alleles 3 \\
        --genomicsdb-max-alternate-alleles 4 \\
        --dbsnp ${dbsnp} \\
        -L ${interval_file} \\
        --tmp-dir "\$tmp_dir" \\
        -O joint.${interval_name}.vcf.gz
    """
}


// ============================================================
//  PATH A - GnarlyGenotyper (scalable joint genotyping)
//
//  Drop-in replacement for GENOTYPE_GVCFS_JOINT in the large-cohort path.
//  Consumes a GenomicsDB built from REBLOCKED gVCFs. Runs per interval
//  shard (scatter), identical scatter contract to GENOTYPE_GVCFS_JOINT so
//  the downstream merge/VQSR wiring is unchanged.
//
//  --only-output-calls-starting-in-intervals avoids duplicate emission at
//  shard boundaries (records whose start is outside the shard interval are
//  suppressed), which matters when gathering shards.
// ============================================================
process GNARLY_GENOTYPER {
    tag { "${interval_name}" }
    label 'process_variant'

    scratch false

    maxRetries 2

    input:
    tuple val(interval_name), path(interval_file), path(gendb)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(interval_name),
          path("gnarly.${interval_name}.vcf.gz"),
          path("gnarly.${interval_name}.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    """
    set -euo pipefail

    tmp_dir="\$PWD/gnarly_${interval_name}"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" GnarlyGenotyper \\
        -R ${ref_fasta} \\
        -V gendb://${gendb} \\
        --dbsnp ${dbsnp} \\
        -L ${interval_file} \\
        --only-output-calls-starting-in-intervals \\
        --tmp-dir "\$tmp_dir" \\
        -O gnarly.${interval_name}.vcf.gz
    """
}


// ============================================================
//  PATH B - GLnexus (DeepVariant cohort joint genotyping)  [STUB]
//
//  NOT WIRED INTO THE WORKFLOW YET. Documented scaffold only.
//
//  Rationale: the DeepVariant caller branch in GERMLINE_BLOCK currently
//  performs SINGLE-SAMPLE calling only - it emits per-sample VCFs and its
//  gVCFs are not merged. GLnexus is the scalable joint genotyper designed
//  for DeepVariant gVCFs (validated at 50k-240k exomes). Wiring this in
//  would give the DeepVariant path a true multi-sample cohort callset.
//
//  Before enabling:
//    - Provide a GLnexus container (params.container_glnexus) - this is a
//      NEW tool, not in the GATK image. Verify it can be obtained/run in
//      the target environment (egress-restricted clusters may block it).
//    - GLnexus is single-machine multithreaded (not scatter-gather); size
//      one large node (high mem + CPU) rather than many small shards.
//    - Use the tuned preset: --config DeepVariantWES (or DeepVariantWGS).
//    - Output is BCF; convert to VCF.gz with bcftools if downstream needs it.
//
//  This stub is intentionally minimal and untested. Do not enable without
//  validating GLnexus on the target cluster first.
// ============================================================


//  Gather the per-interval joint VCFs into a single cohort callset.
process MERGE_COHORT_VCF {
    tag "cohort"
    label 'process_medium'

    publishDir { "${pubDir('variant_calling', true)}/joint" }, mode: 'copy'

    input:
    path vcfs
    path tbis

    output:
    tuple path("cohort_joint.vcf.gz"),
          path("cohort_joint.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    def inputs    = vcfs.collect { v -> "-I ${v}" }.join(' \\\n        ')

    """
    set -euo pipefail
    gatk --java-options "${java_opts}" MergeVcfs \\
        ${inputs} \\
        -O cohort_joint.vcf.gz
    """
}


//  Build a VQSR model for one variant class (mode is 'SNP' or 'INDEL').
//  Resources/annotations are selected by mode; DP is dropped for exomes
//  (uneven target coverage), MQ is used for SNPs only.
process VARIANT_RECALIBRATOR {
    tag { "${mode}" }
    label 'process_variant'

    scratch false

    maxRetries    2

    input:
    tuple val(mode), path(vcf), path(tbi)
    path  ref_fasta
    path  ref_fai
    path  ref_dict
    path  hapmap
    path  hapmap_tbi
    path  omni
    path  omni_tbi
    path  thousandg_snps
    path  thousandg_tbi
    path  mills
    path  mills_tbi
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(mode),
          path("cohort.${mode}.recal"),
          path("cohort.${mode}.recal.idx"),
          path("cohort.${mode}.tranches"), emit: model

    script:
    def java_opts = gatkJavaOpts(task)
    def is_wes    = (params.assay == 'wes')
    def dp        = is_wes ? '' : ' -an DP'
    def resources
    def annotations
    if (mode == 'SNP') {
        resources = "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} " +
                    "--resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} " +
                    "--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousandg_snps} " +
                    "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp}"
        annotations = "-an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum${dp}"
    } else {
        resources = "--resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} " +
                    "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp}"
        annotations = "-an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum${dp}"
    }
    def max_gaussians = (mode == 'INDEL') ? '--max-gaussians 4' : '--max-gaussians 6'

    """
    set -euo pipefail

    tmp_dir="\$PWD/vqsr_${mode}"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" VariantRecalibrator \\
        -R ${ref_fasta} \\
        -V ${vcf} \\
        -mode ${mode} \\
        ${resources} \\
        ${annotations} \\
        ${max_gaussians} \\
        -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        --tmp-dir "\$tmp_dir" \\
        -O cohort.${mode}.recal \\
        --tranches-file cohort.${mode}.tranches
    """
}


//  Apply a VQSR model. Run SNP first, then INDEL on the SNP-filtered VCF.
process APPLY_VQSR {
    tag { "${mode}" }
    label 'process_medium'

    maxRetries    2

    input:
    tuple val(mode), path(vcf), path(tbi), path(recal), path(recal_idx), path(tranches)
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(mode),
          path("cohort.${mode}_recalibrated.vcf.gz"),
          path("cohort.${mode}_recalibrated.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)
    def ts        = (mode == 'SNP') ? '99.7' : '99.0'

    """
    set -euo pipefail
    gatk --java-options "${java_opts}" ApplyVQSR \\
        -R ${ref_fasta} \\
        -V ${vcf} \\
        -mode ${mode} \\
        --recal-file ${recal} \\
        --tranches-file ${tranches} \\
        --truth-sensitivity-filter-level ${ts} \\
        -O cohort.${mode}_recalibrated.vcf.gz
    """
}


//  Split the VQSR-filtered cohort callset back to per-sample VCFs so the
//  downstream annotation/report flow keeps its (meta, vcf, tbi) contract.
//  Genotypes and VQSR FILTER flags are preserved; hom-ref-only sites dropped.
process SPLIT_COHORT_VCF {
    tag { "${meta.sample_id}" }
    label 'process_medium'

    publishDir { "${pubDir('variant_calling', true)}/snv" }, mode: 'copy'

    input:
    tuple val(meta), path(cohort_vcf), path(cohort_tbi)
    path  ref_fasta
    path  ref_fai
    path  ref_dict

    output:
    tuple val(meta),
          path("${meta.sample_id}_filtered.vcf.gz"),
          path("${meta.sample_id}_filtered.vcf.gz.tbi"), emit: vcf

    script:
    def java_opts = gatkJavaOpts(task)

    """
    set -euo pipefail
    gatk --java-options "${java_opts}" SelectVariants \\
        -R ${ref_fasta} \\
        -V ${cohort_vcf} \\
        --sample-name ${meta.sample_id} \\
        --exclude-non-variants \\
        --remove-unused-alternates \\
        -O ${meta.sample_id}_filtered.vcf.gz
    """
}

// ============================================================
//  PHASE 2: ROH - bcftools roh  (runs of homozygosity)
//
//  Detects long homozygous stretches from the genotyped SNV VCF.
//  Standard for recessive-disease and consanguinity assessment.
//  Uses bcftools (already in the annotation image) - no new container.
//
//  --AF-file (population allele frequencies) sharply improves ROH
//  accuracy. It is OPTIONAL here: when params.roh_af_file is unset the
//  empty-channel default ([]) is passed and we fall back to --estimate-AF,
//  which derives AFs from the callset itself (fine for cohorts, weaker
//  for single samples). We emit the region-level (RG) records only.
// ============================================================
//     // Prefer a provided population-AF file; otherwise self-estimate from the VCF.
//     def af_flag = af_file ? "--AF-file ${af_file}"

//     """
//     set -euo pipefail

//     bcftools roh \\
//         ${af_flag} \\
//         --output-type rz \\
//         --output ${meta.sample_id}_roh.txt.gz \\
//         ${vcf}

//     # Keep only the region-level (RG) records + header for the report.
//     # -G was not used, so bcftools emits both ST and RG lines; RG is the
//     # per-segment summary we want. Decompress to a plain TSV.
//     zcat ${meta.sample_id}_roh.txt.gz \\
//         | awk 'NR==1 || /^# *RG/ || /^RG\\t/' \\
//         > ${meta.sample_id}_roh.txt
//     """
// }


process BCFTOOLS_ROH {
    tag { "${meta.sample_id}" }
    label 'process_low'

    publishDir { "${pubDir('variant_calling', true)}/roh" }, mode: 'copy'

    maxRetries 2

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  af_file
    path  af_file_tbi

    output:
    tuple val(meta), path("${meta.sample_id}_roh.txt"), emit: roh

    script:
    // A population-AF file is mandatory: single-sample self-estimated AF produces
    // unreliable ROH. Fail loudly rather than emit misleading segments.
    if (!af_file) {
        error "BCFTOOLS_ROH requires a population-AF file (params.roh_af_file). " +
              "Self-estimated single-sample AF is not reliable for ROH."
    }
    """
    set -euo pipefail

    bcftools roh \\
        --AF-file ${af_file} \\
        --output-type rz \\
        --output ${meta.sample_id}_roh.txt.gz \\
        ${vcf}

    # Keep the RG (region) header + per-segment RG data rows for the report.
    zcat ${meta.sample_id}_roh.txt.gz \\
        | awk '/^# RG/ || /^RG\\t/' \\
        > ${meta.sample_id}_roh.txt
    """
}

// ============================================================
//  PHASE 1: CONSOLIDATED GERMLINE REPORT
//
//  One sectioned TSV per sample, mirroring CONSOLIDATE_PATIENT_REPORT
//  (somatic) but for the germline variant classes:
//    [1] SNV/indel   (VEP+ClinVar annotated VCF)
//    [2] SV          (AnnotSV, Manta diploidSV)
//    [3] CNV         (AnnotSV, gCNV segments)
//    [4] ROH         (bcftools roh region records)
//
//  Keyed on sample_id (germline is per-sample; no tumor/normal pairing).
//  Every section is defensive: missing/empty inputs yield a "(no ...)"
//  placeholder instead of failing, so partial runs still report.
// ============================================================
process CONSOLIDATE_GERMLINE_REPORT {
    tag { "${sample_id}" }
    label 'process_low'

    publishDir { "${pubDir('classification', true)}/reports" }, mode: 'copy'

    input:
    tuple val(sample_id),
          path(snv_vcf,   stageAs: 'snv.vcf.gz'),
          path(snv_tbi,   stageAs: 'snv.vcf.gz.tbi'),
          path(sv_tsv,    stageAs: 'sv_annotsv.tsv'),
          path(cnv_tsv,   stageAs: 'cnv_annotsv.tsv'),
          path(roh_txt,   stageAs: 'roh.txt')
    path  cgc_tsv,   stageAs: 'cgc/*'        // params.cosmic_cgc_tsv (GERMLINE column -> reportable genes) or []
    path  gene_list, stageAs: 'gene_list/*'  // params.germline_report_genes (one symbol per line) or []

    output:
    tuple val(sample_id), path("${sample_id}_germline_consolidated.tsv"), emit: report

    script:
    def clin_min   = params.germline_report_clin_min ?: 'pathogenic'
    def lof_max_af = (params.germline_lof_max_af != null) ? params.germline_lof_max_af : 0.01
    """
    set -uo pipefail
    REPORT=${sample_id}_germline_consolidated.tsv

    {
      echo "## ============================================================"
      echo "## CONSOLIDATED GERMLINE REPORT - SAMPLE: ${sample_id}"
      echo "## Generated: \$(date '+%Y-%m-%d %H:%M:%S')"
      echo "## NOTE: computational summary for review; not a final ACMG clinical classification."
      echo "## ============================================================"
    } > \$REPORT

    # -- [1] SNV/indel - pull ClinVar-flagged calls from the VEP VCF ----------
    printf '\\n## [1] GERMLINE SNV/INDEL  (VEP+ClinVar; CLIN_SIG matching "%s")\\n' "${clin_min}" >> \$REPORT
    if [ -s snv.vcf.gz ]; then
      printf 'CHROM\\tPOS\\tREF\\tALT\\tGENE\\tCONSEQUENCE\\tIMPACT\\tMAX_AF\\tClinVar\\tGT\\n' >> \$REPORT
      # split-vep flattens the per-transcript CSQ; keep the worst consequence.
      # Filter to rows whose ClinVar significance contains the requested term
      # (case-insensitive); if split-vep lacks the field, emit a note.
      csq_fields=\$(bcftools +split-vep -l snv.vcf.gz 2>/dev/null)
      if grep -qiw CLIN_SIG <<< "\$csq_fields"; then
        # PASS-only; 'conflicting...pathogenicity' contains the substring
        # 'pathogenic' and must not be reported as pathogenic.
        bcftools view -f PASS,. -Ou snv.vcf.gz 2>/dev/null \\
        | bcftools +split-vep -s worst -d \\
          -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%IMPACT\\t%MAX_AF\\t%CLIN_SIG\\t[%GT]\\n' \\
          - 2>/dev/null \\
          | awk -F'\\t' -v term="${clin_min}" 'tolower(\$9) ~ tolower(term) && tolower(\$9) !~ /conflicting/' \\
          | sort -u >> \$REPORT \\
          || echo "(SNV query failed)" >> \$REPORT
      else
        echo "(VEP VCF has no CLIN_SIG field; run ANNOTATE_VCF with ClinVar)" >> \$REPORT
      fi
    else echo "(no annotated SNV/indel VCF)" >> \$REPORT; fi

    # -- [1b] LoF in reportable genes (tumour-agnostic, parameter-driven gene set) -----
    #   gene set = CGC genes flagged GERMLINE (params.cosmic_cgc_tsv) + params.germline_report_genes.
    #   One consequence per gene (VEP PICK) when present. Not LOFTEE-filtered (last-exon /
    #   NMD-escape not assessed) - every call needs review.
    : > report_genes.txt
    if [ -n "${cgc_tsv}" ]; then
        if echo "${cgc_tsv}" | grep -q '[.]gz\$'; then CATC=zcat; else CATC=cat; fi
        "\$CATC" "${cgc_tsv}" | awk -F'\\t' '
            NR == 1 { for (i = 1; i <= NF; i++) { h = toupper(\$i); gsub(/[ ]+/, "_", h)
                          if (h == "GENE_SYMBOL") g = i; if (h == "GERMLINE") m = i }
                      next }
            g && m && tolower(\$m) ~ /^y/ && \$g != "" { print toupper(\$g) }' >> report_genes.txt || true
    fi
    if [ -n "${gene_list}" ]; then
        grep -v -E '^[[:space:]]*(#|\$)' "${gene_list}" | awk '{ print toupper(\$1) }' >> report_genes.txt || true
    fi
    sort -u -o report_genes.txt report_genes.txt
    n_genes=\$(grep -c . report_genes.txt || true)
    printf '\\n## [1b] GERMLINE LoF (HIGH impact) IN REPORTABLE GENES  (%s genes; MAX_AF < %s; not LOFTEE-filtered)\\n' "\$n_genes" "${lof_max_af}" >> \$REPORT
    if [ "\$n_genes" -eq 0 ]; then
        echo "(no reportable gene set: set params.germline_report_genes and/or a CGC file with a GERMLINE column)" >> \$REPORT
    elif [ -s snv.vcf.gz ]; then
        printf 'CHROM\\tPOS\\tREF\\tALT\\tGENE\\tCONSEQUENCE\\tIMPACT\\tMAX_AF\\tClinVar\\tGT\\n' >> \$REPORT
        fields=\$(bcftools +split-vep -l snv.vcf.gz 2>/dev/null | awk '{ print \$2 }')
        clin=\$(echo "\$fields" | grep -qx CLIN_SIG && echo '%CLIN_SIG' || echo '.')
        maf=\$(echo "\$fields"  | grep -qx MAX_AF   && echo '%MAX_AF'   || echo '.')
        fmt="%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%IMPACT\\t\${maf}\\t\${clin}\\t[%GT]\\n"
        lof_filter() { awk -F'\\t' -v maxaf="${lof_max_af}" '
            BEGIN { while ((getline g < "report_genes.txt") > 0) G[g] = 1 }
            \$7 == "HIGH" && (toupper(\$5) in G) && (\$8 == "" || \$8 == "." || \$8 + 0 < maxaf + 0)'; }
        if echo "\$fields" | grep -qx PICK; then
            bcftools view -f PASS,. -Ou snv.vcf.gz 2>/dev/null \\
            | bcftools +split-vep -d -f "%PICK\\t\${fmt}" - 2>/dev/null \\
            | awk -F'\\t' '\$1 == "1"' | cut -f2- | lof_filter | sort -u >> \$REPORT \\
            || echo "(LoF query failed)" >> \$REPORT
        else
            bcftools view -f PASS,. -Ou snv.vcf.gz 2>/dev/null \\
            | bcftools +split-vep -s worst -d -f "\${fmt}" - 2>/dev/null \\
            | lof_filter | sort -u >> \$REPORT \\
            || echo "(LoF query failed)" >> \$REPORT
        fi
    else echo "(no annotated SNV/indel VCF)" >> \$REPORT; fi

    # AnnotSV rows worth reading: one row per SV (Annotation_mode=full), kept when
    # ACMG class 4/5 OR overlapping a reportable gene. Unknown layout -> all rows.
    filter_annotsv() {
        awk -F'\\t' 'BEGIN { while ((getline g < "report_genes.txt") > 0) G[g] = 1 }
            NR == 1 { for (i = 1; i <= NF; i++) h[\$i] = i; print; next }
            { if (("Annotation_mode" in h) && \$(h["Annotation_mode"]) != "full") next
              keep = !(("ACMG_class" in h) || ("Gene_name" in h))
              if (("ACMG_class" in h) && \$(h["ACMG_class"]) ~ /(^|=)(4|5)\$/) keep = 1
              if (!keep && ("Gene_name" in h)) {
                  n = split(toupper(\$(h["Gene_name"])), gs, /[;,]/)
                  for (j = 1; j <= n; j++) if (gs[j] in G) { keep = 1; break } }
              if (keep) print }' "\$1"
    }

    # -- [2] SV - AnnotSV (Manta diploidSV), ACMG class 4/5 highlighted ------
    printf '\\n## [2] GERMLINE SV  (AnnotSV, Manta diploidSV; full rows, ACMG 4/5 or reportable gene)\\n' >> \$REPORT
    if [ -s sv_annotsv.tsv ] && [ \$(wc -l < sv_annotsv.tsv) -gt 1 ]; then
      filter_annotsv sv_annotsv.tsv >> \$REPORT
    else echo "(no annotated SV output)" >> \$REPORT; fi

    # -- [3] CNV - AnnotSV (gCNV segments) -----------------------------------
    printf '\\n## [3] GERMLINE CNV  (AnnotSV, GATK gCNV segments; full rows, ACMG 4/5 or reportable gene)\\n' >> \$REPORT
    if [ -s cnv_annotsv.tsv ] && [ \$(wc -l < cnv_annotsv.tsv) -gt 1 ]; then
      filter_annotsv cnv_annotsv.tsv >> \$REPORT
    else echo "(no annotated CNV output)" >> \$REPORT; fi

    # -- [4] ROH - bcftools roh region records -------------------------------
    printf '\\n## [4] RUNS OF HOMOZYGOSITY  (bcftools roh, region-level)\\n' >> \$REPORT
    if [ -s roh.txt ]; then cat roh.txt >> \$REPORT; else echo "(no ROH output)" >> \$REPORT; fi
    """
}
