#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Shared JVM-options helper (single source of truth in utils), same import
// idiom as components_common / _somatic / _germline.
include { gatkJavaOpts; pubDir } from './utils'

// ============================================================
//  components_rna.nf  -  RNA Fusion Process Definitions
//
//  STAGE MAP (reuses the shared --start_at/--stop_at vocabulary):
//    mapping         - STAR_ALIGN         (one shared 2-pass STAR run)
//    metrics         - MARKDUP_RNA (QC-only) + RNASEQ_METRICS + samtools stats
//    variant_calling - ARRIBA + STARFUSION + FUSIONCATCHER
//                      - FUSION_CONSENSUS - FUSIONINSPECTOR
//    annotation      - ANNOTATE_FUSIONS
//    classification  - CONSOLIDATE_FUSION_REPORT
//
//  KEY EFFICIENCY DECISION - "run STAR once":
//    A single STAR alignment emits BOTH (a) a coordinate-sorted BAM that
//    Arriba reads directly and (b) a STAR-Fusion-compatible Chimeric.out
//    .junction for STAR-Fusion "kickstart" (no second STAR). FusionCatcher
//    runs its own aligners on the raw FASTQs (independent by design - its
//    sensitivity comes from a different alignment strategy + curated DB).
//    This replaces the naive 3-STAR layout with 1-STAR + FusionCatcher,
//    the single largest CPU/RAM/disk saving in the arm. A fidelity escape
//    hatch (params.rna_shared_star=false) reverts STAR-Fusion to its own
//    bundled STAR run.
//
//    IMPORTANT invariant for the shared path: the STAR index MUST be built
//    from the SAME genome+GTF as the CTAT genome library (params.ctat_lib)
//    and use the SAME chromosome naming, or STAR-Fusion kickstart silently
//    yields zero fusions. See STAR_GENOME_INDEX + README.
//
//  DEDUP NOTE: RNA duplicate marking here is QC-ONLY. Fusion callers are
//  fed the RAW STAR BAM / FASTQs, never the dup-marked BAM - deduplication
//  removes true fusion-supporting reads at low input and every caller
//  models duplicates internally.
// ============================================================


// ============================================================
//  STAGE 1: STAR ALIGNMENT  (one shared run - Arriba BAM + STAR-Fusion junction)
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
    path  star_index
    path  gtf

    output:
    tuple val(meta),
          path("${meta.sample_id}.Aligned.sortedByCoord.out.bam"),
          path("${meta.sample_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sample_id}.Chimeric.out.junction"),  emit: chimeric
    path  "${meta.sample_id}.Log.final.out",                           emit: log
    path  "${meta.sample_id}.SJ.out.tab",                              emit: sj
    // Gene counts (cols 2/3/4 = unstranded / fwd / rev strand) from the GTF STAR loaded.
    tuple val(meta), path("${meta.sample_id}.ReadsPerGene.out.tab"),  emit: gene_counts



    script:
    // Reserve 2 threads for the concurrent samtools sort/index; STAR gets the rest.
    def star_threads = Math.max(task.cpus - 2, 4)
    def r1_list      = (reads1 instanceof List ? reads1 : [reads1]).join(',')
    def r2_list      = (reads2 && (reads2 instanceof List ? reads2.size() : 1))
                       ? (reads2 instanceof List ? reads2 : [reads2]).join(',') : ''
    def read_files   = r2_list ? "${r1_list} ${r2_list}" : "${r1_list}"
    // Auto-handle .gz vs plain FASTQ.
    def rf_cmd       = "--readFilesCommand zcat"
    def twopass      = (params.star_twopass == false) ? 'None' : 'Basic'
    def rg_line      = "ID:${meta.sample_id} SM:${meta.sample_id} LB:${meta.batch} PL:ILLUMINA"

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
        --sjdbGTFfile ${gtf} \\
        --twopassMode ${twopass} \\
        --quantMode GeneCounts \\
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 \\
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

// --outSAMtype BAM SortedByCoordinate \\
//         --outSAMunmapped Within --outBAMcompression 6 \\
//         --limitBAMsortRAM ${(task.memory.toBytes() * 0.75) as long} \\


// ============================================================
//  STAGE 3A: RNA duplicate marking  (QC-ONLY - NOT fed to fusion callers)
// ============================================================
process MARKDUP_RNA {
    tag { "${meta.sample_id}" }
    label 'process_markdup'

    maxRetries 2

    // -- 03 metrics / QC ------------------------------------------------
    publishDir { "${pubDir('metrics', true)}/rnaseq_qc" }, mode: 'copy', pattern: '*_dup_metrics.txt'
    
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta),
          path("${meta.sample_id}_rna_markdup.bam"),
          path("${meta.sample_id}_rna_markdup.bam.bai"), emit: bam
    path  "${meta.sample_id}_dup_metrics.txt",           emit: metrics

    script:
    def java_opts = gatkJavaOpts(task)
    """
    set -euo pipefail

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_rna_md_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    gatk --java-options "${java_opts}" MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.sample_id}_rna_markdup.bam \\
        -M ${meta.sample_id}_dup_metrics.txt \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR "\$tmp_dir"

    samtools index -@ ${task.cpus} ${meta.sample_id}_rna_markdup.bam
    """
}


// ============================================================
//  STAGE 3B: RNA-seq metrics - Picard CollectRnaSeqMetrics
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
    tuple val(meta), path("${meta.sample_id}.rna_metrics"),      emit: metrics
    tuple val(meta), path("${meta.sample_id}.strandedness.txt"), emit: strandedness

    script:
    def java_opts = gatkJavaOpts(task)
    // STRAND per library prep. 'auto' = infer from STAR gene counts (col3 = read-1 sense,
    // col4 = read-1 antisense); without counts (re-entry) 'auto' falls back to NONE.
    def strand = (params.rna_strandedness ?: 'auto').toString().toLowerCase()
    def strand_spec = strand in ['reverse','rf','fr-firststrand'] ? 'SECOND_READ_TRANSCRIPTION_STRAND'
                    : strand in ['forward','fr','fr-secondstrand'] ? 'FIRST_READ_TRANSCRIPTION_STRAND'
                    : 'NONE'
    def infer = (strand == 'auto' && gene_counts) ? 'yes' : 'no'
    def rrna_flag = rrna_interval_list ? "--RIBOSOMAL_INTERVALS ${rrna_interval_list}" : ''

    """

    set -euo pipefail

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_rnaseqc_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir"' EXIT

    # Prefilter to properly-paired reads: Picard's CollectRnaSeqMetrics crashes on
    # interchromosomal discordant pairs. QC-only copy; callers use the raw BAM.

    samtools view -@ ${task.cpus} -f 2 -b ${bam} -o qc_input.bam
    samtools index -@ ${task.cpus} qc_input.bam

    # Library strandedness: inferred from STAR ReadsPerGene when rna_strandedness='auto'
    # (>= 80% of gene-assigned reads on one strand), else the configured value.
    strand_spec="${strand_spec}"; basis="rna_strandedness=${strand}"
    if [ "${infer}" = "yes" ]; then
        read -r fwd rev < <(awk -F'\\t' '!/^N_/ { f += \$3; r += \$4 } END { print f + 0, r + 0 }' ${gene_counts})
        tot=\$((fwd + rev))
        if   [ "\$tot" -gt 0 ] && [ \$((100 * fwd / tot)) -ge 80 ]; then strand_spec="FIRST_READ_TRANSCRIPTION_STRAND"
        elif [ "\$tot" -gt 0 ] && [ \$((100 * rev / tot)) -ge 80 ]; then strand_spec="SECOND_READ_TRANSCRIPTION_STRAND"
        else strand_spec="NONE"; fi
        basis="inferred from STAR gene counts: sense=\$fwd antisense=\$rev"
    fi
    printf 'sample_id\\t%s\\nstrand_specificity\\t%s\\nbasis\\t%s\\n' "${meta.sample_id}" "\$strand_spec" "\$basis" > ${meta.sample_id}.strandedness.txt

    gatk --java-options "${java_opts}" CollectRnaSeqMetrics \\
        -I qc_input.bam \\
        -O ${meta.sample_id}.rna_metrics \\
        --REF_FLAT ${refflat} \\
        --STRAND_SPECIFICITY \$strand_spec \\
        ${rrna_flag} \\
        --VALIDATION_STRINGENCY SILENT --TMP_DIR "\$tmp_dir"
    """
}



// ============================================================
//  STAGE 5A: FUSION CALLER 1 - Arriba
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
    def strand   = (params.rna_strandedness ?: 'auto').toString().toLowerCase()
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
//  STAGE 5B: FUSION CALLER 2 - STAR-Fusion
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

    tmp_dir="\${TMPDIR:-/tmp}/${meta.sample_id}_sf_tmp"
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
//  STAGE 5C: FUSION CALLER 3 - FusionCatcher  (independent; heavy, optional)
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

    """
    set -euo pipefail

    tmp_dir="\$PWD/fc_tmp"
    mkdir -p "\$tmp_dir"
    trap 'rm -rf "\$tmp_dir" fc_input fc_out' EXIT

    mkdir -p fc_input
    cat ${r1_files} > fc_input/${meta.sample_id}_R1.fastq.gz
    ${ reads2 ? "cat ${r2_files} > fc_input/${meta.sample_id}_R2.fastq.gz" : "" }

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
//  STAGE 5D: Multi-caller fusion consensus  (analog of ENSEMBLE_CONSENSUS)
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
    path  known_fusions_tsv

    output:
    tuple val(meta), path("${meta.sample_id}_fusion_consensus.tsv"),  emit: consensus
    tuple val(meta), path("${meta.sample_id}_fusion_candidates.txt"), emit: candidates
    tuple val(meta), path("${meta.sample_id}_fusion_consensus_stats.txt"), emit: stats

    script:
    def min_callers  = params.fusion_min_callers      ?: 2
    def readthru_kb  = params.fusion_readthrough_kb   ?: 200
    def known_arg    = known_fusions_tsv ? "known_fusions.tsv" : 'NONE'
    """
    set -euo pipefail
    ${ known_fusions_tsv ? "cp ${known_fusions_tsv} known_fusions.tsv" : "" }

    cat > fusion_consensus.py <<'PYEOF'
#!/usr/bin/env python3
import sys, csv, re

MIN_CALLERS = int(sys.argv[1])
READTHRU_BP = int(sys.argv[2]) * 1000
KNOWN_PATH  = sys.argv[3]
SAMPLE      = sys.argv[4]

def norm_chrom(c):
    if c is None: return ""
    c = str(c).strip()
    c = re.sub(r'^chr', '', c, flags=re.I)
    if c in ('M', 'MT'): return 'MT'
    return c

def norm_gene(g):
    if g is None: return ""
    g = str(g).strip().upper()
    # drop STAR-Fusion ENSG suffix "SYMBOL^ENSG000..." and parenthetical notes
    g = g.split('^')[0]
    g = re.sub(r'\\(.*?\\)', '', g)
    return g.strip()

def canon_key(g5, g3):
    a, b = norm_gene(g5), norm_gene(g3)
    return "--".join(sorted([a, b]))       # order-insensitive for matching

def parse_bp(s):
    # "chr7:55019278:+" or "7:55019278" -> (chrom, pos)
    if not s: return ("", "")
    parts = str(s).split(':')
    if len(parts) >= 2:
        return (norm_chrom(parts[0]), parts[1])
    return ("", "")

def add(store, key, g5, g3, bp5, bp3, caller, support, extra):
    d = store.setdefault(key, {
        "genes": (g5, g3), "bp5": bp5, "bp3": bp3,
        "callers": set(), "support": 0, "arriba_conf": "", "type": "", "known": False
    })
    d["callers"].add(caller)
    try: d["support"] += int(float(support))
    except Exception: pass
    # keep first non-empty breakpoints / directional gene order (prefer Arriba, then STAR-Fusion)
    if caller == "arriba" or (not d["bp5"] and bp5):
        d["genes"] = (g5, g3); d["bp5"] = bp5 or d["bp5"]; d["bp3"] = bp3 or d["bp3"]
    if extra.get("arriba_conf"): d["arriba_conf"] = extra["arriba_conf"]
    if extra.get("type"):        d["type"] = extra["type"]

store = {}

# ---- Arriba (tab; header starts with '#gene1') ----
try:
    with open("arriba.tsv") as fh:
        rows = [l.rstrip("\\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].lstrip("#").split("\\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in rows[1:]:
            if line.startswith("#"): continue
            f = line.split("\\t")
            def g(col, dflt=""):
                return f[idx[col]] if col in idx and idx[col] < len(f) else dflt
            g5, g3 = norm_gene(g("gene1")), norm_gene(g("gene2"))
            if not g5 or not g3: continue
            bp5 = f"{norm_chrom(parse_bp(g('breakpoint1'))[0])}:{parse_bp(g('breakpoint1'))[1]}"
            bp3 = f"{norm_chrom(parse_bp(g('breakpoint2'))[0])}:{parse_bp(g('breakpoint2'))[1]}"
            sup = 0
            for c in ("split_reads1", "split_reads2", "discordant_mates"):
                try: sup += int(g(c, "0") or 0)
                except Exception: pass
            add(store, canon_key(g5, g3), g5, g3, bp5, bp3, "arriba", sup,
                {"arriba_conf": g("confidence"), "type": g("type")})
except FileNotFoundError:
    pass

# ---- STAR-Fusion (tab; header starts with '#FusionName') ----
try:
    with open("starfusion.tsv") as fh:
        rows = [l.rstrip("\\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].lstrip("#").split("\\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in rows[1:]:
            if line.startswith("#"): continue
            f = line.split("\\t")
            def g(col, dflt=""):
                return f[idx[col]] if col in idx and idx[col] < len(f) else dflt
            fname = g("FusionName")
            if "--" in fname:
                g5, g3 = fname.split("--")[0], fname.split("--")[1]
            else:
                g5, g3 = g("LeftGene"), g("RightGene")
            g5, g3 = norm_gene(g5), norm_gene(g3)
            if not g5 or not g3: continue
            lb, rb = parse_bp(g("LeftBreakpoint")), parse_bp(g("RightBreakpoint"))
            bp5, bp3 = f"{lb[0]}:{lb[1]}", f"{rb[0]}:{rb[1]}"
            sup = 0
            for c in ("JunctionReadCount", "SpanningFragCount"):
                try: sup += int(float(g(c, "0") or 0))
                except Exception: pass
            add(store, canon_key(g5, g3), g5, g3, bp5, bp3, "starfusion", sup, {})
except FileNotFoundError:
    pass

# ---- FusionCatcher (tab; long header) ----
try:
    with open("fusioncatcher.txt") as fh:
        rows = [l.rstrip("\\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].split("\\t")
        def find(sub):
            for i, h in enumerate(hdr):
                if sub in h: return i
            return -1
        i_g1, i_g2 = find("Gene_1_symbol"), find("Gene_2_symbol")
        i_b1, i_b2 = find("Fusion_point_for_gene_1"), find("Fusion_point_for_gene_2")
        i_sp, i_su = find("Spanning_pairs"), find("Spanning_unique_reads")
        for line in rows[1:]:
            f = line.split("\\t")
            if i_g1 < 0 or i_g2 < 0 or max(i_g1, i_g2) >= len(f): continue
            g5, g3 = norm_gene(f[i_g1]), norm_gene(f[i_g2])
            if not g5 or not g3: continue
            b1 = parse_bp(f[i_b1]) if 0 <= i_b1 < len(f) else ("", "")
            b2 = parse_bp(f[i_b2]) if 0 <= i_b2 < len(f) else ("", "")
            bp5, bp3 = f"{b1[0]}:{b1[1]}", f"{b2[0]}:{b2[1]}"
            sup = 0
            for i_c in (i_sp, i_su):
                if 0 <= i_c < len(f):
                    try: sup += int(float(f[i_c] or 0))
                    except Exception: pass
            add(store, canon_key(g5, g3), g5, g3, bp5, bp3, "fusioncatcher", sup, {})
except FileNotFoundError:
    pass

# ---- known oncogenic fusions (order-insensitive pair set) ----
known = set()
if KNOWN_PATH != "NONE":
    try:
        with open(KNOWN_PATH) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = re.split(r'[\\t,;\\- ]+|--', line)
                parts = [p for p in parts if p]
                if len(parts) >= 2:
                    known.add("--".join(sorted([norm_gene(parts[0]), norm_gene(parts[1])])))
    except FileNotFoundError:
        pass

def same_chrom_close(bp5, bp3):
    try:
        c5, p5 = bp5.split(":"); c3, p3 = bp3.split(":")
        return c5 == c3 and c5 != "" and abs(int(p5) - int(p3)) <= READTHRU_BP
    except Exception:
        return False

# ---- verdicts ----
out_rows = []
n_pass = n_known = n_review = 0
for key, d in store.items():
    g5, g3 = d["genes"]
    ncall = len(d["callers"])
    is_known = key in known
    arriba_type = (d.get("type") or "").lower()
    readthrough = ("read_through" in arriba_type or "read-through" in arriba_type
                   or same_chrom_close(d["bp5"], d["bp3"]))
    if ncall >= MIN_CALLERS:
        verdict = "PASS_CONSENSUS"; n_pass += 1
    elif is_known and ncall >= 1:
        verdict = "PASS_KNOWN"; n_known += 1
    else:
        verdict = "REVIEW"; n_review += 1
    out_rows.append([
        f"{g5}--{g3}", g5, g3, d["bp5"], d["bp3"],
        ",".join(sorted(d["callers"])), str(ncall), str(d["support"]),
        d["arriba_conf"] or ".", "yes" if is_known else "no",
        "yes" if readthrough else "no", verdict
    ])

# sort: PASS_CONSENSUS, PASS_KNOWN, REVIEW; then by caller count / support desc
order = {"PASS_CONSENSUS": 0, "PASS_KNOWN": 1, "REVIEW": 2}
out_rows.sort(key=lambda r: (order.get(r[11], 3), -int(r[6]), -int(r[7])))

with open(f"{SAMPLE}_fusion_consensus.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\\t")
    w.writerow(["FusionName","Gene5","Gene3","Breakpoint5","Breakpoint3",
                "Callers","NumCallers","TotalSupport","Arriba_Confidence",
                "KnownFusion","ReadThrough","Verdict"])
    w.writerows(out_rows)

# FusionInspector candidate list = confident set only (PASS_*), directional "A--B"
with open(f"{SAMPLE}_fusion_candidates.txt", "w") as out:
    seen = set()
    for r in out_rows:
        if r[11].startswith("PASS") and r[0] not in seen:
            out.write(r[0] + "\\n"); seen.add(r[0])

with open(f"{SAMPLE}_fusion_consensus_stats.txt", "w") as out:
    out.write(f"sample_id\\t{SAMPLE}\\n")
    out.write(f"total_unique_fusions\\t{len(store)}\\n")
    out.write(f"pass_consensus\\t{n_pass}\\n")
    out.write(f"pass_known_single_caller\\t{n_known}\\n")
    out.write(f"review_only\\t{n_review}\\n")
    out.write(f"min_callers\\t{MIN_CALLERS}\\n")
PYEOF

    PY=\$(command -v python3 || echo /opt/conda/bin/python3)
    "\$PY" fusion_consensus.py ${min_callers} ${readthru_kb} ${known_arg} ${meta.sample_id}

    # never let an empty candidate list crash FusionInspector downstream
    if [ ! -s ${meta.sample_id}_fusion_candidates.txt ]; then
        printf '' > ${meta.sample_id}_fusion_candidates.txt
    fi
    """
}


// ============================================================
//  STAGE 5E: In-silico validation - FusionInspector  (analog of a "filter")
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
    set -euo pipefail

    # Early exit if input candidates file is missing or empty
    if [ ! -s "${candidates}" ]; then
        touch "${meta.sample_id}_filtered_candidates.txt"
        exit 0
    fi

    # 1. Extract all valid gene names/IDs directly from CTAT GTF
    gtf_file=\$(find -L "${ctat_lib}" -name "ref_annot.gtf" | head -n1)
    awk -F'\\t' '\$3 == "gene" || \$3 == "exon" { print \$9 }' "\$gtf_file" \\
        | grep -oE '(gene_id|gene_name|gene_symbol) "[^"]+"' \\
        | cut -d'"' -f2 \\
        | sort -u > valid_genes.txt

    # 2. Filter out rows containing unannotated genes
    awk -F'\\t' -v gtf="valid_genes.txt" '
        BEGIN { while ((getline line < gtf) > 0) valid[line] = 1 }
        /^#/ { print; next }
        {
            split(\$1, genes, "--")
            # If line is gene pair (GENE1--GENE2), ensure both exist in GTF
            if (length(genes) > 1) {
                if ((genes[1] in valid) && (genes[2] in valid)) print \$0
            } else if (\$1 in valid) {
                print \$0
            }
        }
    ' "${candidates}" > "${meta.sample_id}_filtered_candidates.txt"
    """
}

// ============================================================
//  STAGE 5E: In-silico validation - FusionInspector  (analog of a "filter")
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

    if [ -s fi_out/FI.fusion_predictions.final.abridged.tsv ]; then
        cp fi_out/FI.fusion_predictions.final.abridged.tsv ${meta.sample_id}_fusioninspector.tsv
    elif [ -s fi_out/finspector.FusionInspector.fusions.abridged.tsv ]; then
        cp fi_out/finspector.FusionInspector.fusions.abridged.tsv ${meta.sample_id}_fusioninspector.tsv
    else
        printf '#FusionName\\tJunctionReadCount\\tSpanningFragCount\\tFFPM\\tannots\\n' \\
            > ${meta.sample_id}_fusioninspector.tsv
    fi

    rm -rf "\$TMPDIR"
    """
}

// ============================================================
//  STAGE 6: Fusion annotation + pre-classification
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
          path(fusioninspector, stageAs: 'fi.tsv')
    path  cgc_tsv
    path  known_fusions_tsv

    output:
    tuple val(meta), path("${meta.sample_id}_fusions_annotated.tsv"), emit: annotated

    script:
    def cgc_arg    = cgc_tsv ? "cgc.tsv" : 'NONE'
    def known_arg  = known_fusions_tsv ? "known.tsv" : 'NONE'
    // Resolve .gz-ness in Groovy (avoids fragile shell filename tests).
    def cgc_prep   = !cgc_tsv ? '' :
                     (cgc_tsv.name.endsWith('.gz') ? "zcat ${cgc_tsv} > cgc.tsv" : "cp ${cgc_tsv} cgc.tsv")
    def known_prep = known_fusions_tsv ? "cp ${known_fusions_tsv} known.tsv" : ''
    def min_known  = (params.fusion_known_min_reads != null) ? params.fusion_known_min_reads : 3
    """
    set -euo pipefail
    ${cgc_prep}
    ${known_prep}

    cat > annotate_fusions.py <<'PYEOF'
#!/usr/bin/env python3
import sys, re

CGC_PATH   = sys.argv[1]
KNOWN_PATH = sys.argv[2]
SAMPLE     = sys.argv[3]

def norm_gene(g):
    g = str(g or "").strip().upper().split('^')[0]
    return re.sub(r'\\(.*?\\)', '', g).strip()

# CGC gene -> role (col layout matches the DNA CLASSIFY_VARIANTS usage:
# col1=symbol, col15=role_in_cancer). Robust to header presence.
role = {}
if CGC_PATH != "NONE":
    # Columns located by HEADER NAME (CGC layouts differ between releases).
    with open(CGC_PATH) as fh:
        lines = [l.rstrip("\\n").split("\\t") for l in fh if l.strip()]
    if lines:
        hdr_cgc = [re.sub(r"\\s+", "_", x.strip().upper()) for x in lines[0]]
        if "GENE_SYMBOL" not in hdr_cgc or "ROLE_IN_CANCER" not in hdr_cgc:
            sys.exit("ERROR: CGC header lacks GENE_SYMBOL / ROLE_IN_CANCER: " + " | ".join(lines[0]))
        gi, ri = hdr_cgc.index("GENE_SYMBOL"), hdr_cgc.index("ROLE_IN_CANCER")
        for f in lines[1:]:
            if len(f) > max(gi, ri) and f[gi]:
                role[norm_gene(f[gi])] = f[ri]

known = set()
if KNOWN_PATH != "NONE":
    try:
        with open(KNOWN_PATH) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = [p for p in re.split(r'[\\t,;\\- ]+|--', line) if p]
                if len(parts) >= 2:
                    known.add("--".join(sorted([norm_gene(parts[0]), norm_gene(parts[1])])))
    except FileNotFoundError:
        pass

# FusionInspector confirmation set (directional A--B, upper). fi.tsv is staged only
# when FusionInspector ran for this sample.
import os
FI_RAN    = os.path.exists("fi.tsv")
MIN_KNOWN = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
confirmed = {}
try:
    with open("fi.tsv") as fh:
        rows = [l.rstrip("\\n") for l in fh if l.strip()]
    if rows:
        hdr = rows[0].lstrip("#").split("\\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in rows[1:]:
            f = line.split("\\t")
            name = f[idx["FusionName"]] if "FusionName" in idx and idx["FusionName"] < len(f) else ""
            name = "--".join(norm_gene(x) for x in name.split("--")) if "--" in name else norm_gene(name)
            jrc = f[idx["JunctionReadCount"]] if "JunctionReadCount" in idx and idx["JunctionReadCount"] < len(f) else "0"
            confirmed[name] = jrc
except FileNotFoundError:
    pass

IG_HLA = re.compile(r'^(IGH|IGK|IGL|TRA|TRB|TRG|TRD|HLA)')

out = []
with open("consensus.tsv") as fh:
    rows = [l.rstrip("\\n") for l in fh if l.strip()]
hdr = rows[0].split("\\t") if rows else []
cidx = {h: i for i, h in enumerate(hdr)}
def col(f, name, d=""):
    return f[cidx[name]] if name in cidx and cidx[name] < len(f) else d

for line in rows[1:]:
    f = line.split("\\t")
    g5, g3 = norm_gene(col(f, "Gene5")), norm_gene(col(f, "Gene3"))
    key = "--".join(sorted([g5, g3]))
    dirname = f"{g5}--{g3}"
    verdict = col(f, "Verdict")
    readthrough = col(f, "ReadThrough") == "yes"
    is_known = key in known or col(f, "KnownFusion") == "yes"
    r5, r3 = role.get(g5, ""), role.get(g3, "")
    in_cgc = bool(r5) or bool(r3)
    fi_reads = confirmed.get(dirname, confirmed.get(key, ""))
    fi_confirmed = "yes" if (fi_reads not in ("", None) and str(fi_reads).strip() not in ("", "0")) else "no"
    ig_hla = bool(IG_HLA.match(g5) or IG_HLA.match(g3))

    # tiering (computational; for pathologist review - not a final AMP/CAP tier)
    # A catalogued pair needs real support: multi-caller consensus, OR a single-caller
    # known rescue with >= MIN_KNOWN reads AND FusionInspector confirmation when FI ran
    # (1-read single-caller hits are typical of index hopping / cross-sample bleed).
    try:
        support = float(col(f, "TotalSupport", "0") or 0)
    except ValueError:
        support = 0.0
    known_supported = (verdict == "PASS_CONSENSUS") or (
        verdict.startswith("PASS") and support >= MIN_KNOWN and (not FI_RAN or fi_confirmed == "yes"))
    if is_known and known_supported:
        tier = "KNOWN_ONCOGENIC"
    elif is_known:
        tier = "KNOWN_LOW_SUPPORT"
    elif readthrough or ig_hla:
        tier = "LIKELY_ARTIFACT"
    elif in_cgc and verdict == "PASS_CONSENSUS":
        tier = "LIKELY_ONCOGENIC"
    elif verdict == "PASS_CONSENSUS":
        tier = "VUS"
    else:
        tier = "LOW_CONFIDENCE"

    out.append([
        dirname, g5, g3, col(f, "Breakpoint5"), col(f, "Breakpoint3"),
        col(f, "Callers"), col(f, "NumCallers"), col(f, "TotalSupport"),
        col(f, "Arriba_Confidence"), "yes" if is_known else "no",
        "yes" if in_cgc else "no", (r5 or ".") + "|" + (r3 or "."),
        "yes" if readthrough else "no", "yes" if ig_hla else "no",
        fi_confirmed, fi_reads or ".", verdict, tier
    ])

order = {"KNOWN_ONCOGENIC":0,"LIKELY_ONCOGENIC":1,"KNOWN_LOW_SUPPORT":2,
         "VUS":3,"LOW_CONFIDENCE":4,"LIKELY_ARTIFACT":5}
out.sort(key=lambda r: (order.get(r[17], 9), r[0]))

with open(f"{SAMPLE}_fusions_annotated.tsv", "w") as o:
    o.write("# COMPUTATIONAL FUSION PRE-CLASSIFICATION - NOT A FINAL AMP/ASCO/CAP CLINICAL TIER\\n")
    o.write("# Final tiering requires a curated knowledge base (OncoKB/CIViC/COSMIC) + expert review.\\n")
    o.write("\\t".join([
        "FusionName","Gene5","Gene3","Breakpoint5","Breakpoint3","Callers","NumCallers",
        "TotalSupport","Arriba_Confidence","KnownFusion","In_CGC","CGC_Roles(5|3)",
        "ReadThrough","IG_HLA_Locus","FI_Confirmed","FI_JunctionReads","Verdict","TIER"
    ]) + "\\n")
    for r in out:
        o.write("\\t".join(r) + "\\n")
PYEOF

    PY=\$(command -v python3 || echo /opt/conda/bin/python3)
    "\$PY" annotate_fusions.py ${cgc_arg} ${known_arg} ${meta.sample_id} ${min_known}
    """
}


// ============================================================
//  STAGE 7: Per-patient consolidated fusion report
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
    set -uo pipefail
    REPORT=${patient_id}_fusion_report.tsv

    {
      echo "## ============================================================"
      echo "## CONSOLIDATED RNA FUSION REPORT - PATIENT: ${patient_id}   (RNA sample: ${sample_id})"
      echo "## Generated: \$(date '+%Y-%m-%d %H:%M:%S')"
      echo "## NOTE: computational summary for review; not a final clinical tier."
      echo "## ============================================================"
    } > \$REPORT

    printf '\\n## [1] CALLING SUMMARY\\n' >> \$REPORT
    if [ -s consensus_stats.txt ]; then cat consensus_stats.txt >> \$REPORT; else echo "(no stats)" >> \$REPORT; fi

    printf '\\n## [2] PRIORITISED FUSIONS (KNOWN_ONCOGENIC / LIKELY_ONCOGENIC first)\\n' >> \$REPORT
    if [ -s annotated.tsv ]; then
        # header + high-tier calls first, then the rest, artifacts last
        grep -E '^FusionName' annotated.tsv >> \$REPORT || true
        grep -E 'KNOWN_ONCOGENIC|LIKELY_ONCOGENIC' annotated.tsv >> \$REPORT || true
    else
        echo "(no fusion calls)" >> \$REPORT
    fi

    printf '\\n## [3] ALL CONSENSUS FUSIONS (full annotated table)\\n' >> \$REPORT
    if [ -s annotated.tsv ]; then grep -v '^#' annotated.tsv >> \$REPORT; else echo "(none)" >> \$REPORT; fi
    """
}
