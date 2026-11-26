#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
//  workflows.nf  -  High-level Workflow Orchestration
// ============================================================


//  ------------------- stage 0 -------------------------------
//  ------------------------ utils ----------------------------
include { GENERATE_META_INFO          } from './utils'
include { stageRankMap; stageRuns; naturalCompare } from './utils'

//  ---------------- components_common ------------------------
include { MANIFEST_ROW                } from './components_common'
include { PREFLIGHT_REFERENCES        } from './components_common'
include { WRITE_MANIFEST              } from './components_common'

//  --------------- multistage processes ----------------------
include { MULTIQC as MULTIQC_TRIMMING } from './components_common'
include { MULTIQC as MULTIQC_MARKDUP  } from './components_common'
include { MULTIQC as MULTIQC_METRICS  } from './components_common'
include { MULTIQC as MULTIQC_SOMALIER } from './components_common'
include { MULTIQC as MULTIQC_ANNOT    } from './components_common'
include { CRAM_CONVERT as CRAM_CONVERT_DEDUP } from './components_common'
include { CRAM_CONVERT as CRAM_CONVERT_BQSR  } from './components_common'
//  ------------------- stage 1 -------------------------------
//  ---------------- components_common ------------------------
include { PREPROCESS_AND_ALIGN_DNA    } from './components_common'
include { ALIGN_DNA_UMI               } from './components_common'

//  ------------------- stage 2 -------------------------------
//  ---------------- components_common ------------------------
include { MERGE_AND_MARK_DUPLICATES   } from './components_common'
include { MERGE_AND_CONSENSUS         } from './components_common'
include { MERGE_AND_MARKDUP_UMI       } from './components_common'

//  ------------------- stage 3 -------------------------------
//  ---------------- components_common ------------------------
include { MOSDEPTH                    } from './components_common'
include { INFER_SEX                   } from './components_common'
include { SAMTOOLS_STATS              } from './components_common'
include { METRICS                     } from './components_common'
include { VERIFYBAMID2                } from './components_common'
include { QC_GATE;QC_GATE_SUMMARY     } from './components_common'
include { SOMALIER_EXTRACT            } from './components_common'
include { SOMALIER_ANCESTRY           } from './components_common'
include { SOMALIER_RELATE             } from './components_common'
include { CHECK_SAMPLE_SWAPS          } from './components_common'

//  ------------------- stage 4 -------------------------------
//  ---------------- components_common ------------------------
include { BASE_RECALIBRATOR           } from './components_common'
include { GATHER_BQSR_REPORTS         } from './components_common'
include { APPLY_BQSR                  } from './components_common'

//  ------------------- stage 5 -------------------------------
//  ---------------- components_common ------------------------
include { MUTATIONAL_SIGNATURES       } from './components_common'
include { CNVKIT                      } from './components_common'
include { CNVKIT_CALL                 } from './components_common'
include { MSISENSOR_PRO;CLASSIFY_MSI  } from './components_common'
include { MSISENSOR_GENE_MAP          } from './components_common'
include { GETPILEUPSUMMARIES          } from './components_common'
include { CALCULATE_CONTAMINATION     } from './components_common'
include { LEARN_READ_ORIENTATION_MODEL } from './components_common'
include { EVIDENCE_BAM                } from './components_common'
//  --------------- components_germline -----------------------
include { HAPLOTYPECALLER             } from './components_germline'
include { MERGE_GVCFS                 } from './components_germline'
include { GENOTYPE_GVCFS              } from './components_germline'
include { NV_SCORE_VARIANTS           } from './components_germline'
include { FILTER_VARIANT_TRANCHES     } from './components_germline'
include { GENOMICSDB_IMPORT           } from './components_germline'
include { GENOTYPE_GVCFS_JOINT        } from './components_germline'
include { REBLOCK_GVCF                } from './components_germline'
include { GNARLY_GENOTYPER            } from './components_germline'
include { MERGE_COHORT_VCF            } from './components_germline'
include { VARIANT_RECALIBRATOR as VQSR_SNP   } from './components_germline'
include { VARIANT_RECALIBRATOR as VQSR_INDEL } from './components_germline'
include { APPLY_VQSR as APPLY_VQSR_SNP       } from './components_germline'
include { APPLY_VQSR as APPLY_VQSR_INDEL     } from './components_germline'
include { SPLIT_COHORT_VCF            } from './components_germline'
include { DEEPVARIANT                 } from './components_germline'
include { MANTA_GERMLINE              } from './components_germline'
include { GATK_COLLECT_READ_COUNTS    } from './components_germline'
include { GATK_PLOIDY_MODEL           } from './components_germline'
include { GATK_GCNV_CALLER            } from './components_germline'
include { GATK_GCNV_POSTPROCESS       } from './components_germline'
include { BCFTOOLS_ROH                } from './components_germline'
include { OPTITYPE                    } from './components_germline'
//  --------------- components_somatic ------------------------
include { MANTA                       } from './components_somatic'
include { STRELKA2;TAG_STRELKA_AF     } from './components_somatic'
include { MUSE_CALL;MUSE_SUMP         } from './components_somatic'
include { TAG_MUSE_AF                 } from './components_somatic'
include { VARSCAN;TAG_VARSCAN_AF      } from './components_somatic'
include { VARSCAN_MERGE               } from './components_somatic'
include { MUTECT2;MERGE_MUTECT2       } from './components_somatic'
include { ENSEMBLE_CONSENSUS          } from './components_somatic'
include { FILTER_SOMATIC_VCF          } from './components_somatic'
include { SEQUENZA                    } from './components_somatic'
include { SEQUENZA_GENE_LOH           } from './components_somatic'

//  ------------------- stage 6 -------------------------------
//  ---------------- components_common ------------------------
include { ANNOTATE_VCF                } from './components_common'
include { ANNOTATE_VCF as ANNOTATE_VCF_GERMLINE } from './components_common'
include { COMPUTE_TMB                 } from './components_common'
include { ANNOTATE_SV as ANNOTATE_SV_SOMATIC } from './components_common'
include { ANNOTATE_SV as ANNOTATE_SV_MANTA } from './components_common'
include { ANNOTATE_SV as ANNOTATE_SV_GCNV  } from './components_common'

//  ------------------- stage 7 -------------------------------
//  ---------------- components_common ------------------------
include { CLASSIFY_VARIANTS           } from './components_common'
include { VCF2MAF                     } from './components_common'

//  --------------- components_somatic ------------------------
include { CONSOLIDATE_PATIENT_REPORT  } from './components_somatic'

//  --------------- components_germline ------------------------
include { CONSOLIDATE_GERMLINE_REPORT } from './components_germline'

// ============================================================
//  TOP-LEVEL HELPERS
// ============================================================
def stripChunkMeta(Map meta) {
    return meta.findAll { k, _v -> !(k in ['uid', 'rg']) }
}

def toFileChannel(String path) {
    channel.value(file(path, checkIfExists: true))
}

def toOptionalFileChannel(path) {
    path ? channel.value(file(path, checkIfExists: true)) : channel.value([])
}


// ============================================================
//  SUB-WORKFLOW: RUN_BQSR
// ============================================================
workflow RUN_BQSR {
    take:
    ch_bam_bai
    ch_ref_fasta
    ch_ref_fai
    ch_ref_dict
    ch_dbsnp
    ch_dbsnp_tbi
    ch_mills
    ch_mills_tbi
    ch_known_indels
    ch_known_indels_tbi
    ch_interval_bed
    ch_scatter_intervals

    main:
    // WGS: one BaseRecalibrator task per scatter shard, tables gathered (GATK WGS
    // practice). WES: one task on the capture BED, exactly as before.
    ch_bqsr_shards = (params.assay == 'wgs')
        ? ch_scatter_intervals
        : channel.of(tuple('all', []))

    BASE_RECALIBRATOR(
        ch_bam_bai.combine(ch_bqsr_shards),
        ch_ref_fasta, ch_ref_fai, ch_ref_dict,
        ch_dbsnp, ch_dbsnp_tbi,
        ch_mills, ch_mills_tbi,
        ch_known_indels, ch_known_indels_tbi,
        ch_interval_bed
    )

    GATHER_BQSR_REPORTS(
        BASE_RECALIBRATOR.out.recal_table
            .map { meta, tbl -> tuple(meta.sample_id, meta, tbl) }
            .groupTuple(by: 0)
            .map { _sid, metas, tbls -> tuple(metas[0], tbls.sort { t -> t.name }) }
    )

    ch_bqsr_input = ch_bam_bai
        .map { meta, bam, bai -> tuple(meta.sample_id, meta, bam, bai) }
        .join(GATHER_BQSR_REPORTS.out.recal_table.map { meta, tbl -> tuple(meta.sample_id, tbl) }, by: 0)
        .map { _sid, meta, bam, bai, tbl -> tuple(meta, bam, bai, tbl) }

    APPLY_BQSR(
        ch_bqsr_input,
        ch_ref_fasta, ch_ref_fai, ch_ref_dict,
        ch_interval_bed
    )

    emit:
    bam_bai_pair = APPLY_BQSR.out.bam_bai_pair
}


// ============================================================
//  SUB-WORKFLOW: QC_BLOCK
// ============================================================
workflow QC_BLOCK {
    take:
    ch_bam_bai
    ch_ref_fasta
    ch_ref_fai
    ch_ref_dict
    ch_interval_bed
    ch_bait_interval_list
    ch_target_interval_list
    ch_verifybamid_db
    ch_coverage_report_bed

    main:
    MOSDEPTH(ch_bam_bai, ch_interval_bed, ch_coverage_report_bed)
    INFER_SEX(MOSDEPTH.out.summary)
    SAMTOOLS_STATS(ch_bam_bai, ch_ref_fasta, ch_ref_fai)
    METRICS(ch_bam_bai, ch_ref_fasta, ch_ref_fai, ch_ref_dict, ch_bait_interval_list, ch_target_interval_list)

    ch_qc_mixed = MOSDEPTH.out.summary.map          { _m, f -> f }
        .mix( MOSDEPTH.out.global_dist.map           { _m, f -> f } )
        .mix( MOSDEPTH.out.region_dist.map           { _m, f -> f } )
        .mix( SAMTOOLS_STATS.out.stats.map           { _m, f -> f } )
        .mix( METRICS.out.alignment.map              { _m, f -> f } )
        .mix( METRICS.out.insert_size.map            { _m, f -> f } )
        .mix( METRICS.out.extra_metrics.flatMap      { _m, fs -> fs } )
        .mix(
            METRICS.out.hs_metrics
                .filter { meta, _f -> meta.assay == 'wes' }
                .map    { _m, f -> f }
        )
    // Contamination QC is OPTIONAL - only runs if the VerifyBamID2
    // resource bundle is provided (params.verifybamid_db). Without it,
    // this step is skipped so stages 1-4 run with reference + BQSR files only.
    ch_contamination = channel.empty()
    if (params.verifybamid_db) {
        // Normals only: tumour aneuploidy/LOH biases FREEMIX. Tumour contamination
        // comes from CalculateContamination (somatic block).
        VERIFYBAMID2(ch_bam_bai.filter { m, _b, _i -> m.type == 'DN' },
                     ch_verifybamid_db, ch_ref_fasta, ch_ref_fai)
        ch_contamination = VERIFYBAMID2.out.contamination
        ch_qc_mixed = ch_qc_mixed.mix( VERIFYBAMID2.out.contamination.map { _m, f -> f } )
    } else {
        log.info "Contamination QC skipped (params.verifybamid_db not set)."
    }

    if (params.qc_gate) {
        // Use ch_contamination (empty-channel default) NOT VERIFYBAMID2.out - the
        // latter does not exist when params.verifybamid_db is unset.

        ch_qc_gate_input = MOSDEPTH.out.summary
            .map { meta, summ -> tuple(meta.sample_id, meta, summ) }
            .join( ch_contamination.map { meta, sm -> tuple(meta.sample_id, sm) }, by: 0, remainder: true )
            .join( METRICS.out.hs_metrics.map { meta, hs -> tuple(meta.sample_id, hs) }, by: 0, remainder: true )
            .join( SAMTOOLS_STATS.out.stats.map { meta, st -> tuple(meta.sample_id, st) }, by: 0, remainder: true )
            .join( MOSDEPTH.out.global_dist.map { meta, gd -> tuple(meta.sample_id, gd) }, by: 0, remainder: true )
            .map { row ->
                def meta = row[1]; def summ = row[2]
                def sm = row.size()>3 ? (row[3]?:[]) : []
                def hs = row.size()>4 ? (row[4]?:[]) : []
                def st = row.size()>5 ? (row[5]?:[]) : []
                def gd = row.size()>6 ? (row[6]?:[]) : []
                tuple(meta, summ, sm, hs, st, gd)
            }
        QC_GATE(ch_qc_gate_input)
        QC_GATE_SUMMARY( QC_GATE.out.report.map { _meta, tsv -> tsv }.collect(sort: true) )
    }

    ch_qc_mixed
        .collect(sort: true)
        .map { files -> tuple("03_metrics", files) }
        | MULTIQC_METRICS

    emit:
    mosdepth_summary = MOSDEPTH.out.summary
    depth_quantized  = MOSDEPTH.out.quantized
    contamination    = ch_contamination
}


// ============================================================
//  SUB-WORKFLOW: GERMLINE_BLOCK
//  SNP/INDEL + SV (Manta) + CNV (gCNV)
// ============================================================
workflow GERMLINE_BLOCK {
    take:
    ch_bam_bai
    ch_ref_fasta
    ch_ref_fai
    ch_ref_dict
    ch_dbsnp
    ch_dbsnp_tbi
    ch_interval_bed
    //ch_interval_list
    ch_call_regions_gz       // value: bgzipped+tabixed BED tuple, or []
    ch_call_regions_tbi
    ch_hapmap
    ch_hapmap_tbi
    ch_omni
    ch_omni_tbi
    ch_thousandg_snps
    ch_thousandg_tbi
    ch_mills
    ch_mills_tbi
    ch_preprocessed_intervals
    ch_contig_ploidy_priors
    ch_gcnv_scatter_intervals
    ch_gcnv_annotated_intervals
    ch_germ_scatter
    ch_roh_af_file
    ch_roh_af_file_tbi

    main:
    ch_snv_vcf = channel.empty()
    ch_gvcf    = channel.empty()
    ch_roh     = channel.empty()

    //DN = DNA normal
    ch_germ_bam_bai = ch_bam_bai.filter { meta, _b, _i -> meta.type == 'DN' }


    // -- SNP / INDEL -------------------------------------------
    if (params.germline_caller == 'deepvariant') {
        DEEPVARIANT(ch_germ_bam_bai, ch_ref_fasta, ch_ref_fai, ch_interval_bed)
        ch_snv_vcf = DEEPVARIANT.out.vcf_gvcf.map { meta, vcf, tbi, _g, _gt -> tuple(meta, vcf, tbi) }
        ch_gvcf    = DEEPVARIANT.out.vcf_gvcf.map { meta, _v, _t, gvcf, gtbi -> tuple(meta, gvcf, gtbi) }

    } else {
        ch_hc_input = ch_germ_bam_bai.combine(ch_germ_scatter)

        HAPLOTYPECALLER(ch_hc_input, ch_ref_fasta, ch_ref_fai, ch_ref_dict, ch_dbsnp, ch_dbsnp_tbi)

        ch_gvcf_chunks = HAPLOTYPECALLER.out.gvcf_chunk
            .map { meta, _iname, gvcf, tbi -> tuple(meta.sample_id, meta, gvcf, tbi) }
            .groupTuple(by: 0)
            .map { _sid, metas, gvcfs, tbis -> tuple(metas[0], gvcfs.sort { g -> g.name }, tbis.sort { t -> t.name }) }

        MERGE_GVCFS(ch_gvcf_chunks, ch_ref_fasta, ch_ref_fai, ch_ref_dict)
        ch_gvcf = MERGE_GVCFS.out.gvcf



        if (params.germline_filter == 'vqsr') {

            // Cohort-size guard (side channel, resume-safe): VQSR cannot train a
            // reliable model on a small callset. Abort with guidance rather than
            // emit mis-calibrated tranches.
            ch_germ_bam_bai
                .count()
                .map { n ->
                    def min_vqsr = params.vqsr_min_samples ?: 30
                    if (n < min_vqsr) {
                        error "germline_filter='vqsr' needs >= ${min_vqsr} samples to train a " +
                              "reliable model, but only ${n} were provided. Use " +
                              "params.germline_filter='nv' (single-sample, robust at any size) " +
                              "for this cohort, or add samples."
                    }
                    log.info "VQSR: ${n} samples - training cohort model."
                    return n
                }

            // Reblock per sample (GnarlyGenotyper
            // REQUIRES the -do-qual-approx annotation ReblockGVCF adds).
            def use_gnarly = (params.joint_genotyper ?: 'classic') == 'gnarly'

            def do_reblock = use_gnarly || (params.reblock_gvcf != false)

            if (do_reblock) {
                REBLOCK_GVCF(ch_gvcf, ch_ref_fasta, ch_ref_fai, ch_ref_dict)
                ch_gvcf_for_db = REBLOCK_GVCF.out.gvcf
            } else {
                ch_gvcf_for_db = ch_gvcf
            }

            ch_all_gvcfs = ch_gvcf_for_db.map { _m, g, _t -> g }.collect(sort: true)
            ch_all_tbis  = ch_gvcf_for_db.map { _m, _g, t -> t }.collect(sort: true)
            ch_gdb_shards = ch_germ_scatter

            GENOMICSDB_IMPORT(
                ch_gdb_shards,
                ch_all_gvcfs, ch_all_tbis,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict
            )

            if (use_gnarly) {
                GNARLY_GENOTYPER(
                    GENOMICSDB_IMPORT.out.gendb,
                    ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                    ch_dbsnp, ch_dbsnp_tbi
                )
                ch_joint_vcf = GNARLY_GENOTYPER.out.vcf
            } else {
                GENOTYPE_GVCFS_JOINT(
                    GENOMICSDB_IMPORT.out.gendb,
                    ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                    ch_dbsnp, ch_dbsnp_tbi
                )
                ch_joint_vcf = GENOTYPE_GVCFS_JOINT.out.vcf
            }

            MERGE_COHORT_VCF(
                ch_joint_vcf.map { _n, vcf, _t -> vcf }.collect(sort: true),
                ch_joint_vcf.map { _n, _v, tbi -> tbi }.collect(sort: true)
            )

            // VQSR: build SNP + INDEL models on the cohort VCF, then apply
            // SNP first and INDEL on the SNP-recalibrated output.
            ch_cohort_snp   = MERGE_COHORT_VCF.out.vcf.map { vcf, tbi -> tuple('SNP',   vcf, tbi) }
            ch_cohort_indel = MERGE_COHORT_VCF.out.vcf.map { vcf, tbi -> tuple('INDEL', vcf, tbi) }

            VQSR_SNP(
                ch_cohort_snp,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_hapmap, ch_hapmap_tbi, ch_omni, ch_omni_tbi,
                ch_thousandg_snps, ch_thousandg_tbi, ch_mills, ch_mills_tbi,
                ch_dbsnp, ch_dbsnp_tbi
            )
            VQSR_INDEL(
                ch_cohort_indel,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_hapmap, ch_hapmap_tbi, ch_omni, ch_omni_tbi,
                ch_thousandg_snps, ch_thousandg_tbi, ch_mills, ch_mills_tbi,
                ch_dbsnp, ch_dbsnp_tbi
            )

            APPLY_VQSR_SNP(
                ch_cohort_snp.join(VQSR_SNP.out.model, by: 0),
                ch_ref_fasta, ch_ref_fai, ch_ref_dict
            )
            // Feed the SNP-recalibrated cohort VCF into the INDEL apply step.
            ch_indel_apply_in = APPLY_VQSR_SNP.out.vcf
                .map { _mode, vcf, tbi -> tuple('INDEL', vcf, tbi) }
                .join(VQSR_INDEL.out.model, by: 0)
            APPLY_VQSR_INDEL(
                ch_indel_apply_in,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict
            )

            // Split the final cohort callset back to per-sample VCFs.
            ch_cohort_final = APPLY_VQSR_INDEL.out.vcf.map { _mode, vcf, tbi -> tuple(vcf, tbi) }.first()
            ch_split_in = ch_gvcf
                .map { meta, _g, _t -> meta }
                .combine(ch_cohort_final)
                .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }

            SPLIT_COHORT_VCF(ch_split_in, ch_ref_fasta, ch_ref_fai, ch_ref_dict)
            ch_snv_vcf = SPLIT_COHORT_VCF.out.vcf

        } else if (params.germline_filter == 'none') {
            // -- SINGLE-SAMPLE, UNFILTERED: GenotypeGVCFs output as the final callset.
            //    On WES the WGS-tuned CNN tranches lowered F1 in our HG002 benchmark.
            GENOTYPE_GVCFS(
                ch_gvcf,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_dbsnp, ch_dbsnp_tbi,
                ch_interval_bed
            )
            ch_snv_vcf = GENOTYPE_GVCFS.out.vcf

        } else {
            // -- SINGLE-SAMPLE PATH: per-sample genotype + NV filter --
            GENOTYPE_GVCFS(
                ch_gvcf,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_dbsnp, ch_dbsnp_tbi,
                ch_interval_bed
            )

            ch_nv_input = GENOTYPE_GVCFS.out.vcf
                .join(ch_germ_bam_bai, by: 0)
                .map { meta, vcf, tbi, bam, bai -> tuple(meta, vcf, tbi, bam, bai) }

            NV_SCORE_VARIANTS(
                ch_nv_input.map { meta, vcf, tbi, _b, _bi -> tuple(meta, vcf, tbi) },
                ch_nv_input.map { _m, _v, _t, bam, _bi -> bam },
                ch_nv_input.map { _m, _v, _t, _b, bai  -> bai },
                ch_ref_fasta, ch_ref_fai, ch_ref_dict
            )

            FILTER_VARIANT_TRANCHES(
                NV_SCORE_VARIANTS.out.vcf,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_hapmap, ch_hapmap_tbi,
                ch_omni, ch_omni_tbi,
                ch_thousandg_snps, ch_thousandg_tbi,
                ch_mills, ch_mills_tbi
            )
            ch_snv_vcf = FILTER_VARIANT_TRANCHES.out.vcf
        }
    }


    // -- ROH - runs of homozygosity (Phase 2) ------------------
    //   Runs on the final per-sample SNV/indel callset.
    //   Reliable only on dense (WGS) callsets WITH a population-AF file.
    //   WES marker density is too sparse for the HMM, and self-estimated AF
    //   from a single sample is not meaningful - so ROH is gated on both.
    ch_roh = channel.empty()
    if (params.run_roh) {
        def roh_ok      = true
        def roh_reason  = ''
        if (params.assay == 'wes' && !params.roh_allow_wes) {
            roh_ok = false
            roh_reason = "assay=wes (exome marker density is too sparse for reliable ROH; " +
                         "set params.roh_allow_wes=true to override)"
        } else if (!params.roh_af_file) {
            roh_ok = false
            roh_reason = "no roh_af_file provided (single-sample self-estimated AF is unreliable; " +
                         "supply a population-AF file, e.g. gnomAD)"
        }

        if (roh_ok) {
            BCFTOOLS_ROH(ch_snv_vcf, ch_roh_af_file, ch_roh_af_file_tbi)
            ch_roh = BCFTOOLS_ROH.out.roh
        } else {
            log.warn "Germline ROH skipped: ${roh_reason}."
        }
    } else {
        log.info "Germline ROH skipped (params.run_roh = false)."
    }


    // -- HLA class I typing (OptiType) on normals ---------------
    if (params.run_optitype) {
        OPTITYPE(ch_germ_bam_bai, ch_ref_fasta, ch_ref_fai)
    }

    // -- SV - Manta germline -----------------------------------
    MANTA_GERMLINE(
        ch_germ_bam_bai,
        ch_ref_fasta, ch_ref_fai,
        ch_call_regions_gz, ch_call_regions_tbi
    )

    // -- CNV - GATK gCNV (cohort mode; requires a minimum number of normals) --
    ch_cnv_vcf = channel.empty()

    // Somatic runs call normals too (germline_on_normals); the gCNV cohort model
    // needs >= gcnv_min_samples normals, so it is opt-in there (gcnv_in_somatic).
    def gcnv_on = params.run_gcnv && (params.analysis == 'germline' || params.gcnv_in_somatic)
    if (gcnv_on) {

        // Cohort-size guard as a SIDE channel: it counts normals and aborts if the
        // cohort is too small, WITHOUT sitting between ch_germ_bam_bai and the gCNV
        // processes - so process input hashes are unchanged and -resume still works.
        ch_germ_bam_bai
            .count()
            .map { n ->
                if (n < params.gcnv_min_samples) {
                    error "gCNV cohort mode requires >= ${params.gcnv_min_samples} normals, " +
                          "but only ${n} were provided. Add normals, lower " +
                          "params.gcnv_min_samples (not recommended below ~20-30), or set " +
                          "params.run_gcnv = false to skip germline CNV."
                }
                log.info "gCNV: ${n} normals - cohort model will be built."
                return n
            }

        GATK_COLLECT_READ_COUNTS(
            ch_germ_bam_bai,
            ch_ref_fasta, ch_ref_fai, ch_ref_dict,
            ch_preprocessed_intervals
        )

        // sort counts by sample_id for deterministic order.
        // The SAME sorted order is reused to assign --sample-index.
        ch_counts_sorted = GATK_COLLECT_READ_COUNTS.out.counts
            .toSortedList { a, b -> a[0].sample_id <=> b[0].sample_id }
            // emits a single list: [ [meta,counts], [meta,counts], ... ]

        // Value channel of just the count files, in sorted order.
        //  value channel - broadcasts to every shard task.
        ch_all_counts = ch_counts_sorted
            .map { rows -> rows.collect { _meta, counts -> counts } }

        GATK_PLOIDY_MODEL(
            ch_all_counts,
            ch_ref_fasta, ch_ref_fai, ch_ref_dict,
            ch_preprocessed_intervals,
            ch_contig_ploidy_priors
        )

        // ploidy_calls and counts as value channels so the
        // shard scatter (queue) drives N executions, broadcasting the
        // single-value inputs to each.
        ch_ploidy_calls_val = GATK_PLOIDY_MODEL.out.ploidy_calls.first()
        ch_all_counts_val   = ch_all_counts.first()

        ch_gcnv_shards = ch_gcnv_scatter_intervals

        GATK_GCNV_CALLER(
            ch_gcnv_shards.map { shard_name, _f -> shard_name },
            ch_all_counts_val,
            ch_ploidy_calls_val,
            ch_preprocessed_intervals,
            ch_gcnv_shards.map { _n, f -> f },
            ch_gcnv_annotated_intervals
        )

        // Calls and models taken from ONE list sorted by shard name, so the i-th
        // --calls-shard-path and i-th --model-shard-path are the same shard, in genomic
        // (natural shard-name) order, and the POSTPROCESS input hash is stable for -resume.
        ch_gcnv_out_sorted = GATK_GCNV_CALLER.out.shard_output
            .toSortedList { a, b -> naturalCompare(a[0].toString(), b[0].toString()) }
        ch_gcnv_calls  = ch_gcnv_out_sorted.map { rows -> rows.collect { r -> r[2] } }.first()
        ch_gcnv_models = ch_gcnv_out_sorted.map { rows -> rows.collect { r -> r[1] } }.first()

        // build (meta, index) from the SAME sorted list used
        // for ch_all_counts - guarantees --sample-index matches the
        // positional order GermlineCNVCaller saw the count files.
        ch_sample_indexed = ch_counts_sorted
            .flatMap { rows ->
                rows.withIndex().collect { row, idx -> tuple(row[0], idx) }
            }

        GATK_GCNV_POSTPROCESS(
            ch_sample_indexed,
            ch_gcnv_calls,
            ch_gcnv_models,
            ch_ploidy_calls_val,
            ch_ref_fasta, ch_ref_fai, ch_ref_dict
        )
        ch_cnv_vcf = GATK_GCNV_POSTPROCESS.out.cnv_vcf

    } else {
        log.info "Germline gCNV skipped (run_gcnv=${params.run_gcnv}, analysis=${params.analysis}, " +
                 "gcnv_in_somatic=${params.gcnv_in_somatic}). gCNV needs >=${params.gcnv_min_samples} normals."
    }

    emit:
    snv_vcf = ch_snv_vcf
    gvcf    = ch_gvcf
    sv_vcf  = MANTA_GERMLINE.out.sv_vcf
    cnv_vcf = ch_cnv_vcf
    roh     = ch_roh
}


// ============================================================
//  SUB-WORKFLOW: SOMATIC_BLOCK
//  SNP/INDEL ensemble + SV + CNV + MSI
// ============================================================
workflow SOMATIC_BLOCK {
    take:
    ch_bam_bai
    ch_ref_fasta
    ch_ref_fai
    ch_ref_dict
    ch_germline_resource
    ch_germline_resource_tbi
    ch_panel_of_normals
    ch_panel_of_normals_tbi
    ch_interval_bed
    ch_call_regions_gz
    ch_call_regions_tbi
    ch_scatter_intervals
    ch_cnvkit_pon
    ch_msisensor_list
    ch_pileup_variants
    ch_pileup_variants_tbi
    ch_sequenza_gc_wiggle
    ch_gene_bed
    ch_dbsnp
    ch_dbsnp_tbi

    main:
    // DT = DNA tumor, DN = DNA normal, RN = RNAseq (excluded from DNA somatic).
    ch_tumor = ch_bam_bai.filter { meta, _b, _i -> meta.type == 'DT' }
    ch_normal = ch_bam_bai.filter { meta, _b, _i -> meta.type == 'DN' }

    // combine (not join): every tumour of a patient pairs with that patient's normal
    // (multi-region / longitudinal). join would pair only one tumour per key.
    // normal_id = the normal's read-group SM, used by Mutect2 -normal and the filters.
    ch_paired = ch_tumor
        .map { meta, bam, bai -> tuple(meta.patient_id, meta, bam, bai) }
        .combine(
            ch_normal.map { meta, bam, bai -> tuple(meta.patient_id, meta.sample_id, bam, bai) },
            by: 0
        )
        .map { _pid, t_meta, t_bam, t_bai, n_id, n_bam, n_bai ->
            tuple(t_meta + [normal_id: n_id], t_bam, t_bai, n_bam, n_bai)
        }

    // Side channel (does not feed any process): the join above drops a tumour whose
    // patient has no DN sample WITHOUT any message. Make that visible.
    ch_tumor.map { meta, _b, _i -> tuple(meta.patient_id, 'DT') }
        .mix(ch_normal.map { meta, _b, _i -> tuple(meta.patient_id, 'DN') })
        .groupTuple(by: 0)
        .map { pid, types ->
            if (!('DN' in types))
                log.warn "Patient ${pid}: tumour (DT) has no matched normal (DN) - EXCLUDED from somatic calling."
            pid
        }

    // SV first - feeds Strelka2 indel candidates
    MANTA(
        ch_paired,
        ch_ref_fasta, ch_ref_fai,
        ch_call_regions_gz, ch_call_regions_tbi
    )

    // Strelka2 only runs when ensemble calling is enabled. With
    // use_strelka2_in_consensus=false the final SNV VCF comes from Mutect2 alone,
    // so Strelka2 is skipped entirely (avoids running it needlessly).

    // Join Manta indel candidates to each pair BY sample_id (scalar key - robust;
    // joining on the whole meta Map can silently fail if object identity differs).
    // All STRELKA2 inputs are then derived from this ONE joined channel, so the
    // BAMs and the matching indel candidates can never cross-pair between samples.
    if (params.use_strelka2_in_consensus) {
        ch_strelka_joined = ch_paired
            .map { meta, t_bam, t_bai, n_bam, n_bai -> tuple(meta.sample_id, meta, t_bam, t_bai, n_bam, n_bai) }
            .join(
                MANTA.out.indel_candidates.map { meta, indel, tbi -> tuple(meta.sample_id, indel, tbi) },
                by: 0
            )
            .map { _sid, meta, t_bam, t_bai, n_bam, n_bai, indel, tbi ->
                tuple(meta, t_bam, t_bai, n_bam, n_bai, indel, tbi)
            }

        STRELKA2(
            ch_strelka_joined.map { meta, t_bam, t_bai, n_bam, n_bai, _indel, _tbi ->
                tuple(meta, t_bam, t_bai, n_bam, n_bai)
            },
            ch_ref_fasta, ch_ref_fai,
            ch_call_regions_gz, ch_call_regions_tbi,
            ch_strelka_joined.map { _m, _tb, _ti, _nb, _ni, indel, _tbi -> indel },
            ch_strelka_joined.map { _m, _tb, _ti, _nb, _ni, _indel, tbi -> tbi }
        )
        // Strelka2 doesn't emit AF/AD natively; inject them so downstream
        // consensus/filter/vcf2maf see fields comparable to Mutect2's.
        TAG_STRELKA_AF(STRELKA2.out.vcf)
    }
    // MuSE - independent SNV caller, runs in parallel with Mutect2/Strelka2.
    if (params.use_muse_in_consensus) {
        MUSE_CALL(ch_paired, ch_ref_fasta, ch_ref_fai)
        MUSE_SUMP(MUSE_CALL.out.call_output, ch_dbsnp, ch_dbsnp_tbi)
        // compute AF from AD (which MuSE reliably has), defensively. If AD is (ref,alt), AF = alt/(ref+alt)
        TAG_MUSE_AF(MUSE_SUMP.out.vcf, ch_ref_fasta, ch_ref_fai)
    }
    // VarScan2 - deep-data-safe indel caller (pileup-based, no realigner to
    // crash on high-depth data). Runs on the paired BAMs directly, like MuSE.
    if (params.use_varscan_in_consensus) {
        // WGS: scatter by shard; WES: one task on the capture BED (unchanged behaviour).
        ch_vs_shards = (params.assay == 'wgs') ? ch_scatter_intervals : channel.of(tuple('all', []))
        VARSCAN(
            ch_paired.combine(ch_vs_shards),
            ch_ref_fasta, ch_ref_fai,
            ch_interval_bed
        )
        VARSCAN_MERGE(
            VARSCAN.out.vcf_chunk
                .map { meta, snv, indel -> tuple(meta.sample_id, meta, snv, indel) }
                .groupTuple(by: 0)
                .map { _sid, metas, snvs, indels -> tuple(metas[0], snvs.sort { v -> v.name }, indels.sort { v -> v.name }) },
            ch_ref_fai
        )
        // VarScan's AD is alt-only and FREQ is a percent string; normalize to
        // standard AF/AD so downstream consensus/filter/vcf2maf are consistent.
        TAG_VARSCAN_AF(VARSCAN_MERGE.out.vcf, ch_ref_fasta, ch_ref_fai)
    }

    ch_mutect_input = ch_paired.combine(ch_scatter_intervals)
    MUTECT2(
        ch_mutect_input,
        ch_ref_fasta, ch_ref_fai, ch_ref_dict,
        ch_germline_resource, ch_germline_resource_tbi,
        ch_panel_of_normals, ch_panel_of_normals_tbi
    )

    // GROUPING: key on scalar sample_id. Shard lists sorted by name so the gather
    // task's input hash does not depend on shard completion order (-resume).
    ch_mutect_chunks = MUTECT2.out.vcf_chunk
        .map { meta, _iname, vcf, tbi, stats -> tuple(meta.sample_id, meta, vcf, tbi, stats) }
        .groupTuple(by: 0)
        .map { _sid, metas, vcfs, tbis, stats ->
            tuple(metas[0], vcfs.sort { v -> v.name }, tbis.sort { t -> t.name }, stats.sort { s -> s.name })
        }

    // -- Orientation-bias model (FFPE) - gather per-sample f1r2, learn priors -
    if (params.orientation_bias) {
        ch_orient_input = MUTECT2.out.f1r2
            .map { meta, _iname, f -> tuple(meta.sample_id, meta, f) }
            .groupTuple(by: 0)
            .map { _sid, metas, fs -> tuple(metas[0], fs.sort { f -> f.name }) }
        LEARN_READ_ORIENTATION_MODEL(ch_orient_input)
        ch_ob_priors = LEARN_READ_ORIENTATION_MODEL.out.priors.map { meta, p -> tuple(meta.sample_id, p) }
    } else {
        ch_ob_priors = channel.empty()
    }

    // -- Contamination - pileups (tumor + normal) - CalculateContamination --
    //   Runs only if a pileup-variants resource is provided; else skipped.
    if (params.pileup_variants) {
        GETPILEUPSUMMARIES(
            ch_paired.map { meta, t_bam, t_bai, _nb, _ni -> tuple(meta, 'tumor', t_bam, t_bai) }
                .mix(ch_paired.map { meta, _tb, _ti, n_bam, n_bai -> tuple(meta, 'normal', n_bam, n_bai) }),
            ch_ref_fasta, ch_ref_fai, ch_ref_dict,
            ch_pileup_variants, ch_pileup_variants_tbi,
            ch_interval_bed
        )
        ch_pile_t = GETPILEUPSUMMARIES.out.table
            .filter { _m, label, _t -> label == 'tumor' }
            .map    { meta, _l, tbl -> tuple(meta.sample_id, meta, tbl) }
        ch_pile_n = GETPILEUPSUMMARIES.out.table
            .filter { _m, label, _t -> label == 'normal' }
            .map    { meta, _l, tbl -> tuple(meta.sample_id, tbl) }

        ch_contam_input = ch_pile_t.join(ch_pile_n, by: 0)
            .map { _sid, meta, t_tbl, n_tbl -> tuple(meta, t_tbl, n_tbl) }

        CALCULATE_CONTAMINATION(ch_contam_input)

        ch_contam = CALCULATE_CONTAMINATION.out.contamination
            .map { meta, c, seg -> tuple(meta.sample_id, c, seg) }
        ch_contam_report = CALCULATE_CONTAMINATION.out.contamination
            .map { meta, c, _seg -> tuple(meta.sample_id, c) }
    } else {
        ch_contam = channel.empty()
        ch_contam_report = channel.empty()
    }

    // -- Join chunks + (optional) contamination + (optional) ob-priors -------
    //   feature on: inner join on the optional channel; feature off: [] placeholders.
    //   c/seg/p come from one channel, so samples cannot cross-pair.
    ch_merge_keyed = ch_mutect_chunks.map { meta, v, t, s -> tuple(meta.sample_id, meta, v, t, s) }

    ch_with_contam = params.pileup_variants
        ? ch_merge_keyed.join(ch_contam, by: 0)
            .map { sid, meta, v, t, s, c, seg -> tuple(sid, meta, v, t, s, c, seg) }
        : ch_merge_keyed.map { sid, meta, v, t, s -> tuple(sid, meta, v, t, s, [], []) }

    ch_with_ob = params.orientation_bias
        ? ch_with_contam.join(ch_ob_priors, by: 0)
            .map { _sid, meta, v, t, s, c, seg, p -> tuple(meta, v, t, s, c, seg, p) }
        : ch_with_contam.map { _sid, meta, v, t, s, c, seg -> tuple(meta, v, t, s, c, seg, []) }

    MERGE_MUTECT2(
        ch_with_ob.map { m, v, t, s, _c, _seg, _p -> tuple(m, v, t, s) },
        ch_ref_fasta, ch_ref_fai, ch_ref_dict,
        ch_with_ob.map { _m, _v, _t, _s, c,  _seg, _p -> c   },
        ch_with_ob.map { _m, _v, _t, _s, _c, seg,  _p -> seg },
        ch_with_ob.map { _m, _v, _t, _s, _c, _seg, p  -> p   }
    )

// -- Somatic SNV/indel final call set ----------------------------------
    if (params.somatic_consensus) {
        // Variant-type-aware N-of-M consensus across enabled callers.
        // Mutect2 is the always-on backbone; MuSE/Strelka2 optional per param.
        // Each caller channel is keyed by sample_id; disabled or sample-missing
        // callers are padded with [] so no sample is dropped.


        ch_cons_m2 = MERGE_MUTECT2.out.vcf
            .map { meta, v, t -> tuple(meta.sample_id, meta, v, t) }

        ch_cons_strelka = (params.use_strelka2_in_consensus)
            ? TAG_STRELKA_AF.out.vcf.map { meta, s, st, i, it -> tuple(meta.sample_id, s, st, i, it) }
            : channel.empty()

        ch_cons_muse = (params.use_muse_in_consensus)
            ? TAG_MUSE_AF.out.vcf.map { meta, v, t -> tuple(meta.sample_id, v, t) }
            : channel.empty()

        ch_cons_varscan = (params.use_varscan_in_consensus)
            ? TAG_VARSCAN_AF.out.vcf.map { meta, s, st, i, it -> tuple(meta.sample_id, s, st, i, it) }
            : channel.empty()

        // Left-join onto the Mutect2 backbone (remainder:true keeps every sample,
        // null-filling absent callers), then map nulls to [] placeholders.

        // after m2 (sid,meta,v,t)  strelka(sid,ssnv,ssnvt,sind,sindt): row = (sid,meta,v,t,ssnv,ssnvt,sind,sindt) [0..7]
        ch_consensus_input = ch_cons_m2
            .join(ch_cons_strelka, by: 0, remainder: true)
            .map { row ->
                def ssnv  = row.size() > 4 ? (row[4] ?: []) : []
                def ssnvt = row.size() > 5 ? (row[5] ?: []) : []
                def sind  = row.size() > 6 ? (row[6] ?: []) : []
                def sindt = row.size() > 7 ? (row[7] ?: []) : []
                tuple(row[0], row[1], row[2], row[3], ssnv, ssnvt, sind, sindt)
                // -> (sid, meta, v, t, ssnv, ssnvt, sind, sindt)  [0..7]
            }
            //  muse(sid,musev,muset): row = (...,musev,muset) [8,9]
            .join(ch_cons_muse, by: 0, remainder: true)
            .map { row ->
                def musev = row.size() > 8 ? (row[8] ?: []) : []
                def muset = row.size() > 9 ? (row[9] ?: []) : []
                tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], musev, muset)
                // -> (sid, meta, v, t, ssnv, ssnvt, sind, sindt, musev, muset)  [0..9]
            }
            //  varscan(sid,vsnv,vsnvt,vind,vindt): row = (...,vsnv,vsnvt,vind,vindt) [10..13]
            .join(ch_cons_varscan, by: 0, remainder: true)
            .map { row ->
                def vsnv  = row.size() > 10 ? (row[10] ?: []) : []
                def vsnvt = row.size() > 11 ? (row[11] ?: []) : []
                def vind  = row.size() > 12 ? (row[12] ?: []) : []
                def vindt = row.size() > 13 ? (row[13] ?: []) : []
                tuple(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], vsnv, vsnvt, vind, vindt)
                // final 13 elements: (meta, v, t, ssnv, ssnvt, sind, sindt, musev, muset, vsnv, vsnvt, vind, vindt)
            }

        ENSEMBLE_CONSENSUS(ch_consensus_input, ch_ref_fasta, ch_ref_fai)
        ch_final_snv_vcf = ENSEMBLE_CONSENSUS.out.vcf

    } else {
        ch_final_snv_vcf = MERGE_MUTECT2.out.vcf
    }

    // Optional variant-level soft-filtering (VAF/depth/tumor-in-normal).
    if (params.somatic_variant_filter) {
        FILTER_SOMATIC_VCF(ch_final_snv_vcf)
        ch_final_snv_vcf = FILTER_SOMATIC_VCF.out.vcf
    }

    if (params.mutational_signatures) {
        MUTATIONAL_SIGNATURES(ch_final_snv_vcf)
    }

    // CNV - matched normal BAM wired explicitly
    ch_cnvkit_input = ch_paired.map { meta, t_bam, t_bai, n_bam, n_bai ->
        tuple(meta, t_bam, t_bai, n_bam, n_bai)
    }
    CNVKIT(
        ch_cnvkit_input,
        ch_ref_fasta,
        ch_interval_bed,
        ch_cnvkit_pon
    )
    // LOH + allele-specific CN + purity/ploidy (SEQUENZA) - runs only if
    // sequenza_gc_wiggle are provided; else skipped.
    ch_sequenza = channel.empty()
    if (params.sequenza_gc_wiggle) {
        SEQUENZA(
            ch_paired,
            ch_ref_fasta, ch_ref_fai,
            ch_sequenza_gc_wiggle
        )
        ch_sequenza = SEQUENZA.out.sequenza
    }
    // CNVkit call: purity-rescaled (-m clonal) with Sequenza cellularity when available.
    // ch_sequenza is empty when Sequenza did not run.
    ch_seq_confints = ch_sequenza.map { row -> tuple(row[0].sample_id, row[2]) }
    CNVKIT_CALL(
        CNVKIT.out.cnv
            .map { meta, cns, cnr -> tuple(meta.sample_id, meta, cns, cnr) }
            .join(ch_seq_confints, by: 0, remainder: true)
            .filter { row -> row[1] instanceof Map }   // drop Sequenza-only rows (CNVkit missing)
            .map { row -> tuple(row[1], row[2], (row.size() > 4 && row[4]) ? row[4] : []) }
    )
    // (meta, cns, cnr, call.cns) - same shape the report consumed before
    ch_cnv_out = CNVKIT.out.cnv
        .map { meta, cns, cnr -> tuple(meta.sample_id, meta, cns, cnr) }
        .join(CNVKIT_CALL.out.call_cns.map { meta, call -> tuple(meta.sample_id, call) }, by: 0)
        .map { _sid, meta, cns, cnr, call -> tuple(meta, cns, cnr, call) }

    ch_gene_loh = channel.empty()
    if (params.sequenza_gc_wiggle && params.gene_bed) {
        SEQUENZA_GENE_LOH( SEQUENZA.out.sequenza.map { row -> tuple(row[0], row[1], row[2]) }, ch_gene_bed,
                           toOptionalFileChannel(params.arm_bed) )
        ch_gene_loh = SEQUENZA_GENE_LOH.out.gene_loh
    } else if (!params.sequenza_gc_wiggle) {
        log.info "SEQUENZA (LOH/purity) skipped (params.sequenza_gc_wiggle not set)."
    } else {
        log.info "SEQUENZA gene-level LOH skipped (params.gene_bed not set); segment-level Sequenza still runs."
    }

    // MSI - needs the pre-scanned site list; without it msisensor-pro would fail on '-d '
    // and stop the whole run, so the step is skipped instead (main.nf warns).
    ch_msi        = channel.empty()
    ch_msi_report = channel.empty()
    if (params.msisensor_list) {
        MSISENSOR_PRO(
            ch_paired.map { meta, t_bam, t_bai, _n, _ni -> tuple(meta, t_bam, t_bai) },
            ch_paired.map { _m, _t, _ti, n_bam, _ni     -> n_bam },
            ch_paired.map { _m, _t, _ti, _n, n_bai      -> n_bai },
            ch_msisensor_list,
            ch_interval_bed
        )
        ch_msi = MSISENSOR_PRO.out.msi

        if (params.gene_bed) {
            ch_msi_somatic = MSISENSOR_PRO.out.msi
                .map    { meta, msi_files ->
                    def fl = (msi_files instanceof List) ? msi_files : [msi_files]
                    // paired mode lists the unstable sites in <prefix>_unstable (MSIsensor-pro 1.x;
                    // <prefix>_somatic in older releases)
                    tuple(meta, fl.find { f -> f.name.endsWith('_unstable') } ?: fl.find { f -> f.name.endsWith('_somatic') })
                }
                .filter { _m, f -> f != null }                 // no per-site table (tumor-only) - skip
            MSISENSOR_GENE_MAP(ch_msi_somatic, ch_gene_bed)
        }

        // The patient report shows the CLASSIFIED status (cutoff + minimum sites) when available.
        ch_msi_report = MSISENSOR_PRO.out.msi
        if (params.classify_msi) {
            CLASSIFY_MSI(MSISENSOR_PRO.out.msi)
            ch_msi_report = CLASSIFY_MSI.out.status
        }
    }

    emit:
    snv_vcf  = ch_final_snv_vcf
    sv_vcf   = MANTA.out.sv_vcf
    cnv      = ch_cnv_out
    loh      = ch_sequenza
    gene_loh = ch_gene_loh
    msi      = ch_msi
    msi_report    = ch_msi_report
    contamination = ch_contam_report
}

// ============================================================
//  SUB-WORKFLOW: SOMATIC_REPORT_BLOCK
//  Consolidated per-patient multi-variant-type report.
// ============================================================
workflow SOMATIC_REPORT_BLOCK {
    take:
    ch_classification   // CLASSIFY_VARIANTS.out.classification
    ch_cnv              // SOMATIC_BLOCK.out.cnv
    ch_loh              // SOMATIC_BLOCK.out.loh
    ch_sv               // SOMATIC_BLOCK.out.sv_vcf
    ch_msi              // SOMATIC_BLOCK.out.msi
    ch_tmb              // COMPUTE_TMB.out.tmb
    ch_annotsv_dir      // AnnotSV DB dir (or [] when disabled)
    ch_contam           // SOMATIC_BLOCK.out.contamination (sample_id, table) or empty

    main:
    // Key = tumour sample_id. Report name = patient_id for the classic <patient>_T tumour
    // (unchanged file names), the sample_id when the samplesheet 'sample' column is used.
    ch_rep_class = ch_classification
        .map { meta, tsv ->
            def report_id = (meta.sample_id == "${meta.patient_id}_T") ? meta.patient_id : meta.sample_id
            tuple(meta.sample_id, report_id, (meta.sex ?: 'NA').toString(), tsv) }
    ch_rep_cnv = ch_cnv
        .map { meta, _mk, _mkr, callcns -> tuple(meta.sample_id, callcns) }
    ch_rep_loh = ch_loh
        .map { row -> tuple(row[0].sample_id, row[1]) }

    ch_rep_sv = channel.empty()
    if (params.annotsv_dir) {
        ANNOTATE_SV_SOMATIC(
            ch_sv.map { meta, vcf, tbi -> tuple(meta, 'manta_somatic', vcf, tbi) },
            ch_annotsv_dir
        )
        ch_rep_sv = ANNOTATE_SV_SOMATIC.out.annotated_sv
            .map { meta, _tag, tsv -> tuple(meta.sample_id, tsv) }
    }

    ch_rep_msi = ch_msi
        .map { meta, msi ->
            // MSISENSOR_PRO emits a glob (_msi, _msi_dis, ...); pick the main summary file (basename == "<sample>_msi").
            def files = (msi instanceof List) ? msi : [msi]
            def main  = files.find { f -> f.name == "${meta.sample_id}_msi_status.txt" } ?:
                        files.find { f -> f.name == "${meta.sample_id}_msi" } ?: files[0]
            tuple(meta.sample_id, main)
        }
    ch_rep_tmb = ch_tmb
        .map { meta, tmb -> tuple(meta.sample_id, tmb) }

    ch_report_in = ch_rep_class
        .join(ch_rep_cnv, by: 0, remainder: true)
        .join(ch_rep_loh, by: 0, remainder: true)
        .join(ch_rep_sv,  by: 0, remainder: true)
        .join(ch_rep_msi, by: 0, remainder: true)
        .join(ch_rep_tmb, by: 0, remainder: true)
        .join(ch_contam, by: 0, remainder: true)
        .map { row ->
            // key row: (sid, report_id, sex, cls, cnv, loh, sv, msi, tmb, contamination)
            // emitted:  (report_id, sid, sex, ...) - process input order unchanged
            def sid = row[0]
            def pid = row.size() > 1 ? (row[1] ?: sid) : sid
            def sex = row.size() > 2 ? (row[2] ?: 'NA') : 'NA'
            def cls = row.size() > 3 ? (row[3] ?: []) : []
            def cnv = row.size() > 4 ? (row[4] ?: []) : []
            def loh = row.size() > 5 ? (row[5] ?: []) : []
            def sv  = row.size() > 6 ? (row[6] ?: []) : []
            def msi = row.size() > 7 ? (row[7] ?: []) : []
            def tmb = row.size() > 8 ? (row[8] ?: []) : []
            def con = row.size() > 9 ? (row[9] ?: []) : []
            tuple(pid, sid, sex, cls, cnv, loh, sv, msi, tmb, con)
        }


    CONSOLIDATE_PATIENT_REPORT(ch_report_in)

    emit:
    report = CONSOLIDATE_PATIENT_REPORT.out.report
}


// ============================================================
//  SUB-WORKFLOW: GERMLINE_REPORT_BLOCK   (Phase 1)
//  Annotate germline SV (Manta) + CNV (gCNV) with AnnotSV, then
//  consolidate SNV/indel + SV + CNV + ROH into one per-sample report.
//
//  Keyed on sample_id throughout (germline is per-sample). Every input
//  channel is optional: SV/CNV/ROH may be empty, and the final join uses
//  remainder:true + [] padding so a sample with only SNVs still reports.
// ============================================================
workflow GERMLINE_REPORT_BLOCK {
    take:
    ch_snv_vcf      // GERMLINE_BLOCK.out.snv_vcf  (meta, vcf, tbi) - VEP-annotated upstream
    ch_sv_vcf       // GERMLINE_BLOCK.out.sv_vcf   (meta, vcf, tbi)
    ch_cnv_vcf      // GERMLINE_BLOCK.out.cnv_vcf  (meta, iv, ivtbi, seg, segtbi) or empty
    ch_roh          // GERMLINE_BLOCK.out.roh      (meta, txt) or empty
    ch_annotsv_dir
    ch_cgc_tsv      // CGC (GERMLINE column) or []
    ch_report_genes // params.germline_report_genes or []

    main:
    // SV/CNV annotation needs the AnnotSV DB; without it those report sections are
    // empty placeholders and the SNV / LoF / ROH sections are still produced.
    ch_r_sv  = channel.empty()
    ch_r_cnv = channel.empty()
    if (params.annotsv_dir) {
        // -- Annotate SV (Manta diploidSV) -------------------------
        ch_sv_in = ch_sv_vcf.map { meta, vcf, tbi -> tuple(meta, 'manta', vcf, tbi) }
        ANNOTATE_SV_MANTA(ch_sv_in, ch_annotsv_dir)
        ch_r_sv  = ANNOTATE_SV_MANTA.out.annotated_sv.map { meta, _s, tsv -> tuple(meta.sample_id, tsv) }

        // -- Annotate CNV (gCNV segments - the 2nd VCF pair) -------
        ch_cnv_in = ch_cnv_vcf.map { meta, _iv, _ivt, seg, segt -> tuple(meta, 'gcnv', seg, segt) }
        ANNOTATE_SV_GCNV(ch_cnv_in, ch_annotsv_dir)
        ch_r_cnv = ANNOTATE_SV_GCNV.out.annotated_sv.map  { meta, _s, tsv -> tuple(meta.sample_id, tsv) }
    } else {
        log.info "Germline report: SV/CNV sections empty (params.annotsv_dir not set)."
    }

    // -- Assemble per-sample report inputs (key = sample_id) ---
    ch_r_snv = ch_snv_vcf.map { meta, vcf, tbi -> tuple(meta.sample_id, vcf, tbi) }
    ch_r_roh = ch_roh.map { meta, txt -> tuple(meta.sample_id, txt) }

    // SNV backbone (always present); left-join optional SV/CNV/ROH.
    ch_report_in = ch_r_snv
        .join(ch_r_sv,  by: 0, remainder: true)
        .map { row ->
            def sid = row[0]
            def vcf = row[1]; def tbi = row[2]
            def sv  = row.size() > 3 ? (row[3] ?: []) : []
            tuple(sid, vcf, tbi, sv)
        }
        .join(ch_r_cnv, by: 0, remainder: true)
        .map { row ->
            def cnv = row.size() > 4 ? (row[4] ?: []) : []
            tuple(row[0], row[1], row[2], row[3], cnv)
        }
        .join(ch_r_roh, by: 0, remainder: true)
        .map { row ->
            def roh = row.size() > 5 ? (row[5] ?: []) : []
            tuple(row[0], row[1], row[2], row[3], row[4], roh)
        }

    CONSOLIDATE_GERMLINE_REPORT(ch_report_in, ch_cgc_tsv, ch_report_genes)

    emit:
    report = CONSOLIDATE_GERMLINE_REPORT.out.report
}

// ============================================================
//  DNA WORKFLOW  -  Top-level orchestration
// ============================================================
workflow DNA_WORKFLOW {

    take:
    samplesheet

    main:
    // Stage ordering comes from the shared map in utils_v3 (same one main's
    // validation uses) - no local duplicate that could silently desync.
    int start_rank = stageRankMap()[params.start_at]
    // int stop_rank  = stageRankMap()[params.stop_at]

    // -- Reference channels (value - broadcast to all tasks) ---
    ch_ref_fasta = toFileChannel(params.ref_fasta)
    ch_ref_fai   = toFileChannel(params.ref_fai)
    ch_ref_dict  = toFileChannel(params.ref_dict)

    // Stage ONLY the bwa-mem2 index sidecar files - never the .fa/.fai/.dict.
    // Globbing the whole dir would re-stage the FASTA (already passed as
    // ch_ref_fasta) into the SAME work dir as the index, causing a Nextflow
    // "input file name collision". The sidecar suffixes are tool-defined and
    // fixed; the .fa/.fai/.dict are naturally excluded by this allowlist.
    // files() always returns a list (file() returns a bare Path when exactly one entry
    // matches, and a Path iterates over its name components).
    def bwamem2_suffixes = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
    def ref_index_files  = params.ref_index_dir
        ? files("${params.ref_index_dir}/*").findAll { f -> bwamem2_suffixes.any { suf -> f.name.endsWith(suf) } }
        : []
    // bwa-mem2 runs in mapping (all modes) and in umi_consensus re-alignment (markdup stage).
    if (stageRuns('mapping') || (stageRuns('markdup') && params.dedup_mode == 'umi_consensus')) {
        def missing_idx = bwamem2_suffixes.findAll { suf -> !ref_index_files.any { f -> f.name.endsWith(suf) } }
        if (missing_idx)
            error "ERROR: bwa-mem2 index incomplete in params.ref_index_dir='${params.ref_index_dir}': " +
                  "no file ending in ${missing_idx.join(', ')}"
    }
    ch_ref_index = channel.value(ref_index_files)

    // dbSNP is optional at channel level; main.nf requires it where it is consumed
    // (HaplotypeCaller, BQSR, MuSE sump).
    ch_dbsnp            = toOptionalFileChannel(params.dbsnp)
    ch_dbsnp_tbi        = params.dbsnp ? toFileChannel("${params.dbsnp}.tbi") : channel.value([])
    ch_mills            = params.mills        ? toFileChannel(params.mills)               : channel.value([])
    ch_mills_tbi        = params.mills        ? toFileChannel("${params.mills}.tbi")      : channel.value([])
    ch_known_indels     = params.known_indels ? toFileChannel(params.known_indels)        : channel.value([])
    ch_known_indels_tbi = params.known_indels ? toFileChannel("${params.known_indels}.tbi") : channel.value([])

    ch_interval_bed  = toOptionalFileChannel(params.interval_bed)
    //ch_interval_list = toOptionalFileChannel(params.interval_list)
    ch_bait_interval_list   = toOptionalFileChannel(params.bait_interval_list)
    ch_target_interval_list = toOptionalFileChannel(params.target_interval_list)

    // Interval shards only: any other file in the folder (checksums, logs, README) would
    // otherwise become an extra Mutect2/HaplotypeCaller/VarScan shard and fail it.
    ch_scatter_intervals = channel.fromPath("${params.scatter_intervals}/*.{interval_list,bed}", checkIfExists: true)
        .flatten().map { f -> tuple(f.baseName, f) }

    ch_germline_resource     = toOptionalFileChannel(params.germline_resource)
    ch_germline_resource_tbi = params.germline_resource ? toFileChannel("${params.germline_resource}.tbi") : channel.value([])
    ch_panel_of_normals      = toOptionalFileChannel(params.panel_of_normals)
    ch_panel_of_normals_tbi  = params.panel_of_normals ? toFileChannel("${params.panel_of_normals}.tbi") : channel.value([])

    ch_verifybamid_db = params.verifybamid_db ? channel.value(file(params.verifybamid_db, checkIfExists: true)) : channel.value([])

    ch_hapmap         = toOptionalFileChannel(params.hapmap)
    ch_hapmap_tbi     = params.hapmap         ? toFileChannel("${params.hapmap}.tbi")         : channel.value([])
    ch_omni           = toOptionalFileChannel(params.omni)
    ch_omni_tbi       = params.omni           ? toFileChannel("${params.omni}.tbi")           : channel.value([])
    ch_thousandg_snps = toOptionalFileChannel(params.thousandg_snps)
    ch_thousandg_tbi  = params.thousandg_snps ? toFileChannel("${params.thousandg_snps}.tbi") : channel.value([])

    ch_preprocessed_intervals = toOptionalFileChannel(params.preprocessed_intervals)
    ch_contig_ploidy_priors   = toOptionalFileChannel(params.contig_ploidy_priors)
    ch_gcnv_scatter_intervals = params.gcnv_scatter_intervals ? channel.fromPath("${params.gcnv_scatter_intervals}/*.{interval_list,bed}", checkIfExists: true)
        .flatten().map { f -> tuple(f.baseName, f) } : channel.empty()
    ch_gcnv_annotated_intervals = toOptionalFileChannel(params.gcnv_annotated_intervals)

    ch_germ_scatter = params.germline_scatter_intervals
                ? channel.fromPath("${params.germline_scatter_intervals}/*.interval_list", checkIfExists: true)
                    .flatten().map { f -> tuple(f.baseName, f) }
                : ch_scatter_intervals

    ch_cnvkit_pon     = toOptionalFileChannel(params.cnvkit_pon)
    ch_msisensor_list = toOptionalFileChannel(params.msisensor_list)
    ch_pileup_variants     = toOptionalFileChannel(params.pileup_variants)
    ch_pileup_variants_tbi = params.pileup_variants ? toFileChannel("${params.pileup_variants}.tbi") : channel.value([])

    ch_sequenza_gc_wiggle  = toOptionalFileChannel(params.sequenza_gc_wiggle)
    ch_gene_bed            = toOptionalFileChannel(params.gene_bed)

    ch_somalier_sites      = toOptionalFileChannel(params.somalier_sites)
    ch_evidence_bed    = toOptionalFileChannel(params.evidence_bed)

    ch_vep_cache_dir   = params.vep_cache_dir   ? channel.value(file(params.vep_cache_dir,   checkIfExists: true)) : channel.value([])
    ch_vep_plugins_dir = params.vep_plugins_dir ? channel.value(file(params.vep_plugins_dir, checkIfExists: true)) : channel.value([])

    ch_clinvar         = toOptionalFileChannel(params.clinvar_vcf)
    ch_clinvar_tbi     = params.clinvar_vcf ? toFileChannel("${params.clinvar_vcf}.tbi") : channel.value([])
    ch_cosmic          = toOptionalFileChannel(params.cosmic_vcf)
    ch_cosmic_tbi      = params.cosmic_vcf ? toFileChannel("${params.cosmic_vcf}.tbi") : channel.value([])
    ch_gnomad_faf      = toOptionalFileChannel(params.gnomad_faf_vcf)
    ch_gnomad_faf_tbi  = params.gnomad_faf_vcf ? toFileChannel("${params.gnomad_faf_vcf}.tbi") : channel.value([])

    // -- Germline Phase 1/2: ROH AF-file (optional) + AnnotSV DB dir --
    ch_roh_af_file     = toOptionalFileChannel(params.roh_af_file)
    ch_roh_af_file_tbi = params.roh_af_file ? toFileChannel("${params.roh_af_file}.tbi") : channel.value([])
    ch_annotsv_dir     = params.annotsv_dir
        ? channel.value(file(params.annotsv_dir, checkIfExists: true)) : channel.value([])
    // Pre-bgzipped + tabixed call-regions BED (prepared ONCE outside the
    // pipeline; see references README). Pass params.call_regions_gz pointing at
    // the .bed.gz; the .tbi is derived as <gz>.tbi. Reused across all runs.
    if (params.assay == 'wes' && params.call_regions_gz) {
        ch_call_regions_gz  = toFileChannel(params.call_regions_gz)
        ch_call_regions_tbi = toFileChannel("${params.call_regions_gz}.tbi")
    } else {
        ch_call_regions_gz  = channel.value([])
        ch_call_regions_tbi = channel.value([])
    }

    ch_active_bam = channel.empty()
    ch_active_vcf = channel.empty()

    // -- Reference preflight (build / contig naming), in parallel with mapping --
    if (params.preflight_references) {
        def pf_vcfs = ['dbsnp', 'mills', 'known_indels', 'germline_resource', 'panel_of_normals',
                       'pileup_variants', 'clinvar_vcf', 'cosmic_vcf', 'gnomad_faf_vcf', 'somalier_sites']
                      .findAll { k -> params[k] }.collect { k -> file(params[k], checkIfExists: true) }
        def pf_beds = (params.assay == 'wes' ? ['interval_bed', 'call_regions_gz'] : [])
                      .findAll { k -> params[k] }.collect { k -> file(params[k], checkIfExists: true) }
        PREFLIGHT_REFERENCES(ch_ref_dict, channel.value(pf_vcfs), channel.value(pf_beds))
    }

    // -- Manifest ----------------------------------------------
    GENERATE_META_INFO(samplesheet)
    ch_manifest = GENERATE_META_INFO.out.reads

    // All rows (lanes/chunks) of one sample must agree on sample-level fields; otherwise
    // grouping splits the sample and its outputs are produced (and published) twice.
    ch_manifest
        .map { meta, _reads -> tuple(meta.sample_id, [patient: meta.patient_id, sex: meta.sex, batch: meta.batch, type: meta.type]) }
        .unique()
        .groupTuple()
        .map { sid, vals -> if (vals.size() > 1) error "Samplesheet rows of ${sid} disagree on patient/sex/batch/sample_type: ${vals} - every row (lane/chunk) of one sample needs the same values (two different samples need different 'sample' labels)"; sid }

    // Read-group IDs (<sample>_L<lane>_C<chunk>) must be unique: two rows with the same lane
    // and chunk would give two read groups one ID, and duplicate marking/BQSR would merge them.
    ch_manifest
        .map { meta, _reads -> tuple(meta.uid, 1) }
        .groupTuple()
        .map { uid, n -> if (n.size() > 1) error "Samplesheet has ${n.size()} rows for read group ${uid} (same sample, lane_num and chunk); give each FASTQ pair its own lane_num/chunk"; uid }

    def umi_mode = (params.dedup_mode in ['umi_consensus', 'umi_markdup'])

    // Re-entry works per SAMPLE, not per samplesheet row (chunk): one item per
    // sample_id, chunk-level meta (uid, rg) stripped - same meta as the mapping path.
    ch_samples = ch_manifest
        .map { meta, _reads -> stripChunkMeta(meta) }
        .unique { m -> m.sample_id }

    // Somatic pairing needs exactly one normal per patient (tumours may be several).
    ch_samples
        .filter { m -> params.analysis == 'somatic' && m.type == 'DN' }
        .map    { m -> tuple(m.patient_id, m.sample_id) }
        .groupTuple()
        .map    { pid, ids -> if (ids.size() > 1) error "Patient ${pid} has ${ids.size()} normal (DN) samples ${ids}; tumour-normal pairing needs exactly one"; pid }

    // VCF re-entry (annotation/classification): only samples that have a
    // callset in this analysis - tumors (somatic) or normals (germline).
    // Somatic tumours carry normal_id (the matched normal's SM) for VCF2MAF on cold start.
    def vcf_type = (params.analysis == 'germline') ? 'DN' : 'DT'
    // A tumour listed without its normal keeps the classic '<patient>_N' name (the normal
    // column name in its VCF) instead of being dropped from re-entry.
    ch_normal_ids = ch_samples
        .filter { m -> m.type == 'DN' }
        .map    { m -> [m.patient_id, m.sample_id] }
        .toList()
        .map    { pairs -> pairs.collectEntries { pr -> [(pr[0]): pr[1]] } }
    ch_vcf_samples = (vcf_type == 'DT')
        ? ch_samples.filter { m -> m.type == 'DT' }
            .combine(ch_normal_ids)
            .map { m, nmap -> m + [normal_id: nmap[m.patient_id] ?: "${m.patient_id}_N".toString()] }
        : ch_samples.filter { m -> m.type == vcf_type }

    if (params.write_manifest) {
        // One row per sample (input file is the first read / BAM), then consolidate.
        ch_manifest_input = ch_manifest.map { meta, reads ->
            // 'reads' may be a list (paired FASTQ) or a single file; take the first.
            def f = (reads instanceof List) ? reads[0] : reads
            tuple(meta, f)
        }
        MANIFEST_ROW(ch_manifest_input)
        WRITE_MANIFEST(MANIFEST_ROW.out.row.collect(sort: true))
    }

    // -- Entry-point routing (keyed on the stage name = params.start_at) ------
    // Stage-2 BAM name follows dedup_mode (umi_consensus writes <sid>_consensus.bam);
    // variant_calling re-entry reads <sid>_bqsr.bam only when bqsr=true.
    def dedup_sfx = (params.dedup_mode == 'umi_consensus') ? 'consensus' : 'markdup'
    if (params.start_at == 'mapping') {
        ch_alignment_inputs = ch_manifest
    }
    else if (params.start_at == 'markdup') {
        // chunk BAM names: <uid>.fixmate.bam (markdup) | <uid>.umi.mapped.bam (umi_*),
        // uid = <sample_id>_L<lane>_C<chunk>
        def chunk_suffix = umi_mode ? 'umi.mapped.bam' : 'fixmate.bam'
        ch_active_bam = ch_samples.map { meta ->
            def chunks = files("${params.input_bam_dir}/${meta.sample_id}_L*_C*.${chunk_suffix}")
            if (!chunks) error "start_at='markdup': no *.${chunk_suffix} chunks for ${meta.sample_id} in ${params.input_bam_dir}"
            tuple(meta, chunks)
        }
    }
    else if (params.start_at in ['metrics', 'bqsr']) {
        ch_active_bam = ch_samples.map { meta ->
            def cram = file("${params.input_bam_dir}/${meta.sample_id}.cram")
            def bam  = file("${params.input_bam_dir}/${meta.sample_id}_${dedup_sfx}.bam")
            def final_bam = cram.exists() ? cram : bam
            def final_idx = cram.exists()
                ? file("${params.input_bam_dir}/${meta.sample_id}.cram.crai",       checkIfExists: true)
                : file("${params.input_bam_dir}/${meta.sample_id}_${dedup_sfx}.bai", checkIfExists: true)
            return tuple(meta, final_bam, final_idx)
        }
    }
    else if (params.start_at == 'variant_calling') {
        def vc_sfx = params.bqsr ? 'bqsr' : dedup_sfx
        ch_active_bam = ch_samples.map { meta ->
            tuple(meta,
                file("${params.input_bam_dir}/${meta.sample_id}_${vc_sfx}.bam", checkIfExists: true),
                file("${params.input_bam_dir}/${meta.sample_id}_${vc_sfx}.bai", checkIfExists: true))
        }
    }
    else if (params.start_at == 'annotation') {
        // <sample>.vcf.gz, or the final call set as this pipeline publishes it in
        // 05_variant_calling (<sample>_filtered.vcf.gz; germline_filter none: _genotyped)
        ch_active_vcf = ch_vcf_samples.map { meta ->
            def names = ["${meta.sample_id}.vcf.gz", "${meta.sample_id}_filtered.vcf.gz", "${meta.sample_id}_genotyped.vcf.gz"]
            def hit   = names.collect { n -> file("${params.input_vcf_dir}/${n}") }.find { f -> f.exists() }
            if (!hit) error "start_at='annotation': none of ${names.join(', ')} in ${params.input_vcf_dir}"
            tuple(meta, hit, file("${hit}.tbi", checkIfExists: true))
        }
    }
    else if (params.start_at == 'classification') {
        // Cold-start at classification: load the annotated VCF (like stage 6 loads
        // its VCF). The consolidated report does NOT cold-start - it needs upstream
        // emit channels, so use a full -resume run to regenerate it.
        ch_active_vcf = ch_vcf_samples.map { meta ->
            tuple(meta,
                file("${params.input_annot_dir}/${meta.sample_id}_vep.vcf.gz",     checkIfExists: true),
                file("${params.input_annot_dir}/${meta.sample_id}_vep.vcf.gz.tbi", checkIfExists: true))
        }
    }

    // -- STAGE 1 mapping: alignment (UMI extraction in the umi_* modes)
    //   dedup_mode == 'markdup'     - fastp(+5bp front-trim) | bwa-mem2 | fixmate
    //   dedup_mode == 'umi_*'       - fgbio FastqToBam (3M2S+T) - bwa-mem2 - ZipperBams

    if (stageRuns('mapping')) {
        if (umi_mode) {
            ALIGN_DNA_UMI(ch_alignment_inputs, ch_ref_fasta, ch_ref_fai, ch_ref_dict, ch_ref_index)
            ch_chunk_bams = ALIGN_DNA_UMI.out.chunk_bams
        } else {
            PREPROCESS_AND_ALIGN_DNA(ch_alignment_inputs, ch_ref_fasta, ch_ref_index)
            PREPROCESS_AND_ALIGN_DNA.out.json_report
                .collect(sort: true)
                .map { files -> tuple("01_trimming", files) }
                | MULTIQC_TRIMMING
            ch_chunk_bams = PREPROCESS_AND_ALIGN_DNA.out.chunk_bams
        }

        // GROUPING: key on scalar sample_id, carry meta as payload,
        // keep first meta (all sample-level fields identical per sample).
        // Chunk order sorted by name: groupTuple emits in task-completion order, which
        // would change the merge task's input hash between runs and defeat -resume.
        ch_active_bam = ch_chunk_bams
            .map { meta, bam -> tuple(meta.sample_id, stripChunkMeta(meta), bam) }
            .groupTuple(by: 0)
            .map { _sid, metas, bams -> tuple(metas[0], bams.sort { b -> b.name }) }
    }

    // -- STAGE 2 markdup: duplicate marking or UMI consensus, CRAM --
    if (stageRuns('markdup')) {
        if (params.dedup_mode == 'umi_consensus') {
            MERGE_AND_CONSENSUS(ch_active_bam, ch_ref_fasta, ch_ref_fai, ch_ref_dict, ch_ref_index)
            ch_dedup_bam     = MERGE_AND_CONSENSUS.out.bam_bai_pair
            ch_dedup_metrics = MERGE_AND_CONSENSUS.out.metrics
            ch_dedup_label   = "02_consensus"
        } else if (params.dedup_mode == 'umi_markdup') {
            MERGE_AND_MARKDUP_UMI(ch_active_bam, ch_ref_fasta, ch_ref_fai, ch_ref_dict)
            ch_dedup_bam     = MERGE_AND_MARKDUP_UMI.out.bam_bai_pair
            ch_dedup_metrics = MERGE_AND_MARKDUP_UMI.out.metrics
            ch_dedup_label   = "02_markdup"
        } else {
            MERGE_AND_MARK_DUPLICATES(ch_active_bam)
            ch_dedup_bam     = MERGE_AND_MARK_DUPLICATES.out.bam_bai_pair
            ch_dedup_metrics = MERGE_AND_MARK_DUPLICATES.out.metrics
            ch_dedup_label   = "02_markdup"
        }

        ch_dedup_metrics
            .collect(sort: true)
            .map { files -> tuple(ch_dedup_label, files) }
            | MULTIQC_MARKDUP

        // CRAM conversion is optional (storage-only). Default off: keep BAM.
        if (params.make_cram) {
            CRAM_CONVERT_DEDUP('dedup', ch_dedup_bam, ch_ref_fasta, ch_ref_fai)
            log.info "CRAM archive written (dedup); BAM retained as active channel."
        }
        ch_active_bam = ch_dedup_bam
    }

    // -- STAGE 3 metrics: QC ----------------------------------------
    ch_depth_quantized = channel.empty()
    if (stageRuns('metrics')) {
        QC_BLOCK(
            ch_active_bam,
            ch_ref_fasta, ch_ref_fai, ch_ref_dict,
            ch_interval_bed, ch_bait_interval_list, ch_target_interval_list, ch_verifybamid_db,
            toOptionalFileChannel(params.coverage_report_bed)
        )
        ch_depth_quantized = QC_BLOCK.out.depth_quantized
        if (params.somalier_sites) {
            SOMALIER_EXTRACT(ch_active_bam, ch_ref_fasta, ch_ref_fai, ch_somalier_sites)
            SOMALIER_RELATE(SOMALIER_EXTRACT.out.counts.collect(sort: true))

            CHECK_SAMPLE_SWAPS(SOMALIER_RELATE.out.pairs)

            ch_ancestry = channel.empty()
            if (params.somalier_ancestry_labels && params.somalier_1kg_dir) {
                SOMALIER_ANCESTRY(
                    SOMALIER_EXTRACT.out.counts_meta.filter { m, _f -> m.type == 'DN' }.map { _m, f -> f }.collect(sort: true),
                    toFileChannel(params.somalier_ancestry_labels),
                    toFileChannel(params.somalier_1kg_dir)
                )
                ch_ancestry = SOMALIER_ANCESTRY.out.ancestry
            }

            SOMALIER_RELATE.out.pairs
                .mix(SOMALIER_RELATE.out.samples)
                .mix(ch_ancestry)
                .collect(sort: true)
                .map { files -> tuple("03_somalier", files) }
                | MULTIQC_SOMALIER
        } else {
            log.info "Sample-swap QC skipped (params.somalier_sites not set)."
        }
    }

    // -- STAGE 4 bqsr: base recalibration (bqsr: true) --------------
    if (stageRuns('bqsr')) {
        if (params.bqsr) {
            RUN_BQSR(
                ch_active_bam,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_dbsnp, ch_dbsnp_tbi,
                ch_mills, ch_mills_tbi,
                ch_known_indels, ch_known_indels_tbi,
                ch_interval_bed,
                ch_scatter_intervals
            )
            ch_active_bam = RUN_BQSR.out.bam_bai_pair
            if (params.make_cram) {
                CRAM_CONVERT_BQSR('bqsr', ch_active_bam, ch_ref_fasta, ch_ref_fai)
                log.info "CRAM archive written (bqsr); BAM retained as active channel."
            }
            // ch_active_bam stays the BQSR BAM (no CRAM reassignment)
        }
    }

    // -- STAGE 5 variant_calling ------------------------------------
    // Somatic runs also call the matched normals' germline (FH, SMARCB1, VHL,
    // FLCN, BAP1, TSC1/2, SDHx, HBB ...) unless germline_on_normals=false.
    def germline_on_normals = (params.analysis == 'somatic' && params.germline_on_normals)
    ch_germline_vcf = channel.empty()

    if (stageRuns('variant_calling')) {
        if (params.evidence_bed) {
            EVIDENCE_BAM(ch_active_bam, ch_evidence_bed, ch_ref_fasta, ch_ref_fai)
        }

        if (params.analysis == 'germline' || germline_on_normals) {

            GERMLINE_BLOCK(
                ch_active_bam,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_dbsnp, ch_dbsnp_tbi,
                ch_interval_bed,
                //ch_interval_list,
                ch_call_regions_gz, ch_call_regions_tbi,
                ch_hapmap, ch_hapmap_tbi,
                ch_omni, ch_omni_tbi,
                ch_thousandg_snps, ch_thousandg_tbi,
                ch_mills, ch_mills_tbi,
                ch_preprocessed_intervals,
                ch_contig_ploidy_priors,
                ch_gcnv_scatter_intervals,
                ch_gcnv_annotated_intervals,
                ch_germ_scatter,
                ch_roh_af_file,
                ch_roh_af_file_tbi
            )
        }

        if (params.analysis == 'germline') {
            ch_active_vcf = GERMLINE_BLOCK.out.snv_vcf

        } else if (params.analysis == 'somatic') {
            SOMATIC_BLOCK(
                ch_active_bam,
                ch_ref_fasta, ch_ref_fai, ch_ref_dict,
                ch_germline_resource, ch_germline_resource_tbi,
                ch_panel_of_normals, ch_panel_of_normals_tbi,
                ch_interval_bed,
                ch_call_regions_gz, ch_call_regions_tbi,
                ch_scatter_intervals,
                ch_cnvkit_pon, ch_msisensor_list,
                ch_pileup_variants, ch_pileup_variants_tbi,
                ch_sequenza_gc_wiggle,
                ch_gene_bed,
                ch_dbsnp, ch_dbsnp_tbi

            )
            ch_active_vcf = SOMATIC_BLOCK.out.snv_vcf
            if (germline_on_normals) ch_germline_vcf = GERMLINE_BLOCK.out.snv_vcf

        } else if (params.analysis == 'tumor_only') {
            error "analysis = tumor_only is not implemented; use analysis = somatic with a matched normal"
        } else {
            error "ERROR: params.analysis must be 'germline', 'somatic', or 'tumor_only'. Got: ${params.analysis}"
        }
    }

    // -- STAGE 6 annotation -----------------------------------------
    if (stageRuns('annotation')) {

        ANNOTATE_VCF(
            ch_active_vcf,
            ch_vep_cache_dir,
            ch_ref_fasta, ch_ref_fai,
            ch_vep_plugins_dir,
            ch_clinvar, ch_clinvar_tbi,
            ch_cosmic, ch_cosmic_tbi,
            ch_gnomad_faf, ch_gnomad_faf_tbi
        )

        // Normals' germline calls from a somatic run (same VEP config, no COSMIC).
        if (germline_on_normals && stageRuns('variant_calling')) {
            ANNOTATE_VCF_GERMLINE(
                ch_germline_vcf,
                ch_vep_cache_dir,
                ch_ref_fasta, ch_ref_fai,
                ch_vep_plugins_dir,
                ch_clinvar, ch_clinvar_tbi,
                channel.value([]), channel.value([]),
                ch_gnomad_faf, ch_gnomad_faf_tbi
            )
        }

        // TMB from the VEP-annotated calls (needs consequence annotations).
        if (params.compute_tmb && params.analysis == 'somatic') {
            // Callable territory from the tumour and normal depth classes when the metrics
            // stage ran in this run; otherwise the BED territory alone (re-entry runs).
            // combine (not join) on the normal: several tumours can share one normal.
            ch_tmb_vcf = ANNOTATE_VCF.out.annotated_vcf.map { meta, vcf, tbi, _html -> tuple(meta, vcf, tbi) }
            if (stageRuns('metrics')) {
                ch_depth   = ch_depth_quantized.map { meta, q -> tuple(meta.sample_id, q) }
                ch_tmb_vcf = ch_tmb_vcf
                    .map     { meta, vcf, tbi -> tuple(meta.sample_id, meta, vcf, tbi) }
                    .join    (ch_depth, by: 0)
                    .map     { _sid, meta, vcf, tbi, tq -> tuple(meta.normal_id, meta, vcf, tbi, tq) }
                    .combine (ch_depth, by: 0)
                    .map     { _nid, meta, vcf, tbi, tq, nq -> tuple(meta, vcf, tbi, tq, nq) }
            } else {
                ch_tmb_vcf = ch_tmb_vcf.map { meta, vcf, tbi -> tuple(meta, vcf, tbi, [], []) }
            }
            COMPUTE_TMB(
                ch_tmb_vcf,
                ch_interval_bed,
                toOptionalFileChannel(params.tmb_coding_bed)
            )
        }
        ANNOTATE_VCF.out.annotated_vcf
            .map { _meta, _vcf, _tbi, html -> html }
            .collect(sort: true)
            .map { files -> tuple("06_annotation", files) }
            | MULTIQC_ANNOT
    }

    // -- STAGE 7 classification: classification and reports ---------
    if (stageRuns('classification')) {

        def cold_start = (start_rank == stageRankMap()['classification'])
        // Consolidated reports consume SOMATIC_BLOCK/GERMLINE_BLOCK outputs, which
        // exist only if variant_calling ran in THIS run (start_at=annotation too).
        def calls_this_run = (start_rank <= stageRankMap()['variant_calling'])

        ch_annot_for_classify = cold_start
            ? ch_active_vcf
            : ANNOTATE_VCF.out.annotated_vcf.map { meta, vcf, tbi, _html -> tuple(meta, vcf, tbi) }
        // Somatic runs with germline_on_normals: the normals' germline calls are classified
        // (germline rules) and converted to MAF too.
        if (germline_on_normals && !cold_start && stageRuns('annotation') && calls_this_run)
            ch_annot_for_classify = ch_annot_for_classify.mix(
                ANNOTATE_VCF_GERMLINE.out.annotated_vcf.map { meta, vcf, tbi, _html -> tuple(meta, vcf, tbi) })

        if (params.classify_variants) {
            CLASSIFY_VARIANTS(ch_annot_for_classify, toOptionalFileChannel(params.cosmic_cgc_tsv),
                              toOptionalFileChannel(params.germline_af_thresholds))
        }
        // MAF needs only the annotated VCF - independent of the pre-classification toggle.
        if (params.generate_maf && params.analysis in ['somatic', 'germline']) {
            VCF2MAF(ch_annot_for_classify, ch_ref_fasta, ch_ref_fai)
        }

        // TMB is optional in the report: with compute_tmb=false the TMB section reads
        // "(no TMB output)" instead of the whole report being dropped.
        def can_report = calls_this_run && params.consolidate_somatic_reports &&
                         params.classify_variants && params.analysis == 'somatic'
        if (can_report) {
            SOMATIC_REPORT_BLOCK(
                CLASSIFY_VARIANTS.out.classification.filter { m, _t -> m.type == 'DT' },
                SOMATIC_BLOCK.out.cnv,
                SOMATIC_BLOCK.out.loh,
                SOMATIC_BLOCK.out.sv_vcf,
                SOMATIC_BLOCK.out.msi_report,
                params.compute_tmb ? COMPUTE_TMB.out.tmb : channel.empty(),
                ch_annotsv_dir,
                SOMATIC_BLOCK.out.contamination
            )
        } else if (!calls_this_run && params.consolidate_somatic_reports) {
            log.info "Report skipped on cold-start; use -resume for the consolidated report."
        } else if (params.consolidate_somatic_reports && params.analysis == 'somatic' && !params.classify_variants) {
            log.info "Somatic consolidated report skipped (params.classify_variants = false)."
        }

        // -- Germline consolidated report (Phase 1) ----------------
        //   Needs SV/CNV/ROH produced this run, so it is gated on !cold_start
        //   exactly like the somatic report. SNV source is the VEP-annotated
        //   VCF (for ClinVar/CLIN_SIG in the report's SNV section).
        def germ_report_mode = (params.analysis == 'germline' || germline_on_normals)
        // AnnotSV is optional: without it only the SV/CNV sections are empty.
        def can_germline_report = calls_this_run && params.consolidate_germline_reports && germ_report_mode
        if (can_germline_report) {
            ch_germ_snv_annot = (params.analysis == 'germline'
                    ? ANNOTATE_VCF.out.annotated_vcf
                    : ANNOTATE_VCF_GERMLINE.out.annotated_vcf)
                .map { meta, vcf, tbi, _html -> tuple(meta, vcf, tbi) }
            GERMLINE_REPORT_BLOCK(
                ch_germ_snv_annot,
                GERMLINE_BLOCK.out.sv_vcf,
                GERMLINE_BLOCK.out.cnv_vcf,
                GERMLINE_BLOCK.out.roh,
                ch_annotsv_dir,
                toOptionalFileChannel(params.cosmic_cgc_tsv),
                toOptionalFileChannel(params.germline_report_genes)
            )
        } else if (params.consolidate_germline_reports && germ_report_mode && !calls_this_run) {
            log.info "Germline report skipped on cold-start; use -resume for the consolidated report."
        }
    }
}