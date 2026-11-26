# PAPIO - Pan-omic Analysis PIpeline for Oncology (v1.0.0)

PAPIO is a Nextflow DSL2 workflow for tumor exome and RNA sequencing. From FASTQ files it produces
germline variants, somatic SNVs and indels (a four-caller ensemble of Mutect2, Strelka2, MuSE and VarScan2, with two
voting modes: Mutect2-anchored, or any two callers), copy number and LOH with purity and ploidy,
structural variants, microsatellite instability, tumor mutational burden and RNA gene fusions, with
per-sample QC and one consolidated report per patient. It runs on an HPC cluster (LSF) with Singularity/Apptainer images.

This README first shows how to install PAPIO and check it on a small test dataset, then how to run it on your own data
and what it produces. Building the full reference set and reproducing the manuscript benchmarks are in Appendices A and B.
The demonstration data use the Agilent SureSelect V6 + UTR exome kit (SEQC2 HCC1395/HCC1395BL, GIAB HG002);
settings for SureSelect V7, used for our own RCC cohort, are shown as commented alternatives.

---

## Contents

1. [Analysis arms](#analysis-arms)
2. [Requirements](#requirements)
3. [Repository layout](#repository-layout)
4. [Quick start: install and test run](#quick-start-install-and-test-run)
5. [Containers](#containers)
6. [Samplesheets](#samplesheets)
7. [Parameter file](#parameter-file)
8. [Running PAPIO](#running-papio)
9. [Processes in each arm](#processes-in-each-arm)
10. [Outputs](#outputs)
11. [Troubleshooting](#troubleshooting)
12. [Appendix A. Reference resources (step by step)](#appendix-a-reference-resources-step-by-step)
13. [Appendix B. Demonstration data and benchmarking](#appendix-b-demonstration-data-and-benchmarking)
14. [Code availability and citation](#code-availability-and-citation)

## Analysis arms

| Arm | Set in the parameter file | Input | Main results |
|---|---|---|---|
| Somatic | `assay: wes` (or `wgs`), `analysis: somatic` | Matched tumor (DT) and normal (DN) DNA; several tumors of one patient are each paired with the patient's normal | Somatic SNVs/indels from the caller ensemble (VCF, MAF, classification table), copy number and LOH with purity and ploidy (Sequenza, CNVkit), somatic SVs (Manta), MSI, TMB (WES only), contamination, optional mutational signatures, one report per tumor |
| Germline | `assay: wes` (or `wgs`), `analysis: germline` | Normal DNA (DN) | Germline SNVs/indels (VCF, MAF, classification table), germline SVs (Manta), optional CNVs (gCNV, at least 30 normals of one kit), optional ROH (WGS) and HLA class I typing, one report per sample |
| Germline within a somatic run | `analysis: somatic`, `germline_on_normals: true` | As somatic | Somatic results plus germline SNVs/indels, classification and a germline report for each matched normal, written in the `somatic/` folders under the normal's sample ID |
| RNA fusion | `assay: rna`, `analysis: rna_fusion` | RNA-seq (RN or RT) | STAR alignment and RNA QC, fusions from Arriba and STAR-Fusion (FusionCatcher optional), FusionInspector evidence, tiered fusion table and report; optional Salmon quantification, TCR/BCR (TRUST4) and HERV (Telescope) |

The QC block (coverage, alignment metrics, sex check, contamination, sample-swap check, QC gate) runs in both DNA arms.
Tumor-only analysis (`analysis: tumor_only`) is not implemented; the pipeline stops with an error.

---

## Requirements

| Item | Version used | Note |
|---|---|---|
| Nextflow | 26.04.1 or later | Strict-syntax parser; earlier releases fail at parse time |
| Java | 17 or later | For Nextflow |
| Singularity / Apptainer | 3.8+ / 1.2+ | All tools run from images |
| Scheduler | LSF (default), SLURM (`-profile slurm`) or one server (`-profile local`) | Queues set with `lsf_queue_*` (SLURM: partition names) |
| Disk | ~250 GB references; ~150 GB work space per tumor-normal exome pair | Work directory on scratch |
| Memory | up to 93 GB for duplicate marking of a deep normal exome | Other steps ≤ 64 GB |

Typical run time for one tumor-normal exome pair: about 3.5 h wall clock and about 25 CPU-hours.

---

Per-process CPU, memory and time, and the resources each analysis needs, are listed in `nextflow/docs/RESOURCES.md`.
The test run below needs one Linux server with 16 CPUs and 64 GB RAM, without a scheduler.

---

## Repository layout

```
genomics_pipeline/
|-- nextflow/                 # pipeline code (main.nf, nextflow.config, workflows_*.nf, components_*.nf, utils.nf)
|   |-- setup.sh              # one-command setup: --test (test bundle) or --full (all public references)
|   |-- bin/                  # the data logic of the processes (Bash, Python, R, awk; one script per task, each with
|   |                         #   a usage header); processes keep only the tool call. Nextflow puts this folder on PATH
|   |                         #   inside every task and container, and a change to a script re-runs its tasks on -resume
|   |-- conf/                 # local.config (one server), slurm.config (SLURM)
|   |-- params/               # ready-made params files: somatic, somatic_sensitive, germline, rna_fusion
|   |-- scripts/              # run by the user, not the pipeline: check_installation.sh, download_demo_fastq.sh,
|   |                         #   compare_to_expected.sh, make_reference_manifest.sh, make_demo_subset.sh,
|   |                         #   make_test_bundle.sh, prepare_gnomad_faf.sh, gtf_to_transcripts.sh
|   |-- docs/                 # RESOURCES.md (resources and compute), DEMO.md (test run and full benchmark)
|   |-- tests/expected/       # expected outputs of the test run
|   |                         # optional PoN builder: main_pon.nf + components_pon.nf + pon.config (Appendix A, step 4.5)
|-- make_figures.py           # Figures 2 and 4 (Appendix B, step 14)
|-- benchmark_plots.py        # Figures 3 and 5 (Appendix B, step 14)
|-- trace_table.py            # Table 1, resource use (Appendix B, step 14)
|-- fusion_truth_match.py     # fusion calls vs a published list (Appendix B, step 13)
|-- make_refflat.py           # refFlat for Picard RNA metrics (Appendix A, section 8.6)
|-- containers/               # Singularity images (*.sif)
|-- references/               # reference data (Appendix A; same layout as setup.sh --full)
|-- demo/                     # demonstration FASTQs + samplesheets (nextflow/scripts/download_demo_fastq.sh ../demo)
|-- runs/                     # work_root / outdir of the runs (nextflow/params/*.yaml)
|-- README.md
```

Run all commands from `genomics_pipeline/`.

```bash
git clone https://github.com/chauhan2987/PAPIO.git genomics_pipeline
cd genomics_pipeline
```

---

## Quick start: install and test run

The test run checks that PAPIO is installed correctly and that every arm runs from FASTQ to report. It uses demo reads
and a small test bundle (tens of MB: the demo windows of one chromosome and a mini CTAT library; setup builds their
bwa-mem2 and STAR indexes), so its calls are **not** comparable with the manuscript's
benchmark numbers; Appendix B reproduces those.

```bash
# 1. containers + test bundle + ready-to-run params files
bash nextflow/setup.sh --test ~/papio_test          # add --bundle <file> for a locally built bundle

# 2. check that every file, image and setting the params file names is present
bash nextflow/scripts/check_installation.sh ~/papio_test/params_test_somatic.yaml

# 3. run the three arms on one server
nextflow run nextflow/main.nf -params-file ~/papio_test/params_test_somatic.yaml  -profile local
nextflow run nextflow/main.nf -params-file ~/papio_test/params_test_germline.yaml -profile local
nextflow run nextflow/main.nf -params-file ~/papio_test/params_test_rna.yaml      -profile local

# 4. compare with the expected outputs
bash nextflow/scripts/compare_to_expected.sh check ~/papio_test/runs/test_somatic/Result  nextflow/tests/expected/somatic
bash nextflow/scripts/compare_to_expected.sh check ~/papio_test/runs/test_germline/Result nextflow/tests/expected/germline
bash nextflow/scripts/compare_to_expected.sh check ~/papio_test/runs/test_rna/Result      nextflow/tests/expected/rna
```

- Use `-profile slurm` on a SLURM cluster; without a profile, tasks go to LSF. Set the CPU and memory limits of your
  machine in `nextflow/conf/local.config`.
- `compare_to_expected.sh` passes when every expected file exists, variant and row counts are within 5 % (or 2 records),
  and the set of PASS fusions is identical.
- Some steps are switched off in the test run because they need the whole exome or licensed data (Sequenza, TMB,
  signatures, COSMIC, AnnotSV, FusionCatcher); CNVkit runs, but one chromosome gives no meaningful copy-number profile.
  The reasons are listed in `nextflow/docs/DEMO.md`.

For a full installation, `bash nextflow/setup.sh --full <dir>` first checks that every download URL answers, then
downloads every public GRCh38 resource of Appendix A (sections 2–8), restricts the GATK and somalier VCFs to the
GRCh38.d1.vd1 contigs, builds the bwa-mem2, STAR, Sequenza and MSIsensor-pro indexes, the Salmon transcript FASTA and the gene, MANE CDS,
chromosome-arm, reportable-gene coverage, refFlat and rRNA files, fetches the somalier ancestry, TRUST4, HERV (Telescope)
and refPhi files, and writes params files under `<dir>/params` with the paths filled in (folder layout of
`nextflow/params/`). With `--kit <ID> --kit-regions <ID>_Regions.bed --kit-probes <ID>_Covered.bed` (the SureDesign
design BEDs) it also builds the capture-kit files of sections 4 and 5.2: interval lists, calling shards, call regions
and gCNV intervals; with `--cosmic` it downloads and prepares COSMIC v104 (section 7.4) with the login in
`COSMIC_EMAIL` / `COSMIC_PASSWORD`, or asks for it. About 250 GB are downloaded; the bwa-mem2 index needs 64 GB of RAM.
A re-run resumes where it stopped. `--no-fusioncatcher`, `--no-annotsv` and `--no-gnomad` skip the large optional
resources and `--no-salmon` the Salmon transcript FASTA (the matching params are then set to null), `--check-urls` only tests the URLs, and
`PAPIO_SETUP_CONF=<file>` can point the URLs at a local mirror. Optional resources that are absent (CNVkit panel of
normals, VerifyBamID2, an in-house Mutect2 panel of normals from `main_pon.nf`) are set to null in the params files.
The reportable-gene list is the example of Appendix A, step 6.2; edit `<dir>/references/annotation/reportable_genes.txt`
and re-run setup to rebuild its coverage BED. With the Salmon, TRUST4 and HERV resources set, the RNA arm also runs
Salmon, TRUST4 and Telescope; set `salmon_transcripts`, `trust4_ref_dir` or `telescope_gtf` to null to skip them.
`bash nextflow/scripts/make_reference_manifest.sh
<params.yaml> <out.tsv>` records the version and checksum of each resource used in a run.

## Containers

The pipeline looks for these files in `containers/` (`container_dir`). The names are set in `nextflow/nextflow.config`.

| Parameter | File | Tools |
|---|---|---|
| `container_align` | `bwa_samtools_fastp_fgbio_bcftools.sif` | fastp 1.3.3, bwa-mem2 2.3, samtools/bcftools 1.23.1, fgbio 4.1.0 |
| `container_gatk` | `gatk.sif` | GATK 4.6.2.0, Picard 3.4.0, bedtools 2.30.0|
| `container_qc` | `multiqc_mosdepth_verifybam2.sif` | MultiQC 1.35, mosdepth 0.3.14, VerifyBamID2 2.0.2 |
| `container_somalier` | `somalier.sif` | somalier 0.3.2 |
| `container_manta` | `manta.sif` | Manta 1.6.0 |
| `container_strelka` | `strelka2.sif` | Strelka2 2.9.10 |
| `container_strelka_af` | `strelka_tag_af.sif` | pysam, bcftools (adds AF/AD to Strelka2 calls) |
| `container_muse` | `muse_bcftools.sif` | MuSE 2.1.2, bcftools 1.23.1 |
| `container_varscan` | `varscan.sif` | VarScan 2.4.6 |
| `container_deepvariant` | `deepvar_freebayes.sif` | DeepVariant 1.10.0, freebayes 1.3.10 (unused optional germline caller) |
| `container_msisensor` | `msisensorpro.sif` | MSIsensor-pro 1.2.0 |
| `container_sequenza` | `sequenza.sif` | sequenza-utils / R sequenza 3.0.0 |
| `container_annotation` | `bcftool_cnvkit_vep_snpeff.sif` | Ensembl VEP 115.2, CNVkit 0.9.13, bcftools 1.23.1, Python 3.10.20|
| `container_annotsv` | `annotsv.sif` | AnnotSV 3.5.3 |
| `container_vcf2maf` | `vcf2maf.sif` | vcf2maf 1.6.22 |
| `container_signatures` | `sigprofiler.sif` | SigProfilerAssignment 1.1.3 |
| `container_optitype` | `optitype_samtools.sif` | Optitype 1.5.0, samtools 1.23.1 |
| `container_star` | `star_samtools.sif` | STAR 2.7.11b, samtools 1.23.1|
| `container_salmon` | `salmon.sif` | Salmon 1.10.3|
| `container_arriba` | `arriba.sif` | Arriba 2.5.1 |
| `container_trust4` | `trust4.sif` | Trust4 1.1.10 |
| `container_telescope` | `telescope.sif` | Telescope 1.0.4.1 |
| `container_starfusion` | `starfusion_fusioninspector_ctat.sif` | STAR-Fusion 1.15.1, FusionInspector 2.10.0, CTAT (GRCh38 GENCODE v44, Oct 2023 plug-n-play)|
| `container_fusioncatcher` | `fusioncatcher.sif` | FusionCatcher 1.33 |
| (benchmarking only) | `rtg_tools.sif` | RTG Tools 3.13 (vcfeval) |

---

## Samplesheets

One CSV for DNA (`dna_samplesheet`) and one for RNA (`rna_samplesheet`), with the same columns.
One row per FASTQ pair (lane or chunk). Rows that give the same sample name (below) are merged into one sample; they
must agree on `sample_type`, `sex` and `batch` (DNA) or `sample_type` and single/paired-end (RNA), and each row needs
its own `lane_num`/`chunk`, otherwise the run stops with an error.

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier; links tumor, normal and RNA of one patient |
| `sample` | Optional. Label for several samples of one patient and type (multi-region or longitudinal tumors, e.g. `T1`, `T2`, `relapse`); letters, digits and `-` only. Leave blank (or omit the column) for one sample per type |
| `sex` | `0` (or `F`, `female`) female, `1` (or `M`, `male`) male, `NA` unknown (NA switches off the sex check; Sequenza then takes the sex from the normal's chrX depth) |
| `tumor_type` | Free-text label, e.g. `BRCA`, `RCC`; `NA` for normal-only samples |
| `batch` | Library or batch identifier (`LB` in the read group) |
| `sample_type` | DNA sheet: `DT` tumor, `DN` normal. RNA sheet: `RN` or `RT` (both get the suffix `_R`, so give two RNA samples of one patient distinct `sample` labels) |
| `run_num`, `flowcell`, `lane_num` | Read-group fields (`PU`) |
| `chunk` | `001`, `002`, … for lanes or replicate libraries of one sample |
| `read1`, `read2` | FASTQ paths (absolute, or relative to `fastq_root`); `read2` is required for DNA and may be left blank for single-end RNA |

Sample names are `<patient_id>_T`, `<patient_id>_N` and `<patient_id>_R`, or `<patient_id>_<sample>` when `sample` is
filled. In the somatic arm every tumor of a patient is paired with that patient's normal, so each patient needs exactly
one normal (DN) sample; reports of a labelled tumor are named after its sample ID.

Demonstration DNA samplesheet (SEQC2; somatic arm):

```
patient_id,sex,tumor_type,batch,sample_type,run_num,flowcell,lane_num,chunk,read1,read2
valid1,0,BRCA,IL_1,DN,SRR7890874,SRR7890874,1,001,SRR7890874_1.fastq.gz,SRR7890874_2.fastq.gz
valid1,0,BRCA,IL_1,DT,SRR7890883,SRR7890883,1,001,SRR7890883_1.fastq.gz,SRR7890883_2.fastq.gz
```

Germline demonstration (GIAB HG002; `analysis: germline`):

```
patient_id,sex,tumor_type,batch,sample_type,run_num,flowcell,lane_num,chunk,read1,read2
hg002,1,NA,NOVO_1,DN,HG002_rep1,HG002_rep1,1,001,HG002_novogene_wes_18g_rep1_R1.fq.gz,HG002_novogene_wes_18g_rep1_R2.fq.gz
```

RNA demonstration (HCC1395 bulk RNA-seq; one sample):

```
patient_id,sex,tumor_type,batch,sample_type,run_num,flowcell,lane_num,chunk,read1,read2
rna,0,BRCA,IL_1,RN,SRR9134735,SRR9134735,1,001,SRR9134735_1.fastq.gz,SRR9134735_2.fastq.gz
```

When DNA and RNA come from the same patient, give them the same `patient_id` so the results can be joined; the
demonstration RNA is a separate study, so it has its own (`rna`).

---

## Parameter file

All site- and run-specific values go in one YAML file passed with `-params-file`.
`container_dir`, `work_root`, `outdir` and `bind_paths` must be set here, not in a `-c` config, because they are read
while the configuration loads. Environment variables such as `${HOME}` are **not** expanded in a params file:
use absolute paths or `${projectDir}`, which points to `nextflow/`.

`params_demo.yaml` (somatic demonstration; V6 + UTR active, own V7 data commented):

```yaml
# ---------------- Run control ----------------
assay:                  wes             # wes | rna
analysis:               somatic         # somatic | germline | rna_fusion
start_at:               mapping         # mapping | markdup | metrics | bqsr | variant_calling | annotation | classification
stop_at:                classification
write_manifest:         true
preflight_references:   true            # check references before any job is submitted

# ---------------- Input / output ----------------
pipeline_root:          ${projectDir}                  # the nextflow/ folder (holds main.nf)
container_dir:          ${projectDir}/../containers
fastq_root:             ${projectDir}/../demo/fastq
work_root:              ${projectDir}/../runs/somatic_demo
outdir:                 ${projectDir}/../runs/somatic_demo/Result
dna_samplesheet:        ${projectDir}/../demo/samplesheets/dna_somatic.csv
rna_samplesheet:        null

# directories the containers may read/write (add the folders of your own FASTQs and runs)
bind_paths:             ['${projectDir}/../references', '${projectDir}/../runs', '${projectDir}/../demo']

# ---------------- Reference genome ----------------
# (GDC GRCh38.d1.vd1)
ref_index_dir:          ${projectDir}/../references/genome/bwa_mem2
ref_fasta:              ${projectDir}/../references/genome/GRCh38.d1.vd1.fa
ref_fai:                ${projectDir}/../references/genome/GRCh38.d1.vd1.fa.fai
ref_dict:               ${projectDir}/../references/genome/GRCh38.d1.vd1.dict

# DNA arms - shared (somatic and germline)
# ---------------- Library ----------------
dedup_mode:             markdup         # umi_markdup| umi_consensus | markdup (demonstration: no UMIs)
umi_read_structure:     null
# dedup_mode:           umi_markdup     # own data: inline UMIs
# umi_read_structure:   '3M2S+T' # per-mate: 3bp UMI, 2bp linker skip, rest template
umi_read_structure_r2:  null     # R2 structure if it differs from R1 (e.g. '+T' for UMI on R1 only); null = same as R1

bqsr:                   false
skip_trimming:          false
make_cram:              false     # BAM → CRAM conversion (storage only). false = keep BAM.

# ---------------- Capture: demonstration = SureSelect V6 + UTR (S07604624) ----------------
interval_bed:           ${projectDir}/../references/capture/S07604624/S07604624_targets.bed
call_regions_gz:        ${projectDir}/../references/capture/S07604624/S07604624_targets_pad100.bed.gz
bait_interval_list:     ${projectDir}/../references/capture/S07604624/S07604624_baits.interval_list
target_interval_list:   ${projectDir}/../references/capture/S07604624/S07604624_targets.interval_list

# ---------------- Capture: own data = SureSelect V7 (S31285117) ----------------
# interval_bed:           ${projectDir}/../references/capture/S31285117/S31285117_targets.bed
# call_regions_gz:        ${projectDir}/../references/capture/S31285117/S31285117_targets_pad100.bed.gz
# bait_interval_list:     ${projectDir}/../references/capture/S31285117/S31285117_baits.interval_list
# target_interval_list:   ${projectDir}/../references/capture/S31285117/S31285117_targets.interval_list

# ---------------- Known sites (BQSR, germline) ----------------
# (restricted to the GRCh38.d1.vd1 contigs, Appendix A step 3.3; need .tbi)
dbsnp:                  ${projectDir}/../references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz
mills:                  ${projectDir}/../references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_indels:           ${projectDir}/../references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz

# ---------------- Sample QC ----------------
somalier_sites:         ${projectDir}/../references/somalier/sites.hg38.vcf.gz
somalier_ancestry_labels: ${projectDir}/../references/somalier/ancestry-labels-1kg.tsv # ancestry-labels-1kg.tsv; null = skip ancestry
somalier_1kg_dir:       ${projectDir}/../references/somalier/1kg-somalier # extracted 1kg-somalier/ directory

verifybamid_db:         null
coverage_report_bed:    ${projectDir}/../references/annotation/ensembl115_reportable_genes_MANE_CDS.bed

# ---------------- Annotation (both DNA arms) ----------------
vep_cache_dir:          ${projectDir}/../references/vep_cache
vep_plugins_dir:        ${projectDir}/../references/vep_plugins   # NMD.pm (Appendix A, step 7.7)
gnomad_faf_vcf:         ${projectDir}/../references/gnomad/variantsGnomADv4.1Joint.vcf.bgz   # both DNA arms (Appendix A, step 7.6)
clinvar_vcf:            ${projectDir}/../references/clinvar/clinvar_YYYYMMDD.vcf.gz   # YYYYMMDD = ClinVar release date (Appendix A, step 7.3)
annotsv_dir:            ${projectDir}/../references/annotsv
cosmic_cgc_tsv:         ${projectDir}/../references/cosmic/Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz

# ---------------- Somatic arm (analysis: somatic) ----------------
scatter_intervals:      ${projectDir}/../references/capture/S07604624/S07604624_scatter_24_pad100
# scatter_intervals:    ${projectDir}/../references/capture/S31285117/S31285117_scatter_24_pad100   # own data (V7)
germline_resource:      ${projectDir}/../references/gatk/af-only-gnomad.hg38.vcf.gz
pileup_variants:        ${projectDir}/../references/gatk/small_exac_common_3.hg38.vcf.gz
# choose either the 1000G PoN or your own PoN
panel_of_normals:       ${projectDir}/../references/gatk/1000g_pon.hg38.vcf.gz
# panel_of_normals:     ${projectDir}/../references/capture/S31285117/S31285117_inhouse_pon.vcf.gz   # optional: main_pon.nf (Appendix A, step 4.5)
orientation_bias:       true
consensus_backbone:     mutect2   # mutect2: vote on Mutect2 PASS calls; union: vote over all callers (Mutect2 is one vote)
cnvkit_pon:             null      # or ${projectDir}/../references/capture/<KIT>/<KIT>_cnvkit_pon.cnn (same kit)
msisensor_list:         ${projectDir}/../references/msi/GRCh38.d1.vd1.msisensor.list
sequenza_gc_wiggle:     ${projectDir}/../references/sequenza/GRCh38.d1.vd1.gc50Base.wig.gz
gene_bed:               ${projectDir}/../references/annotation/ensembl115_genes.bed
arm_bed:                ${projectDir}/../references/annotation/GRCh38_arms.bed
cosmic_vcf:             ${projectDir}/../references/cosmic/Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.chr_dedup.vcf.gz
tmb_coding_bed:         ${projectDir}/../references/annotation/ensembl115_MANE_CDS.bed
germline_on_normals:    false     # true = germline report for each matched normal

# ---------------- Germline arm (analysis: germline) ----------------
germline_scatter_intervals: ${projectDir}/../references/capture/S07604624/S07604624_scatter_150_pad100
# germline_scatter_intervals: ${projectDir}/../references/capture/S31285117/S31285117_scatter_150_pad100   # own data (V7)
germline_caller:        haplotypecaller   # haplotypecaller | deepvariant
germline_filter:        nv                # nv (any sample number NVScoreVariants) | none (unfiltered GenotypeGVCFs) | vqsr (≥ 30 samples)
vqsr_min_samples:       30                # minimum samples before VQSR is allowed (else use 'nv')
run_gcnv:               false             # true only with ≥ 30 normals of one kit
gcnv_min_samples:       30
germline_ba1_af:        0.05
hapmap:                 ${projectDir}/../references/gatk/hapmap_3.3.hg38.vcf.gz
omni:                   ${projectDir}/../references/gatk/1000G_omni2.5.hg38.vcf.gz
thousandg_snps:         ${projectDir}/../references/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz
germline_report_genes:  ${projectDir}/../references/annotation/reportable_genes.txt
# gCNV inputs: only with run_gcnv: true and ≥ 30 normals of one kit (not used for the demonstration);
# the paths below are for own data (SureSelect V7, S31285117)
# contig_ploidy_priors:     ${projectDir}/../references/annotation/GRCh38_contig_ploidy_priors.tsv
# preprocessed_intervals:   ${projectDir}/../references/capture/S31285117/gcnv/S31285117_preprocessed_intervals.interval_list
# gcnv_scatter_intervals:   ${projectDir}/../references/capture/S31285117/gcnv/S31285117_scatter_24_pad250
# gcnv_annotated_intervals: ${projectDir}/../references/capture/S31285117/gcnv/S31285117_annotated_intervals.tsv

run_optitype:           true              # HLA class I typing (OptiType)

# ---------------- RNA fusion arm (analysis: rna_fusion) ----------------
rna_use_ctat_reference: false
gtf:                    ${projectDir}/../references/gencode/gencode.v49.annotation.gtf
star_index:             ${projectDir}/../references/star/StarIndex_gencode_v49
salmon_transcripts:     ${projectDir}/../references/star/gencode.v49.annotation.transcripts.fa   # Salmon (Appendix A, step 8.7)
refflat:                ${projectDir}/../references/rna_qc/gencode.v49.refFlat
rrna_interval_list:     ${projectDir}/../references/rna_qc/gencode.v49.rrna.interval_list
rna_strandedness:       auto              # auto (from STAR gene counts) | no | yes (forward) | reverse (dUTP, most kits)
ctat_lib:               ${projectDir}/../references/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir   # needed only for STAR-Fusion, FusionInspector, rna_use_ctat_reference, or when star_index is not given
arriba_blacklist:       ${projectDir}/../references/arriba/blacklist_hg38_GRCh38_v2.5.1.tsv.gz
arriba_known_fusions:   ${projectDir}/../references/arriba/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz
arriba_protein_domains: ${projectDir}/../references/arriba/protein_domains_hg38_GRCh38_v2.5.1.gff3
known_fusions_tsv:      ${projectDir}/../references/fusion/known_fusions_hg38_GRCh38_v2.5.1.tsv
use_fusioncatcher:      true
fusioncatcher_db:       ${projectDir}/../references/fusioncatcher/human_v102

# immune processes resources (TRUST4 + Telescope)
trust4_ref_dir:         ${projectDir}/../references/trust4
telescope_gtf:          ${projectDir}/../references/telescope/HERV_rmsk.hg38.v2.gtf

# ---------------- LSF ----------------
lsf_queue_short:        short
lsf_queue_medium:       medium
lsf_queue_high:         long
lsf_queue_long:         long
lsf_project:            null
```

For the other arms, copy the file and change only:

| Arm | Change |
|---|---|
| Germline (HG002) | `analysis: germline`, `dna_samplesheet:` the HG002 sheet, new `work_root` and `outdir` |
| RNA fusion | `assay: rna`, `analysis: rna_fusion`, `dna_samplesheet: null`, `rna_samplesheet:` the RNA sheet, new `work_root` and `outdir` |

---
Ready-made params files for each arm, and a sensitive somatic preset (`params_somatic_sensitive.yaml`), are in
`nextflow/params/`.

### Checking a library for inline UMIs

Look at the base composition of the first 15 cycles of R1 and R2, for tumor and normal, in several patients.
A random block followed by a fixed block means an inline UMI. In the RCC cohort: 3 random bases, 2 fixed, then genome,
which is `umi_read_structure: 3M2S+T` (3 bases UMI, skip 2, rest template).

```bash
FQ=<sample>_R1.fq.gz        # any R1 or R2 file of the library
zcat $FQ | head -40000 | awk 'NR%4==2 {n++; for (i=1;i<=15;i++) c[i, substr($0,i,1)]++}
  END {for (i=1;i<=15;i++) printf "pos%2d  A:%5.1f  C:%5.1f  G:%5.1f  T:%5.1f\n", i,
       100*c[i,"A"]/n, 100*c[i,"C"]/n, 100*c[i,"G"]/n, 100*c[i,"T"]/n}'
```

## Running PAPIO

Launch from a login node or a long-running interactive job; Nextflow submits every task to LSF.

```bash
nextflow run ${HOME}/genomics_pipeline/nextflow/main.nf -c ${HOME}/genomics_pipeline/nextflow/nextflow.config \
  -params-file params_demo.yaml -profile production -resume
```

- Profiles: `production` (LSF, default settings), `clinical` (BQSR on, DeepVariant), `local` (one server),
  `slurm` (SLURM cluster), `bignode` (raises the memory ceiling for very deep samples; with `production`) and `test`
  (CI only: local, mapping to markdup). Profiles can be combined, e.g. `-profile production,slurm`.
- Before the first run, `bash nextflow/scripts/check_installation.sh params_demo.yaml` lists anything missing.
  The pipeline itself also stops at start when a required resource is missing (e.g. `dbsnp`, `hapmap`, `omni`,
  `thousandg_snps` and `mills` for HaplotypeCaller; `dbsnp` for MuSE).
- `-resume` reuses finished tasks from `work_root`; always add it when rerunning the same analysis. With Nextflow
  26.04 (strict parser) a changed parameter does **not** make cached tasks run again, so PAPIO stops a resumed run whose
  parameters differ from the run it resumes and lists them. Use a new `work_root` (or drop `-resume`) after changing
  parameters; `allow_param_change_on_resume: true` downgrades the stop to a warning.
- **Partial runs:** set `stop_at` (e.g. `variant_calling`). To restart later from published files, set `start_at`
  and the matching input folder (`input_bam_dir` from `markdup` on; `input_vcf_dir` for `annotation`, the folder
  holding `<S>_filtered.vcf.gz`, `<S>_genotyped.vcf.gz` or `<S>.vcf.gz`; `input_annot_dir` for `classification`, the
  folder holding `<S>_vep.vcf.gz`). BAMs that are already deduplicated (or UMI consensus) elsewhere enter with
  `dedup_mode: external` and `start_at: metrics` or later. The consolidated reports need upstream results of the same run,
  so they are not made when starting at `classification`; use `-resume` instead.
- **Process-specific resources:** put overrides in a small extra config and add it with a second `-c`:

```groovy
// more_time.config
process {
    withName: 'PREPROCESS_AND_ALIGN_DNA' { time = { 12.h * task.attempt }; queue = params.lsf_queue_medium }
    withName: 'METRICS'                  { time = { 5.h  * task.attempt }; queue = params.lsf_queue_medium }
}
```

A run writes `pipeline_trace.txt` and `pipeline_report.html` to `outdir`; `detailed_reports: true` adds
`pipeline_timeline.html` and `pipeline_dag.html`.

---

## Processes in each arm

Kind: **required** - always runs; **toggle** - switched by the parameter shown; **resource** - runs only when that
reference file is set. Stages run in the order listed.

### DNA stages shared by the somatic and germline arms

| Stage | Processes | Kind |
|---|---|---|
| 0 | `PREFLIGHT_REFERENCES` | toggle (`preflight_references`, on) |
| | `MANIFEST_ROW` → `WRITE_MANIFEST` | toggle (`write_manifest`, off) |
| 1 mapping | `PREPROCESS_AND_ALIGN_DNA` (fastp → bwa-mem2 → samtools fixmate) for `dedup_mode: markdup`, or `ALIGN_DNA_UMI` (fgbio FastqToBam → bwa-mem2 → ZipperBams) for the UMI modes; `MULTIQC_TRIMMING` | required |
| 2 markdup | `MERGE_AND_MARK_DUPLICATES` (`markdup`), `MERGE_AND_MARKDUP_UMI` (`umi_markdup`) or `MERGE_AND_CONSENSUS` (`umi_consensus`); `MULTIQC_MARKDUP` | required (one, by `dedup_mode`) |
| | `CRAM_CONVERT_DEDUP` | toggle (`make_cram`, off) |
| 3 metrics | `MOSDEPTH` → `INFER_SEX`; `SAMTOOLS_STATS`; `METRICS` (Picard); `MULTIQC_METRICS` | required |
| | `VERIFYBAMID2` (normals) | resource (`verifybamid_db`) |
| | `QC_GATE` → `QC_GATE_SUMMARY` | toggle (`qc_gate`, on) |
| | `SOMALIER_EXTRACT` → `SOMALIER_RELATE` → `CHECK_SAMPLE_SWAPS`; `MULTIQC_SOMALIER` | resource (`somalier_sites`) |
| | `SOMALIER_ANCESTRY` (normals) | resource (`somalier_ancestry_labels` and `somalier_1kg_dir`) |
| 4 bqsr | `BASE_RECALIBRATOR` (per shard on WGS) → `GATHER_BQSR_REPORTS` → `APPLY_BQSR`; `CRAM_CONVERT_BQSR` | toggle (`bqsr`, off; CRAM with `make_cram`) |
| 5 | `EVIDENCE_BAM` (read slices around given regions) | resource (`evidence_bed`) |

The QC gate reports only; it never stops a run. Tumors are checked against `qc_min_mean_coverage` (100) and
`qc_min_pct_target_30x` (90% of target bases at ≥ 30×); normals and germline samples against
`qc_min_mean_coverage_normal` (50) and `qc_min_pct_target_20x` (90% at ≥ 20×). WGS uses
`qc_wgs_min_mean_coverage_tumor` (60) and `qc_wgs_min_mean_coverage_normal` (30); on WGS the two breadth rows give the
share of the autosomes (chromosome length, N gaps included) at ≥ 30× for tumors and ≥ 20× for normals, checked against
`qc_wgs_min_pct_genome_30x` and `qc_wgs_min_pct_genome_20x` (0.80). The other checks are contamination
(`qc_max_contamination`, 0.03), on-target fraction (`qc_min_on_target`, 0.50), mapping rate (`qc_min_mapping_rate`,
0.95) and fold-80 penalty (`qc_max_fold80`, 3.0). Tumor contamination (GATK CalculateContamination) needs the variant
calls and is reported in the patient report, not in the gate.

The sex check compares chrX and chrY depth with the autosomes (`sex_female_min_xratio` 0.80, `sex_male_max_xratio`
0.65). When chrY has data it must agree: a female call needs chrY ≤ `sex_female_max_yratio` (0.10), a male normal
chrY ≥ `sex_male_min_yratio` (0.20; tumors can lose chrY). Disagreement (X0, XXY, contamination) or a chrX ratio between
the cut-offs is REVIEW.

### Somatic arm

| Stage | Processes | Kind |
|---|---|---|
| 5 variant calling | `MANTA` (SVs; its candidate indels feed Strelka2) | required |
| | `MUTECT2` (per shard) → `MERGE_MUTECT2` (merge + FilterMutectCalls) | required |
| | `LEARN_READ_ORIENTATION_MODEL` | toggle (`orientation_bias`, on) |
| | `GETPILEUPSUMMARIES` → `CALCULATE_CONTAMINATION` | resource (`pileup_variants`) |
| | `STRELKA2` → `TAG_STRELKA_AF` | toggle (`use_strelka2_in_consensus`, on) |
| | `MUSE_CALL` → `MUSE_SUMP` → `TAG_MUSE_AF` | toggle (`use_muse_in_consensus`, on; needs `dbsnp`) |
| | `VARSCAN` (per shard on WGS) → `VARSCAN_MERGE` → `TAG_VARSCAN_AF` | toggle (`use_varscan_in_consensus`, on) |
| | `ENSEMBLE_CONSENSUS` | toggle (`somatic_consensus`, on); voting set by `consensus_backbone` (see "Somatic consensus options") |
| | `FILTER_SOMATIC_VCF` | toggle (`somatic_variant_filter`, on) |
| | `MUTATIONAL_SIGNATURES` | toggle (`mutational_signatures`, off) |
| | `SEQUENZA` → `SEQUENZA_GENE_LOH` | resource (`sequenza_gc_wiggle`; the gene and arm tables need `gene_bed` and `arm_bed`) |
| | `CNVKIT` → `CNVKIT_CALL` | required; integer copy number uses the Sequenza purity (`-m clonal`) when Sequenza ran |
| | `MSISENSOR_PRO` → `CLASSIFY_MSI`; `MSISENSOR_GENE_MAP` | resource (`msisensor_list`); status = toggle (`classify_msi`, on); gene map needs `gene_bed` |
| 6 annotation | `ANNOTATE_VCF` (VEP; ClinVar, COSMIC and gnomAD FAF when given); `MULTIQC_ANNOT` | required |
| | `COMPUTE_TMB` | toggle (`compute_tmb`, on; on WES the capture BED is intersected with `tmb_coding_bed` when set; reported as NA on WGS) |
| 7 classification | `CLASSIFY_VARIANTS` (somatic mode) | toggle (`classify_variants`, on) |
| | `VCF2MAF` | toggle (`generate_maf`, on) |
| | `ANNOTATE_SV_SOMATIC` | resource (`annotsv_dir`) |
| | `CONSOLIDATE_PATIENT_REPORT` | toggle (`consolidate_somatic_reports`, on); needs `classify_variants` and the calls of the same run (or `-resume`) |

With `germline_on_normals: true` the germline processes below also run on each normal, and `ANNOTATE_VCF_GERMLINE`
annotates their calls (VEP; ClinVar and gnomAD FAF when given, no COSMIC).

### Germline arm

| Stage | Processes | Kind |
|---|---|---|
| 5 variant calling | `HAPLOTYPECALLER` (per shard) → `MERGE_GVCFS` → `GENOTYPE_GVCFS` | required with `germline_caller: haplotypecaller` (default) |
| | → `NV_SCORE_VARIANTS` → `FILTER_VARIANT_TRANCHES` | `germline_filter: nv` (default; needs `hapmap`, `omni`, `thousandg_snps`, `mills`) |
| | unfiltered `GENOTYPE_GVCFS` output | `germline_filter: none` |
| | Cohort path: `REBLOCK_GVCF` → `GENOMICSDB_IMPORT` → `GNARLY_GENOTYPER` or `GENOTYPE_GVCFS_JOINT` → `MERGE_COHORT_VCF` → `VQSR_SNP`, `VQSR_INDEL` → `APPLY_VQSR_SNP` → `APPLY_VQSR_INDEL` → `SPLIT_COHORT_VCF` | `germline_filter: vqsr` (needs ≥ `vqsr_min_samples`, default 30) |
| | `DEEPVARIANT` instead of all of the above | `germline_caller: deepvariant` |
| | `MANTA_GERMLINE` | required |
| | `GATK_COLLECT_READ_COUNTS` → `GATK_PLOIDY_MODEL` → `GATK_GCNV_CALLER` → `GATK_GCNV_POSTPROCESS` | toggle (`run_gcnv`, off; ≥ 30 normals of one kit) |
| | `BCFTOOLS_ROH` | toggle (`run_roh`, on) but skipped with a warning on WES and without `roh_af_file` |
| | `OPTITYPE` | toggle (`run_optitype`, off) |
| 6 annotation | `ANNOTATE_VCF` (VEP; ClinVar and gnomAD FAF when given, no COSMIC); `MULTIQC_ANNOT` | required |
| 7 classification | `CLASSIFY_VARIANTS` (germline mode) | toggle (`classify_variants`, on) |
| | `VCF2MAF` | toggle (`generate_maf`, on) |
| | `ANNOTATE_SV_MANTA`, `ANNOTATE_SV_GCNV` | resource (`annotsv_dir`) |
| | `CONSOLIDATE_GERMLINE_REPORT` | toggle (`consolidate_germline_reports`, on); needs the calls of the same run (or `-resume`); written without SV/CNV annotation when `annotsv_dir` is unset |

### RNA fusion arm

| Stage | Processes | Kind |
|---|---|---|
| 0 | `MANIFEST_ROW` → `WRITE_MANIFEST` | toggle (`write_manifest`, off) |
| 1 mapping | `STAR_ALIGN` (two-pass unless `star_twopass: false`, chimeric output) | required |
| | `SALMON_QUANT` (on STAR's transcriptome alignments) | resource (`salmon_transcripts`) |
| | `CRAM_CONVERT` | toggle (`rna_cram_convert`, off) |
| 3 metrics | `RNASEQ_METRICS` (Picard CollectRnaSeqMetrics; strandedness inferred from the STAR gene counts when `rna_strandedness: auto`); `SAMTOOLS_STATS`; `MULTIQC` | required |
| 5 variant calling | `ARRIBA` | toggle (`use_arriba`, on) |
| | `STARFUSION` (its own STAR run on the CTAT library; `rna_shared_star: true` is refused with STAR-Fusion, which needs a chimeric-junction format STAR cannot write alongside Arriba's input) | toggle (`use_starfusion`, on) |
| | `FUSIONCATCHER` | toggle (`use_fusioncatcher`, off; needs `fusioncatcher_db`) |
| | `FUSION_CONSENSUS` | required |
| | `FILTER_FUSION_CANDIDATES` → `FUSIONINSPECTOR` | toggle (`use_fusioninspector`, on) |
| | `RNA_IMMUNE_PREP` → `TRUST4`; `TELESCOPE` | resource (`trust4_ref_dir`; `telescope_gtf`) |
| | `EVIDENCE_BAM` | resource (`evidence_bed`) |
| 6 annotation | `ANNOTATE_FUSIONS` (known fusions + COSMIC CGC → tiers) | required |
| 7 classification | `CONSOLIDATE_FUSION_REPORT` | required |

### Panel-of-normals builder (`main_pon.nf`, optional)

`PON_MUTECT2` (per normal and shard) → `PON_MERGE_NORMAL` → `PON_GENOMICSDB_IMPORT` (per shard) → `PON_CREATE_SHARD`
→ `PON_MERGE_FINAL`. A separate workflow, run once per capture kit (Appendix A, step 4.5); its panel is then given to
the somatic arm as `panel_of_normals`.

### Which process assigns which label

| Labels | Column or field | Process | Arm |
|---|---|---|---|
| ONCOGENIC_STRONG, ONCOGENIC_MODERATE, VUS_CANCER_GENE, VUS, LIKELY_BENIGN, LOW_IMPACT | `PRE_CLASSIFICATION` in `<S>_preclassification.tsv` | `CLASSIFY_VARIANTS` (somatic mode) | Somatic SNVs/indels |
| PATHOGENIC_OR_LIKELY, BENIGN_OR_LIKELY, LIKELY_BENIGN, LOF_CANCER_GENE_REVIEW, VUS, LOW_IMPACT | `PRE_CLASSIFICATION` in `<S>_preclassification.tsv` | `CLASSIFY_VARIANTS` (germline mode) | Germline SNVs/indels |
| SOMATIC; GERMLINE, POSSIBLE_MOSAIC_OR_CH, UNRESOLVED | `ORIGIN` | `CLASSIFY_VARIANTS` | Somatic; germline (VAF-based, from `germline_low_vaf` and `germline_mosaic_min_dp`) |
| LOWCONSENSUS | VCF FILTER | `ENSEMBLE_CONSENSUS` | Somatic SNVs/indels below the caller threshold |
| LOWSUPPORT | VCF FILTER | `FILTER_SOMATIC_VCF` | Somatic SNVs/indels failing the VAF, alt-read, depth or matched-normal limits |
| MSI-H, MSS, MSI-INDETERMINATE, NA_LOW_SITES | `msi_status` | `CLASSIFY_MSI` | Somatic |
| TMB-High, TMB-Low (NA on WGS) | `tmb_status` | `COMPUTE_TMB` | Somatic |
| HET, LOH, HOMDEL, HEMIZYGOUS, SEX_UNKNOWN, UNKNOWN, NO_COVERAGE; CN_neutral_LOH, LOH_with_loss, LOH_with_gain | `LOH_Status`; `LOH_Detail` | `SEQUENZA_GENE_LOH` | Somatic, per gene |
| LOSS, GAIN, NEUTRAL | arm call in `<S>_arm_cn.tsv` | `SEQUENZA_GENE_LOH` | Somatic, per chromosome arm |
| HOMDEL, AMP, GAIN, LOSS | copy-number section of the patient report | `CONSOLIDATE_PATIENT_REPORT` | Somatic (Sequenza copy number, else CNVkit) |
| ACMG classes 1-5 | AnnotSV ranking | `ANNOTATE_SV_SOMATIC`, `ANNOTATE_SV_MANTA`, `ANNOTATE_SV_GCNV` (AnnotSV) | Germline SVs/CNVs; somatic SVs (germline ACMG criteria, for orientation only) |
| PASS_CONSENSUS, PASS_KNOWN, REVIEW | `Verdict` | `FUSION_CONSENSUS` | RNA fusions |
| KNOWN_ONCOGENIC, KNOWN_LOW_SUPPORT, LIKELY_ONCOGENIC, VUS, LIKELY_ARTIFACT, LOW_CONFIDENCE | `TIER` | `ANNOTATE_FUSIONS` | RNA fusions |
| PASS, WARN, INFO, SKIP, UNKNOWN | per-metric flag; `OVERALL` | `QC_GATE`, `QC_GATE_SUMMARY` | DNA QC |
| PASS, MISMATCH, REVIEW_TUMOR, REVIEW, SKIP | `sex_concordance` | `INFER_SEX` | DNA QC |
| PASS, FLAG_SWAP_LOW_RELATEDNESS, REVIEW_LOW_RELATEDNESS_IBS0_OK, FLAG_UNEXPECTED_RELATEDNESS, FLAG_POSSIBLE_DUPLICATE | `verdict` in `sample_swap_report.tsv` | `CHECK_SAMPLE_SWAPS` | DNA QC (somalier) |

All labels are computational pre-classifications for review, not AMP/ASCO/CAP or ACMG/AMP clinical tiers.

**Somatic consensus options.** Mutect2 always runs; Strelka2, MuSE and VarScan2 are added with their
`use_*_in_consensus` switches (all on by default). MuSE calls SNVs only, so indels are voted on by Mutect2, Strelka2
and VarScan2. Before the vote each caller's calls are reduced to its own passing set (Mutect2 and Strelka2 PASS, MuSE
without `muse_exclude_tier`, VarScan2 somatic status SS=2), split into single alleles and left-aligned so that the same
variant matches across callers.

| Parameter | Values (default first) | Effect |
|---|---|---|
| `somatic_consensus` | `true` \| `false` | `false`: the final call set is Mutect2 alone (also switch off the other callers to save their run time) |
| `consensus_backbone` | `mutect2` \| `union` | `mutect2`: only Mutect2 PASS calls are voted on, and a call is kept when enough callers, Mutect2 included, support it. `union`: true N-of-M voting over every caller, with Mutect2 as one vote; a call Mutect2 filtered is kept when the other callers outvote it (its Mutect2 filters are kept in `INFO/M2_FILTER`), and a call Mutect2 never made is taken from Strelka2, VarScan2 or MuSE (named in `INFO/CONS_SOURCE`) |
| `consensus_min_snv_callers` | `2` | callers needed for an SNV |
| `consensus_min_indel_callers` | `2` | callers needed for an indel |
| `consensus_filter_mode` | `tag` \| `remove` | `tag`: calls below the threshold stay in the VCF with FILTER `LOWCONSENSUS`; `remove`: they are dropped |
| `muse_exclude_tier` | `Tier5` | MuSE tier not counted as a vote |

A Mutect2 multi-nucleotide variant is kept as one record so that it is annotated at the codon level; another caller
counts as a vote for it only when it reports every base change of that variant. `INFO/CALLERS` and `INFO/NCALLERS` list
the supporting callers in both modes. `FILTER_SOMATIC_VCF` then applies the tumor VAF, alt-read, depth and matched-normal
limits (`filter_*` parameters) and, with `somatic_filter_mode: tag` (default), labels failures instead of dropping them.
VAF is alt / (ref + alt) from FORMAT/AD for every caller (Mutect2's FORMAT/AF is a model estimate). The MAF columns
`Tumor_VAF` and `Normal_VAF` use the same definition. `Tumor_Sample_Barcode` is the tumor sample ID (germline MAF: the
sample itself), `Matched_Norm_Sample_Barcode` the normal, and `Patient_ID` groups several tumors of one patient.

**TMB.** Nonsynonymous coding PASS variants per megabase of the territory: the capture BED, intersected with
`tmb_coding_bed` (when set) and with the bases at depth ≥ `filter_min_tumor_dp` (20) in the tumor and ≥ `tmb_min_normal_dp` (10)
in the normal (mosdepth depth classes from the metrics stage of the same run; thresholds round up to the next class
boundary 1, 10, 20, 30, 50, 100). The `territory` row records what was used, with the BED size before the depth step;
re-entry runs without the metrics stage use the BED territory. A sample with no base at those depths gets TMB NA.

**CNVkit integer copy number.** `CNVKIT_CALL` applies the Sequenza purity (`-m clonal --purity`) only when Sequenza's
ploidy is 1.6–2.4: CNVkit centres log2 on the median, so rescaling a polyploid tumor would give wrong copy numbers.
Otherwise CNVkit's threshold calling is used. The patient report takes copy number from Sequenza whenever it ran.

**RNA fusion consensus.** Arriba and STAR-Fusion run by default and FusionCatcher on request (`use_fusioncatcher`). A
fusion is PASS_CONSENSUS when at least `fusion_min_callers` callers (default 2) report the same gene pair (matched on
Ensembl gene IDs, so symbol differences between annotation releases do not split a fusion; the threshold is capped
at the number of callers that ran, and with Arriba alone its low-confidence calls stay REVIEW), and
PASS_KNOWN when only one caller reports it but the pair is in the known-fusion list (`known_fusions_tsv`); all other
calls are REVIEW. Read support is counted per caller and never added across callers, which see the same reads:
`TotalSupport` is the count of the best-supported caller and `SupportByCaller` lists each caller's count. A PASS_KNOWN
fusion is tiered KNOWN_ONCOGENIC only with at least `fusion_known_min_reads` supporting reads (default 3) and, when FusionInspector ran, its confirmation; otherwise it is KNOWN_LOW_SUPPORT. A pair is flagged
as read-through when Arriba calls it one, or when both genes lie on the same strand with the 5′ partner upstream of the
3′ partner within `fusion_readthrough_kb` (default 200 kb); duplications and inversions between neighbouring genes
(FGFR3::TACC3) are not.

**Oncogenes and tumor suppressors.** Gene roles come from the COSMIC Cancer Gene Census (`CGC_ROLE` column), and the
somatic rules follow the point system of the ClinGen/CGC/VICC oncogenicity standard (Horak et al., 2022). A null variant
(nonsense, frameshift, canonical splice site or start loss) in a tumor suppressor is labelled ONCOGENIC_MODERATE, and
ONCOGENIC_STRONG when the same change is seen in at least `cosmic_recurrence_min` COSMIC samples. In an oncogene, a
missense or in-frame change seen in at least `cosmic_recurrence_min` COSMIC samples is ONCOGENIC_MODERATE, and
ONCOGENIC_STRONG from `cosmic_strong_min` samples (default 50); null variants in oncogenes are left for review, because
gene activation rarely works by truncation. Stop-loss variants lengthen the protein and are not treated as null. A
nonsense or frameshift variant in the last exon, or in the last 50 nucleotides of the penultimate exon, of a multi-exon
transcript escapes nonsense-mediated decay and may leave a working protein (ClinGen PVS1 guidance), so it is graded one
step lower (VUS_CANCER_GENE, or ONCOGENIC_MODERATE when COSMIC-recurrent). The penultimate-exon case uses the Ensembl
VEP `NMD` plugin, which needs no data files: put `NMD.pm` in `vep_plugins_dir`; without it only the last exon is checked.
Variants without coding or splice-site impact and without other evidence are labelled LOW_IMPACT: they were not assessed
and are not benign calls. LIKELY_BENIGN and BENIGN_OR_LIKELY need evidence (population frequency or ClinVar).

**Population frequency.** With `gnomad_faf_vcf` set, the population frequency used for benign calls is the gnomAD v4.1
grpmax filtering allele frequency (`fafmax_faf95_max_j`), the field ClinGen recommends for the ACMG/AMP BA1 and BS1
criteria (file and download in Appendix A, step 7.6). Without it, VEP's `MAX_AF` is used, which includes
founder and small populations and so over-calls common variants. The column header names the source. The GATK
`af-only-gnomad` file is not used here: it holds only the global AF, which misses variants common in one population.
Germline BA1 and BS1 default to `germline_ba1_af` (0.05) and `classify_common_af` (0.01); gene-specific values from ClinGen
expert panel specifications (cspec.genome.network) can be given in `germline_af_thresholds`, a tab-separated file of
GENE, BA1, BS1.

**Germline report.** `CONSOLIDATE_GERMLINE_REPORT` uses the same sources as the classification: ClinVar from
`clinvar_vcf`, the population frequency named above (the column header gives the field) and the VEP PICK transcript.
Its loss-of-function section lists HIGH-impact variants in the reportable genes below `germline_lof_max_af`, with an
`NMD_ESCAPE` column (a nonsense or frameshift change in the last exon, or in the last 50 nt of the penultimate exon, of
a multi-exon transcript; ClinGen SVI PVS1); these calls are not LOFTEE-filtered and all need review.

**VAF and ORIGIN.** The last two columns give the sample's variant allele fraction (from FORMAT/AD) and an origin flag
that does not change the class. In germline mode, a VAF below `germline_low_vaf` (0.30) is flagged
POSSIBLE_MOSAIC_OR_CH, because a low-fraction variant in blood DNA may be mosaic or clonal haematopoiesis rather than
inherited (TP53 is the usual example); such calls need confirmation before being reported as germline. Sites with fewer than
`germline_mosaic_min_dp` (20) reads are UNRESOLVED, because a heterozygote can read below 0.30 by chance at low depth. Paired somatic
calls are flagged SOMATIC.
ClinVar significance (dated release; entries with conflicting interpretations or below `clinvar_min_stars` are not used)
is shown in its own column but is not used as evidence of oncogenicity, because ClinVar P/LP calls are germline
assertions; a ClinVar Benign or Likely benign call does mark a somatic variant LIKELY_BENIGN. ClinVar's separate
oncogenicity classification (somatic, ClinVar releases from 2024; INFO `ONC`/`ONCREVSTAT`) sets a lower bound when it
has at least `clinvar_min_stars` review stars: Oncogenic gives ONCOGENIC_STRONG, Likely oncogenic at least
ONCOGENIC_MODERATE, and Benign or Likely benign gives LIKELY_BENIGN. The value and its stars are in the `CLINVAR_ONC` and
`CLINVAR_ONC_STARS` columns of the pre-classification table, just before `VAF` and `ORIGIN`. In germline mode, a null variant is flagged
LOF_CANCER_GENE_REVIEW only in tumor suppressors marked as germline genes in the CGC (all CGC tumor suppressors if the
file has no `Germline` column), since loss of function is the predisposition mechanism there. Multi-allelic sites are
split before annotation so each alternate allele is classified on its own.

---

## Outputs

Results are grouped by stage under `outdir`; from variant calling on, and for the run manifest and MultiQC, by
analysis (`somatic/`, `germline/`, `rna_fusion/`). `<S>` is a sample ID (a tumor ID for somatic results). Folders of
optional steps appear only when the step runs.

```
<outdir>/
|-- 00_manifest/<analysis>/            preflight_references.txt (DNA arms), run_manifest.tsv (write_manifest)
|-- 01_mapping/rna_fusion/
|   |-- star/                          <S>.Log.final.out, .SJ.out.tab, .ReadsPerGene.out.tab, .Chimeric.out.junction
|   |-- salmon/                        <S>/quant.sf, Salmon logs and library type (salmon_transcripts)
|-- 02_markdup/                        <S>_dup_metrics.json (markdup), <S>_dup_metrics.txt (umi_markdup),
|                                      <S>_family_sizes.txt (umi_consensus)
|-- 03_metrics/
|   |-- mosdepth/                      <S>.mosdepth.summary.txt, .global.dist.txt, .region.dist.txt, .regions.bed.gz,
|   |                                  .thresholds.bed.gz, .quantized.bed.gz (TMB depth classes),
|   |                                  .gene_coverage.tsv (coverage_report_bed)
|   |-- picard/                        <S>.hs_metrics, .alignment_summary_metrics, .insert_size_metrics, .gc_bias.summary_metrics,
|   |                                  .bait_bias_summary_metrics, .pre_adapter_summary_metrics
|   |-- samtools/                      <S>.stats
|   |-- sex/                           <S>_sex.tsv
|   |-- verifybamid/                   <normal>.selfSM (verifybamid_db)
|   |-- somalier/                      somalier.pairs.tsv, somalier.samples.tsv, somalier.html, sample_swap_report.tsv,
|   |                                  somalier-ancestry.somalier-ancestry.tsv/.html (somalier_ancestry_labels)
|   |-- gate/                          <S>_qc_gate.tsv, cohort_qc_summary.tsv
|   |-- rna_fusion/rnaseq_qc/          <S>.rna_metrics, <S>.strandedness.txt
|-- 04_bqsr/bqsr_bams/                 <S>_bqsr.bam, .bai (bqsr)
|-- 05_variant_calling/
|   |-- somatic/
|   |   |-- snv/mutect2/               <S>_mutect2.vcf.gz
|   |   |-- snv/strelka2/              <S>_strelka2_snvs.vcf.gz, <S>_strelka2_indels.vcf.gz (+ .af.vcf.gz with AD/AF)
|   |   |-- snv/muse/                  <S>.MuSE.vcf.gz, <S>_muse.af.vcf.gz
|   |   |-- snv/varscan/               <S>_varscan_snv.vcf.gz, <S>_varscan_indel.vcf.gz (+ .af.vcf.gz)
|   |   |-- snv/callers_consensus/     <S>_consensus.vcf.gz, <S>_consensus_stats.txt, <S>_filtered.vcf.gz (final call set),
|   |   |                              <S>_filter_stats.txt
|   |   |-- sv/manta/                  <S>_manta_somaticSV.vcf.gz, <S>_manta_candidateSmallIndels.vcf.gz
|   |   |-- cnv/cnvkit/                <S>.cnr, <S>.cns, <S>-call.cns
|   |   |-- cnv/sequenza/              <S>_segments.txt, _mutations.txt, _confints_CP.txt, _alternative_solutions.txt,
|   |   |                              _fit_qc.tsv, _gene_loh.tsv, _arm_cn.tsv, _cn_summary.tsv, _sequenza_log.txt
|   |   |                              (includes the sex model used), plots (PDF), R objects
|   |   |-- msi/                       <S>_msi, _msi_all, _msi_dis, _msi_unstable, <S>_msi_status.txt, <S>_msi_genes.tsv
|   |   |-- contamination/             <S>.contamination.table, <S>.segments.table (pileup_variants)
|   |   |-- signatures/                <S>_signatures/ (mutational_signatures)
|   |   |-- evidence_bam/              <S>.evidence.bam, .bai (evidence_bed)
|   |-- germline/
|   |   |-- gvcf/                      <S>.g.vcf.gz (HaplotypeCaller)
|   |   |-- snv/                       <S>_filtered.vcf.gz (nv, vqsr) or <S>_genotyped.vcf.gz (none);
|   |   |                              snv/deepvariant/<S>.vcf.gz, <S>.g.vcf.gz with germline_caller: deepvariant
|   |   |-- joint/                     cohort_joint.vcf.gz (vqsr)
|   |   |-- sv/manta/                  <S>_manta_diploidSV.vcf.gz, <S>_manta_candidateSV.vcf.gz
|   |   |-- cnv/gcnv/                  <S>_gcnv_segments.vcf.gz, _gcnv_intervals.vcf.gz, _gcnv_denoised_cr.tsv (run_gcnv)
|   |   |-- roh/                       <S>_roh.txt (run_roh, WGS)
|   |   |-- hla/                       <S>_result.tsv, OptiType HLA class I types (run_optitype)
|   |   |-- evidence_bam/              <S>.evidence.bam, .bai (evidence_bed)
|   |-- rna_fusion/
|       |-- arriba/                    <S>_arriba_fusions.tsv, <S>_arriba_fusions.discarded.tsv
|       |-- starfusion/                <S>_starfusion.tsv
|       |-- fusioncatcher/             <S>_fusioncatcher.txt (use_fusioncatcher)
|       |-- fusioninspector/           <S>_fusioninspector.tsv (use_fusioninspector)
|       |-- consensus/                 <S>_fusion_consensus.tsv, <S>_fusion_consensus_stats.txt, <S>_fusion_candidates.txt,
|       |                              <S>_fusion_names.tsv (CTAT name -> reported name)
|       |-- trust4/                    <S>_report.tsv, <S>_airr.tsv (trust4_ref_dir)
|       |-- telescope/                 <S>-telescope_report.tsv (telescope_gtf)
|       |-- evidence_bam/              <S>.evidence.bam, .bai (evidence_bed)
|-- 06_annotation/<analysis>/          <S>_vep.vcf.gz, <S>_vep_summary.html; tmb/<S>_tmb.txt (somatic);
|                                      sv/<S>_<source>_annotsv.tsv (annotsv_dir); <S>_fusions_annotated.tsv (rna_fusion)
|-- 07_classification/<analysis>/      <S>_preclassification.tsv; maf/<S>.maf;
|                                      reports/<patient_id>_consolidated.tsv (somatic; a labelled tumor: <S>_consolidated.tsv),
|                                      <S>_germline_consolidated.tsv (germline), <patient_id>_fusion_report.tsv (rna_fusion)
|-- 0x_multiqc/<analysis>/             MultiQC reports per stage (01_trimming, 02_markdup, 03_metrics, 03_somalier,
|                                      06_annotation; one report for rna_fusion)
|-- 0x_crams/                          dedup/, bqsr/ (make_cram), rna/ (rna_cram_convert)
|-- pipeline_trace.txt, pipeline_report.html (pipeline_timeline.html, pipeline_dag.html with detailed_reports)
```

The deduplicated BAMs are kept in the work directory and are not copied to `outdir`; set `make_cram: true` to archive
them as CRAM (with `bqsr: true` the recalibrated BAMs are published in `04_bqsr/`). With `germline_on_normals: true` the
normals' germline results are written in the `somatic/` folders under the normal's sample ID.

---

## Troubleshooting

| Symptom | Cause | Solution |
|---|---|---|
| Exit 140 | LSF run-time limit | Raise `time` for that process in an extra config (see Running PAPIO) |
| Exit 255 in `PREPROCESS_AND_ALIGN_DNA` | fastp stopped on a damaged FASTQ (bwa-mem2 log ends normally) | `gzip -t` both FASTQs; download again (`prefetch` → `vdb-validate` → `fasterq-dump`) |
| `0 total baseline variants` in vcfeval | Truth FILTER is not PASS and `--all-records` was not used | Set FILTER to PASS in the truth VCF (Appendix B, step 11.1, `awk '{$7="PASS"}'`), or run vcfeval with --all-records |
| Many LOWSUPPORT calls | Tumor depth below `filter_min_tumor_dp` (20), or `call_regions_gz` narrower than the Mutect2 shards | Check coverage; use the padded `call_regions_gz` (Appendix A, step 4.4) |
| `command not found` (or `Rscript`/`awk` cannot open) for a `bin/` script such as `classify_variants.sh` | `nextflow/bin/` missing, scripts not executable (copied without permissions), or the pipeline folder not visible inside Singularity | Use `git clone` (keeps permissions) or `chmod +x nextflow/bin/*`; keep the pipeline under a path in `bind_paths` |
| `Unknown configuration profile: local` | `conf/` folder missing next to `nextflow.config` | Copy `conf/` from the release into `nextflow/` |
| Params ignored | Values in a `-c` config, or `${HOME}` in the params file | Put them in the params YAML with absolute paths |
| `Parameters changed since the run being resumed` | `-resume` after editing the params file | New `work_root`, or run without `-resume`; see Running PAPIO |
| `SEQUENZA ... Error is ignored` in the log; report shows `(no Sequenza output)` | Sequenza's model fit failed (too few heterozygous SNPs: very low purity or small panel; or no segment reaches `sequenza_segment_filter`, 3 Mb by default) | Expected for such samples; CNVkit calls copy number without purity. For a small targeted panel, lower `sequenza_segment_filter` (e.g. 500000). Check purity and coverage if it happens on a normal exome |
| Germline recall far below 0.9 | `germline_scatter_intervals` from another capture kit | Use the shards built for the sample's kit (Appendix A, step 4.3) |

---

## Appendix A. Reference resources (step by step)

Every step is marked **[Download]** or **[Generate]**. Output names match the parameter file in the "Parameter file" section.
Annotation versions: **DNA** - Ensembl 115 (VEP cache, gene BED, MANE CDS); **RNA** - GENCODE v49 (STAR index, Arriba,
RNA QC files, Salmon), the same gene set as Ensembl 115, and the CTAT library v44 (STAR-Fusion, FusionInspector).

### 0. Set-up

```bash
P=$PWD                                  # pipeline root, e.g. $HOME/genomics_pipeline
R=$P/references
E=$P/containers
REF=$R/genome/GRCh38.d1.vd1.fa
DICT=$R/genome/GRCh38.d1.vd1.dict
U=$R/downloads                          # raw downloads land here; delete them when Appendix A is done
SX="singularity exec -B $P"
BCF="$SX $E/bwa_samtools_fastp_fgbio_bcftools.sif"
GATK="$SX $E/gatk.sif"
mkdir -p $U $R/{genome/bwa_mem2,gatk,somalier,clinvar,gnomad,annotation,vep_cache,vep_plugins,cosmic,annotsv} \
         $R/{sequenza,msi,capture,gencode,star,rna_qc,ctat,arriba,fusion,fusioncatcher} \
         $R/{trust4,telescope,refphi,truth} $E
```

The folders and file names are those of `nextflow/setup.sh --full`, so the params files in `nextflow/params/` and the
ones written by setup work with either route:

```
references/
|-- downloads/              raw downloads
|-- genome/                 GRCh38.d1.vd1.fa, .fa.fai, .dict; bwa_mem2/ (index)
|-- gatk/                   GATK resource bundle VCFs restricted to the reference contigs (upstream file names)
|-- somalier/  clinvar/  gnomad/  vep_cache/  vep_plugins/  cosmic/  annotsv/
|-- annotation/             gene, MANE CDS, chromosome-arm and reportable-gene files, gCNV contig-ploidy priors
|-- sequenza/  msi/
|-- capture/<KIT>/          <KIT>_targets.bed, <KIT>_baits.bed, their .interval_list, <KIT>_targets_pad100.bed.gz,
|                           <KIT>_scatter_24_pad100/, <KIT>_scatter_150_pad100/, gcnv/; optional <KIT>_cnvkit_pon.cnn, <KIT>_inhouse_pon.vcf.gz
|-- gencode/  star/  rna_qc/  ctat/  arriba/  fusion/  fusioncatcher/  trust4/  telescope/  refphi/
`-- truth/                  benchmark truth sets (Appendix B)
```

Long steps (bwa-mem2 index, STAR index, VEP cache) need an interactive LSF session:

```bash
bsub -Is -q long -W 24:00 -M 64 -R "rusage[mem=64]" -n 16 /bin/bash
```

After opening the session, `cd` to the pipeline root and run the set-up block again (variables do not carry over).

---

### 1. Containers

All images live in `containers/`. Their file names must match the `container_*` parameters in `nextflow/nextflow.config`
(table in "Containers" above). Every image is also published at `oras://ghcr.io/chauhan2987/papio_containers:<tag>`:

```bash
for t in gatk bwa_samtools_fastp_fgbio_bcftools bcftool_cnvkit_vep_snpeff msisensorpro star_samtools; do
  singularity pull $E/$t.sif oras://ghcr.io/chauhan2987/papio_containers:$t
done
```

Download other images using the same approach. These images are:
annotsv, multiqc_mosdepth_verifybam2, sequenza, starfusion_fusioninspector_ctat, strelka2, strelka_tag_af, arriba, deepvar_freebayes, manta, muse_bcftools, sigprofiler, varscan, fusioncatcher, salmon, trust4, telescope, optitype_samtools, rtg_tools, somalier, vcf2maf

The images below can also be rebuilt from public sources.

#### 1.1 Single-tool images [Download]

```bash
cd $E
apptainer pull strelka2.sif      docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1
apptainer pull annotsv.sif       docker://quay.io/biocontainers/annotsv:3.5.3--py313hdfd78af_0
apptainer pull fusioncatcher.sif docker://clinicalgenomics/fusioncatcher:1.33
cd $P
```

#### 1.2 Multi-tool images built with Seqera Wave [Generate]

Install the Wave CLI once:

```bash
URL=$(curl -fsSL https://api.github.com/repos/seqeralabs/wave-cli/releases/latest | grep "browser_download_url.*linux-x86_64" | cut -d '"' -f 4)
curl -fsSL "$URL" -o $P/wave && chmod +x $P/wave
```

For each image below: write its Dockerfile, build it, then pull the URL that Wave prints.
Each build gives a new URL, so always paste the one printed by your own build.

```bash
# after writing the Dockerfile of one image:
NAME=sequenza                                     # or varscan, strelka_tag_af, vcf2maf, arriba
$P/wave --containerfile Dockerfile --await        # prints docker://wave.seqera.io/wt/.../wave/build:...
WAVE_URL="docker://wave.seqera.io/wt/xxxxxxxxxxxx/wave/build:xxxxxxxxxxxxxxxx"   # paste the printed URL
singularity pull $E/${NAME}.sif "$WAVE_URL"
```

**Sequenza** (`sequenza.sif`)

```bash
cat > Dockerfile <<'EOF'
FROM condaforge/mambaforge:latest
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10.0 r-sequenza=3.0.0 sequenza-utils=3.0.0 \
    bioconductor-copynumber bioconductor-genomeinfodb r-biocmanager r-remotes samtools htslib \
    && mamba clean -afy
RUN R -e 'BiocManager::install("GenomeInfoDbData", update=FALSE, ask=FALSE)'
# Bioconductor copynumber rejects assembly="hg38"; replace it with the hg38-capable fork (pinned commit)
RUN R -e 'remotes::install_github("ShixiangWang/copynumber@2e31d59f558d5e3e44215da2f612707791e3f030", dependencies=FALSE, upgrade="never")'
# stop the build if hg38 is still not supported
RUN Rscript -e 'library(copynumber); p<-sort(runif(200,1e6,2e8)); r<-aspcf(logR=data.frame(chrom=1,pos=p,s=rnorm(200)), BAF=data.frame(chrom=1,pos=p,s=runif(200)), assembly="hg38", verbose=FALSE); cat("copynumber hg38 OK\n")'
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8
EOF
```

**VarScan** (`varscan.sif`)

```bash
cat > Dockerfile <<'EOF'
FROM condaforge/mambaforge:latest
RUN mamba create -y -n varscan -c conda-forge -c bioconda \
    varscan=2.4.6 samtools=1.23.1 bcftools=1.23.1 \
    && mamba clean -afy
ENV PATH="/opt/conda/envs/varscan/bin:${PATH}"
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8
EOF
```

**Strelka2 allele-fraction tagging** (`strelka_tag_af.sif`)

```bash
cat > Dockerfile <<'EOF'
FROM condaforge/mambaforge:latest
RUN mamba create -y -n strelka_af -c conda-forge -c bioconda \
    python=3.11 pysam=0.22 samtools=1.23.1 bcftools=1.23.1 htslib=1.23.1 \
    && mamba clean -afy
ENV PATH="/opt/conda/envs/strelka_af/bin:${PATH}"
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8
EOF
```

**vcf2maf** (`vcf2maf.sif`)

```bash
cat > Dockerfile <<'EOF'
FROM condaforge/mambaforge:latest
RUN mamba create -y -n vcf2maf -c conda-forge -c bioconda \
    vcf2maf=1.6.22 samtools=1.23.1 htslib=1.23.1 \
    && mamba clean -afy
ENV PATH="/opt/conda/envs/vcf2maf/bin:${PATH}"
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8
EOF
```

**Arriba** (`arriba.sif`)

```bash
cat > Dockerfile <<'EOF'
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
      wget ca-certificates zlib1g libbz2-1.0 liblzma5 libcurl4 libgomp1 && \
    wget -qO- https://github.com/suhrig/arriba/releases/download/v2.5.1/arriba_v2.5.1.tar.gz | tar xz && \
    install -m755 arriba_v2.5.1/arriba /usr/local/bin/arriba && \
    rm -rf arriba_v2.5.1 && apt-get purge -y wget && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8
EOF
```

---

### 2. Reference genome (GDC GRCh38.d1.vd1)

#### 2.1 FASTA [Download]

```bash
curl -L -o $U/GRCh38.d1.vd1.fa.tar.gz https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834
tar -xzf $U/GRCh38.d1.vd1.fa.tar.gz -C $R/genome
ls $REF
```

#### 2.2 FASTA index, bwa-mem2 index, sequence dictionary [Generate]

```bash
$BCF bash -c "samtools faidx $REF && bwa-mem2 index -p $R/genome/bwa_mem2/GRCh38.d1.vd1.fa $REF"
$GATK gatk CreateSequenceDictionary -R $REF -O $DICT
```

---

### 3. Known-site and population VCFs

#### 3.1 GATK resource bundle [Download]

```bash
cd $U
B=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
S=https://storage.googleapis.com/gatk-best-practices/somatic-hg38
wget -c $B/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c $B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c $B/Homo_sapiens_assembly38.known_indels.vcf.gz
wget -c $B/hapmap_3.3.hg38.vcf.gz
wget -c $B/1000G_omni2.5.hg38.vcf.gz
wget -c $B/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c $S/af-only-gnomad.hg38.vcf.gz
wget -c $S/small_exac_common_3.hg38.vcf.gz
wget -c $S/1000g_pon.hg38.vcf.gz
$BCF bgzip -f $U/Homo_sapiens_assembly38.dbsnp138.vcf          # dbSNP is shipped uncompressed
cd $P
```

#### 3.2 Somalier sites and ancestry files [Download]

```bash
wget -c -O $U/sites.hg38.vcf.gz https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz   # restricted in 3.3

wget -O $R/somalier/ancestry-labels-1kg.tsv https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv
wget -O $U/1kg.somalier.tar.gz "https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1"
tar -xzf $U/1kg.somalier.tar.gz -C $R/somalier
```

#### 3.3 Restrict every VCF to the d1.vd1 contigs [Generate]

```bash
$BCF bash -c '
set -euo pipefail
U='"$U"'; R='"$R"'
CONTIGS=$(grep "^@SQ" '"$DICT"' | sed "s/.*SN:\([^\t]*\).*/\1/" | paste -sd,)
reconcile() {
  [[ -f "$1" ]] || { echo "SKIP (missing): $1"; return 0; }
  echo ">>> $1 -> $2"
  bcftools view -t "$CONTIGS" "$1" -Oz -o "$2" && tabix -f -p vcf "$2"
}
reconcile $U/Homo_sapiens_assembly38.dbsnp138.vcf.gz          $R/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz
reconcile $U/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz $R/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
reconcile $U/Homo_sapiens_assembly38.known_indels.vcf.gz      $R/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz
reconcile $U/af-only-gnomad.hg38.vcf.gz                       $R/gatk/af-only-gnomad.hg38.vcf.gz
reconcile $U/small_exac_common_3.hg38.vcf.gz                  $R/gatk/small_exac_common_3.hg38.vcf.gz
reconcile $U/1000g_pon.hg38.vcf.gz                            $R/gatk/1000g_pon.hg38.vcf.gz
reconcile $U/hapmap_3.3.hg38.vcf.gz                           $R/gatk/hapmap_3.3.hg38.vcf.gz
reconcile $U/1000G_omni2.5.hg38.vcf.gz                        $R/gatk/1000G_omni2.5.hg38.vcf.gz
reconcile $U/1000G_phase1.snps.high_confidence.hg38.vcf.gz    $R/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz
reconcile $U/sites.hg38.vcf.gz                                $R/somalier/sites.hg38.vcf.gz
'
```

Check one file: the record count drops only slightly, and no contig outside the dictionary remains.

```bash
F=$R/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
zcat $U/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | grep -vc '^#'
zcat $F | grep -vc '^#'
comm -23 <(zcat $F | grep -v '^#' | cut -f1 | sort -u) <(grep '^@SQ' $DICT | sed 's/.*SN:\([^\t]*\).*/\1/' | sort)   # must print nothing
$GATK gatk ValidateVariants -V $F -R $REF --validation-type-to-exclude ALL
```

---

### 4. Capture intervals

Design BEDs come from Agilent SureDesign (login required): https://earray.chem.agilent.com/suredesign/. Put the two
design BEDs of each kit in `references/capture/<KIT>/`, named as downloaded.

| Kit | Design | Design BEDs (targets, probes) | Folder |
|---|---|---|---|
| Demonstration (SEQC2, HG002): SureSelect Human All Exon V6 + UTR | S07604624 | `S07604624_Regions.bed`, `S07604624_Covered.bed` | `references/capture/S07604624/` |
| Own data (RCC cohort): SureSelect Human All Exon V7 | S31285117 | `S31285117_Regions.bed`, `S31285117_MergedProbes.bed` | `references/capture/S31285117/` |

Check that a BED is on hg38 (EGFR exons lie between chr7:55.0 and 55.2 Mb):

```bash
cut -f1 $R/capture/S07604624/S07604624_Regions.bed | sort -u | head
grep -w chr7 $R/capture/S07604624/S07604624_Regions.bed | awk '$2>55019000 && $3<55211000' | head -3
```

Sections 4.1-4.4 and 5.2 are run once per kit; they write the same files as
`setup.sh --full <dir> --kit <ID> --kit-regions <targets BED> --kit-probes <probes BED>`. The demonstration (SEQC2,
HG002) uses the V6 + UTR kit; for your own data, comment the V6 line and uncomment the V7 line, then run 4.1-4.4 again.

```bash
ID=S07604624; REGIONS=S07604624_Regions.bed; PROBES=S07604624_Covered.bed        # demonstration: V6 + UTR
# ID=S31285117; REGIONS=S31285117_Regions.bed; PROBES=S31285117_MergedProbes.bed # own data: V7
K=$R/capture/$ID
```

#### 4.1 Clean, sorted target and bait BEDs [Generate]

Vendor BEDs start with `track`/`browser` lines; the clean copies are used everywhere else (`interval_bed` is
`<KIT>_targets.bed`).

```bash
clean_bed() { awk 'BEGIN{OFS="\t"} !/^#/ && !/^track/ && !/^browser/ && NF>=3 {print $1,$2,$3}' "$1" | LC_ALL=C sort -k1,1 -k2,2n > "$2"; }
clean_bed $K/$REGIONS $K/${ID}_targets.bed
clean_bed $K/$PROBES  $K/${ID}_baits.bed
head -3 $K/${ID}_targets.bed
```

#### 4.2 Bait and target interval lists [Generate]

```bash
$GATK gatk BedToIntervalList -I $K/${ID}_baits.bed   -SD $DICT -O $K/${ID}_baits.interval_list
$GATK gatk BedToIntervalList -I $K/${ID}_targets.bed -SD $DICT -O $K/${ID}_targets.interval_list
```

#### 4.3 Padded calling shards: somatic (24) and germline (150) [Generate]

Targets ± 100 bp; overlapping padded targets are merged; no target is split between shards.

```bash
for n in 24 150; do
  $GATK gatk SplitIntervals -R $REF -L $K/${ID}_targets.interval_list --interval-padding 100 --interval-merging-rule ALL \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION --scatter-count $n -O $K/${ID}_scatter_${n}_pad100
done
```

#### 4.4 Padded call regions for Strelka2 and Manta (`call_regions_gz`) [Generate]

Built from the somatic shards, so all callers see the same territory.

```bash
grep -hv '^@' $K/${ID}_scatter_24_pad100/*.interval_list \
  | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' | LC_ALL=C sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"} {if($1==c && $2<=e){if($3>e)e=$3} else {if(c!="")print c,s,e; c=$1;s=$2;e=$3}} END{if(c!="")print c,s,e}' \
  > $K/${ID}_targets_pad100.bed
$BCF bash -c "bgzip -f $K/${ID}_targets_pad100.bed && tabix -f -p bed $K/${ID}_targets_pad100.bed.gz"
```

Params YAML entries (each kit's files carry its ID and sit in its own folder):

| Parameter | File under `references/capture/<KIT>/` |
|---|---|
| `interval_bed` | `<KIT>_targets.bed` |
| `call_regions_gz` | `<KIT>_targets_pad100.bed.gz` |
| `bait_interval_list` | `<KIT>_baits.interval_list` |
| `target_interval_list` | `<KIT>_targets.interval_list` |
| `scatter_intervals` | `<KIT>_scatter_24_pad100` |
| `germline_scatter_intervals` | `<KIT>_scatter_150_pad100` |


#### 4.5 In-house Mutect2 panel of normals (optional) [Generate]

Run this after steps 4.1-4.4 for the kit, once its normals have been through the pipeline (it uses their final
BAMs). The public `1000g_pon` (step 3.3) is the panel of normals used throughout this protocol. A panel built from your own
normals on the same capture kit can also remove recurrent artifacts specific to that kit. `main_pon.nf` builds a panel
from final normal BAMs: it runs Mutect2 tumor-only on each normal over `scatter_intervals`, then
GenomicsDBImport and CreateSomaticPanelOfNormals with `germline_resource`.

Samplesheet (normal BAMs only):

```
sample_id,bam,bai
S001_N,/path/S001_N.bam,/path/S001_N.bam.bai
S002_N,/path/S002_N.bam,/path/S002_N.bam.bai
```

`params_pon.yaml` (same conventions as `params_demo.yaml` in "Parameter file"; V7 shards shown):

```yaml
pon_samplesheet:    /home/<user>/normals_pon.csv
container_dir:      ${projectDir}/../containers
work_root:          /scratch/<user>/work_pon
outdir:             /home/<user>/Result_pon
bind_paths:         [/home/<user>, /scratch/<user>]
ref_fasta:          ${projectDir}/../references/genome/GRCh38.d1.vd1.fa
ref_fai:            ${projectDir}/../references/genome/GRCh38.d1.vd1.fa.fai
ref_dict:           ${projectDir}/../references/genome/GRCh38.d1.vd1.dict
scatter_intervals:  ${projectDir}/../references/capture/S31285117/S31285117_scatter_24_pad100
germline_resource:  ${projectDir}/../references/gatk/af-only-gnomad.hg38.vcf.gz
lsf_queue_short:    short
lsf_queue_medium:   medium
lsf_queue_high:     long
lsf_queue_long:     long
lsf_project:        null
```

```bash
nextflow run ${HOME}/genomics_pipeline/nextflow/main_pon.nf -c ${HOME}/genomics_pipeline/nextflow/pon.config \
 -params-file params_pon.yaml -profile production
```

- The panel is written to `<outdir>/00_reference/final_pon/panel_of_normals.vcf.gz` (with `.tbi`). Copy both to
  `references/capture/<KIT>/` as `<KIT>_inhouse_pon.vcf.gz` and `<KIT>_inhouse_pon.vcf.gz.tbi`, then set `panel_of_normals` to that file.
- This is a Mutect2 (somatic) panel of normals like `1000g_pon`, but it is valid only for the kit it was built on,
  so it is kept with that kit's files.
- The workflow is optional, is not part of this protocol, and was not used for any result in the paper.
- The panel and the per-normal VCFs in `<outdir>/00_reference/pon/per_normal/` are derived from patients' germline
  data: do not commit or share them without approval.

---

### 5. Copy-number and MSI resources

#### 5.1 Sequenza and MSI list [Generate]

```bash
$SX $E/sequenza.sif sequenza-utils gc_wiggle --fasta $REF -w 50 -o $R/sequenza/GRCh38.d1.vd1.gc50Base.wig.gz
$SX $E/msisensorpro.sif msisensor-pro scan -d $REF -o $R/msi/GRCh38.d1.vd1.msisensor.list
```

#### 5.2 Germline copy-number resources (only with `run_gcnv: true` and ≥ 30 normals) [Generate]

```bash
G=$K/gcnv              # K and ID from section 4 (demonstration data has < 30 normals, so gCNV is not needed)
mkdir -p $G/${ID}_scatter_24_pad250

$GATK gatk PreprocessIntervals -R $REF -L $K/${ID}_targets.interval_list \
  --bin-length 0 --padding 250 --interval-merging-rule OVERLAPPING_ONLY -O $G/${ID}_preprocessed_intervals.interval_list
$GATK gatk AnnotateIntervals -R $REF -L $G/${ID}_preprocessed_intervals.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY -O $G/${ID}_annotated_intervals.tsv
$GATK gatk IntervalListTools --INPUT $G/${ID}_preprocessed_intervals.interval_list \
  --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_COUNT 24 --OUTPUT $U/${ID}_gcnv_scatter_tmp
i=0; for d in $U/${ID}_gcnv_scatter_tmp/temp_*; do i=$((i+1)); printf -v n '%04d' $i; cp "$d/scattered.interval_list" "$G/${ID}_scatter_24_pad250/shard_${n}.interval_list"; done
rm -rf $U/${ID}_gcnv_scatter_tmp
```

Contig-ploidy priors (expected copies of each chromosome) depend on the genome, not the kit; written once for all kits:

```bash
awk 'BEGIN{OFS="\t"; print "CONTIG_NAME","PLOIDY_PRIOR_0","PLOIDY_PRIOR_1","PLOIDY_PRIOR_2","PLOIDY_PRIOR_3"
  for (i=1;i<=22;i++) print "chr"i,"0.01","0.01","0.97","0.01"
  print "chrX","0.01","0.49","0.49","0.01"; print "chrY","0.495","0.495","0.01","0.01"}' > $R/annotation/GRCh38_contig_ploidy_priors.tsv
```

CNVkit runs against the matched normal (`cnvkit_pon: null`); a panel of normals built later for a kit goes to
`references/capture/<KIT>/<KIT>_cnvkit_pon.cnn`.

---

### 6. Gene annotation BEDs (Ensembl 115)

#### 6.1 Ensembl GTF [Download]

The same release (115) as the VEP cache in section 7.1.

```bash
wget -c -P $U https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
GTF=$U/Homo_sapiens.GRCh38.115.gtf.gz
```

#### 6.2 Gene BED, TMB territory, reportable genes [Generate]

```bash
# gene_bed: chr, start (0-based), end, symbol, Ensembl ID; chr names; chr1-22, X, Y, M only
zcat $GTF | awk -F'\t' '$3=="gene" && $1 ~ /^([0-9]+|X|Y|MT)$/ {
    id=""; sym=""; chr=($1=="MT") ? "chrM" : "chr" $1
    if (match($9, /gene_id "[^"]+"/))   id=substr($9,RSTART+9,RLENGTH-10)
    if (match($9, /gene_name "[^"]+"/)) sym=substr($9,RSTART+11,RLENGTH-12)
    if (sym=="") sym=id
    printf "%s\t%d\t%d\t%s\t%s\n", chr, $4-1, $5, sym, id }' \
  | LC_ALL=C sort -k1,1 -k2,2n > $R/annotation/ensembl115_genes.bed

# merge overlapping intervals of a sorted BED
merge3() { awk 'BEGIN{OFS="\t"} {if($1==c && $2<=e){if($3>e)e=$3} else {if(c!="")print c,s,e; c=$1;s=$2;e=$3}} END{if(c!="")print c,s,e}'; }

# tmb_coding_bed: MANE Select coding exons, merged
zcat $GTF | awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS" && $9 ~ /tag "MANE_Select"/ {
    c=$1; if (c=="MT") c="M"; if (c !~ /^([0-9]+|X|Y|M)$/) next
    print "chr" c, $4-1, $5 }' | LC_ALL=C sort -k1,1 -k2,2n | merge3 > $R/annotation/ensembl115_MANE_CDS.bed

# germline_report_genes: your own list, one symbol per line
cat > $R/annotation/reportable_genes.txt <<'EOF'
VHL
PBRM1
SETD2
BAP1
KDM5C
SMARCB1
FH
FLCN
TSC1
TSC2
MTOR
PTEN
TP53
CDKN2A
NF2
MET
SDHA
SDHB
SDHC
SDHD
TFE3
TFEB
EOF

# coverage_report_bed: MANE coding exons of the reportable genes (column 4 = gene)
zcat $GTF | awk -F'\t' -v list=$R/annotation/reportable_genes.txt 'BEGIN{OFS="\t"; while ((getline g < list) > 0) G[toupper(g)]=1}
  $3=="CDS" && $9 ~ /tag "MANE_Select"/ && match($9, /gene_name "[^"]+"/) {
    g=substr($9,RSTART+11,RLENGTH-12); c=$1; if (c=="MT") c="M"
    if (c ~ /^([0-9]+|X|Y|M)$/ && (toupper(g) in G)) print "chr" c, $4-1, $5, g }' \
  | LC_ALL=C sort -k1,1 -k2,2n -u > $R/annotation/ensembl115_reportable_genes_MANE_CDS.bed
```

#### 6.3 Chromosome-arm BED [Download + Generate]

```bash
wget -c -O $R/downloads/cytoBand.hg38.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
zcat $R/downloads/cytoBand.hg38.txt.gz | awk -F'\t' 'BEGIN{OFS="\t"}
  $1 ~ /^chr([0-9]+|X|Y)$/ { k=$1 SUBSEP substr($4,1,1)
    if (!(k in s) || $2<s[k]) s[k]=$2
    if (!(k in e) || $3>e[k]) e[k]=$3 }
  END { for (k in s) { split(k,a,SUBSEP); arm=substr(a[1],4) a[2]
          if (arm ~ /^(13|14|15|21|22)p$/) continue      # acrocentric p arms: no assembled sequence
          print a[1], s[k], e[k], arm } }' \
  | LC_ALL=C sort -k1,1 -k2,2n > $R/annotation/GRCh38_arms.bed
```

---

### 7. Variant annotation databases

#### 7.1 VEP cache [Download]

The cache release must match the VEP version in the container (115).

```bash
$SX $E/bcftool_cnvkit_vep_snpeff.sif vep --help 2>&1 | grep -i "ensembl-vep" | head -2
curl -L -o $U/homo_sapiens_vep_115_GRCh38.tar.gz \
  https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar -xzf $U/homo_sapiens_vep_115_GRCh38.tar.gz -C $R/vep_cache
```

#### 7.2 Chromosome-name map for ClinVar and COSMIC [Generate]

```bash
for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' $c $c; done > $R/clinvar/chr_rename.txt
printf 'MT\tchrM\n' >> $R/clinvar/chr_rename.txt
```

#### 7.3 ClinVar [Download + Generate]

The release date (`##fileDate`) goes into the file name, so the params file records which ClinVar was used.

```bash
wget -c -P $U https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
CV_DATE=$(zcat $U/clinvar.vcf.gz | awk -F= '/^##fileDate=/ {gsub(/-/,"",$2); print $2; exit}')
$BCF bash -c "set -euo pipefail
  bcftools annotate --rename-chrs $R/clinvar/chr_rename.txt $U/clinvar.vcf.gz -Oz -o $R/clinvar/clinvar_${CV_DATE}.vcf.gz
  tabix -f -p vcf $R/clinvar/clinvar_${CV_DATE}.vcf.gz"
ls $R/clinvar/        # clinvar_<YYYYMMDD>.vcf.gz: put this name into clinvar_vcf
```

#### 7.4 COSMIC v104 (license required) [Download + Generate]

COSMIC's scripted download takes two calls: the first, with your COSMIC login, returns a temporary link;
the second downloads the file. File paths are listed at https://cancer.sanger.ac.uk/cosmic/download (GRCh38, v104).

```bash
cd $U
read -p "COSMIC e-mail: " CE; read -s -p "COSMIC password: " CP; echo
AUTH=$(printf '%s:%s' "$CE" "$CP" | base64)
get_cosmic() {   # get_cosmic <path in the COSMIC bucket> <output file>
  link=$(curl -s -H "Authorization: Basic $AUTH" \
         "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=$1&bucket=downloads" \
         | sed 's/.*"url" *: *"\([^"]*\)".*/\1/')
  curl -L -o "$2" "$link"
}
get_cosmic grch38/cosmic/v104/Cosmic_CancerGeneCensus_Tsv_v104_GRCh38.tar          Cosmic_CancerGeneCensus_Tsv_v104_GRCh38.tar
get_cosmic grch38/cosmic/v104/Cosmic_GenomeScreensMutant_Normal_Vcf_v104_GRCh38.tar Cosmic_GenomeScreensMutant_Normal_Vcf_v104_GRCh38.tar
tar -xf Cosmic_CancerGeneCensus_Tsv_v104_GRCh38.tar -C $R/cosmic Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz
tar -xf Cosmic_GenomeScreensMutant_Normal_Vcf_v104_GRCh38.tar  # → Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.vcf.gz
# chr names, sorted, transcript-level duplicate rows collapsed to one record per variant
$BCF bash -c "set -euo pipefail
  bcftools annotate --rename-chrs $R/clinvar/chr_rename.txt Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.vcf.gz -Ou \
    | bcftools sort -T $U -Ou | bcftools norm -d exact -Oz -o $R/cosmic/Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.chr_dedup.vcf.gz
  tabix -f -p vcf $R/cosmic/Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.chr_dedup.vcf.gz"
cd $P
```

#### 7.5 AnnotSV annotations [Download]

```bash
wget -c -P $U https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz
tar -xzf $U/Annotations_Human_3.5.tar.gz -C $R/annotsv
ls $R/annotsv
```

#### 7.6 gnomAD v4.1 filtering allele frequency (`gnomad_faf_vcf`) [Download]

The population frequency used for benign calls (BA1, BS1) in variant classification is the gnomAD v4.1 joint
(exome + genome) grpmax filtering allele frequency: INFO field `fafmax_faf95_max_joint` in gnomAD's own VCFs,
`fafmax_faf95_max_j` in the file below (`_e` exome, `_g` genome, `_j` joint), the name the pipeline reads. A size-reduced copy of the
gnomAD v4.1 joint sites VCF, restricted to MANE transcript regions and to the population-frequency fields (including
the filtering allele frequencies) and prepared for VEP, is published by the Catalan Institute of Oncology on Zenodo
(DOI 10.5281/zenodo.23531816, CC BY 4.0; 33.6 GB). Check that the header lists `fafmax_faf95_max_j`.

```bash
cd $R/gnomad
wget -c https://zenodo.org/records/23531816/files/variantsGnomADv4.1Joint.vcf.bgz
wget -c https://zenodo.org/records/23531816/files/variantsGnomADv4.1Joint.vcf.bgz.tbi
md5sum variantsGnomADv4.1Joint.vcf.bgz      # 92f6f7c9f7f375ad64b171de6fa7d2b0
$BCF bcftools view -h variantsGnomADv4.1Joint.vcf.bgz | grep '^##INFO=<ID=fafmax_faf95_max_j,'
cd $P
```

Variants outside MANE transcript regions get no gnomAD frequency from this file. To cover every PASS site instead,
build the file from the official per-chromosome gnomAD files with `nextflow/scripts/prepare_gnomad_faf.sh <out_dir> [threads]`
(about 900 GB downloaded one chromosome at a time; output `gnomad.joint.v4.1.faf95.vcf.gz`, with the field renamed
`fafmax_faf95_max_j`).

#### 7.7 VEP NMD plugin (`vep_plugins_dir`) [Download]

Marks truncating variants that escape nonsense-mediated decay (last 50 nucleotides of the penultimate exon); used by
variant classification. The plugin needs no data files.

```bash
wget -O $R/vep_plugins/NMD.pm https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/115/NMD.pm
```

---

### 8. RNA fusion resources

#### 8.1 GENCODE v49 annotation [Download]

GENCODE v49 is the gene set of Ensembl 115, so the RNA arm (STAR, Arriba, RNA QC, Salmon) and the DNA arm use the same
gene models. STAR-Fusion and FusionInspector use the gene models of the CTAT library (8.2, GENCODE v44). The fusion
consensus matches callers on Ensembl gene IDs, so a gene renamed between v44 and v49 is still one fusion; fusions are
reported with the v49 symbols, and FusionInspector receives the CTAT (v44) symbols it looks up.

```bash
wget -c -P $U https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
zcat $U/gencode.v49.annotation.gtf.gz > $R/gencode/gencode.v49.annotation.gtf    # the params YAML uses the uncompressed .gtf
RGTF=$R/gencode/gencode.v49.annotation.gtf
```

#### 8.2 CTAT genome library, GENCODE v44 [Download]

Check the file name in https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
```bash
wget -c -P $U https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
tar -xzf $U/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz -C $R/ctat
ls $R/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir | head
```

#### 8.3 Arriba database files [Download + Generate]

```bash
wget -c -P $U https://github.com/suhrig/arriba/releases/download/v2.5.1/arriba_v2.5.1.tar.gz
tar -xzf $U/arriba_v2.5.1.tar.gz -C $U
cp $U/arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz \
   $U/arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
   $U/arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3 $R/arriba/
# known_fusions_tsv: uncompressed copy for the fusion annotation step
zcat $R/arriba/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz > $R/fusion/known_fusions_hg38_GRCh38_v2.5.1.tsv
```

#### 8.4 FusionCatcher database (only with `use_fusioncatcher: true`) [Download]

```bash
cd $R/fusioncatcher
for p in aa ab ac ad; do wget -c -O human_v102.tar.gz.$p https://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.$p/download; done
cat human_v102.tar.gz.* > human_v102.tar.gz && rm -f human_v102.tar.gz.*
tar -xzf human_v102.tar.gz && ls human_v102 | head
cd $P
```

#### 8.5 STAR index on the DNA reference genome [Generate]

Built on the same FASTA as the DNA arm, with GENCODE v49. Use a new folder so an older index is not overwritten.
`--sjdbOverhang` is read length − 1: 149 for 2 × 150 bp reads, 99 for 2 × 100 bp; 100 works well for both.

```bash
$SX $E/star_samtools.sif STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir $R/star/StarIndex_gencode_v49 --genomeFastaFiles $REF --sjdbGTFfile $RGTF --sjdbOverhang 100 \
  --limitGenomeGenerateRAM 60000000000
```

#### 8.6 refFlat and rRNA interval list for Picard CollectRnaSeqMetrics [Generate]

```bash
cat > $P/make_refflat.py <<'EOF'
# GENCODE GTF → Picard refFlat; exons written in ascending order (required by Picard, also for minus-strand genes)
import re, sys
from collections import defaultdict
gtf, out_path = sys.argv[1], sys.argv[2]
TID = re.compile(r'transcript_id "([^"]+)"'); GNM = re.compile(r'gene_name "([^"]+)"')
gene, chrom, strand, exons, cds = {}, {}, {}, defaultdict(list), {}
for line in open(gtf):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('exon', 'CDS'): continue
    tid = TID.search(f[8]).group(1); s, e = int(f[3]) - 1, int(f[4])
    if f[2] == 'exon':
        g = GNM.search(f[8]); gene[tid] = g.group(1) if g else tid
        chrom[tid], strand[tid] = f[0], f[6]; exons[tid].append((s, e))
    else:
        c0, c1 = cds.get(tid, (s, e)); cds[tid] = (min(c0, s), max(c1, e))
with open(out_path, 'w') as out:
    for tid in gene:
        ex = sorted(exons[tid]); tx0, tx1 = ex[0][0], ex[-1][1]
        c0, c1 = cds.get(tid, (tx1, tx1))                  # non-coding: CDS start = end = tx end
        out.write('\t'.join([gene[tid], tid, chrom[tid], strand[tid], str(tx0), str(tx1), str(c0), str(c1),
                             str(len(ex)), ','.join(str(a) for a, _ in ex) + ',', ','.join(str(b) for _, b in ex) + ',']) + '\n')
EOF
$SX $E/bcftool_cnvkit_vep_snpeff.sif python3 $P/make_refflat.py $RGTF $R/rna_qc/gencode.v49.refFlat

# rRNA interval list = sequence dictionary header + rRNA / Mt_rRNA transcripts
{ grep '^@' $DICT
  awk -F'\t' '$3=="transcript" && $9 ~ /gene_type "(Mt_)?rRNA(_pseudogene)?"/ {
      match($9, /transcript_id "[^"]+"/); print $1"\t"$4"\t"$5"\t"$7"\t"substr($9,RSTART+15,RLENGTH-16) }' $RGTF
} > $R/rna_qc/gencode.v49.rrna.interval_list
```

#### 8.7 Optional Salmon quantification of RNA: transcript FASTA [Generate]

Salmon (container `salmon.sif`, Salmon 1.10.3) quantifies STAR's transcriptome alignments
(`Aligned.toTranscriptome.out.bam`) in alignment mode, as in nf-core/rnaseq: one pass over the reads and no Salmon
index. It needs the transcript sequences of the GTF the STAR index was built with (step 8.5): every transcript in STAR's
BAM must be in this FASTA, or Salmon stops. `nextflow/scripts/gtf_to_transcripts.sh` joins the exons of each transcript
(minus strand reverse-complemented; names = `transcript_id`):

```bash
cp $P/nextflow/scripts/gtf_to_transcripts.sh $U/
TMPDIR=$U $BCF bash $U/gtf_to_transcripts.sh --gtf $R/gencode/gencode.v49.annotation.gtf --fasta $REF \
    --out $R/star/gencode.v49.annotation.transcripts.fa
grep -c '^>' $R/star/gencode.v49.annotation.transcripts.fa        # one sequence per GENCODE v49 transcript
```

Set `salmon_transcripts` to this file. With `rna_use_ctat_reference: true` (STAR on the CTAT library), build it from
`ctat_genome_lib_build_dir/ref_annot.gtf` and `ref_genome.fa` instead. `check_installation.sh` tests that the FASTA holds
every transcript of the GTF in use.

#### 8.8 TRUST4 TCR/BCR references [Download]

```bash
curl -fsSL --retry 2 https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/hg38_bcrtcr.fa -o $R/trust4/hg38_bcrtcr.fa

curl -fsSL --retry 2 https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/human_IMGT+C.fa -o $R/trust4/human_IMGT+C.fa
```

#### 8.9 Human endogenous retrovirus (HERV) annotation [Download]

```bash
curl -fsSL --retry 2 https://media.githubusercontent.com/media/mlbendall/telescope_annotation_db/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf -o $R/telescope/HERV_rmsk.hg38.v2.gtf
```

#### 8.10 refPhi annotation [Download]

```bash
curl -fsSL --retry 2 https://raw.githubusercontent.com/humengying0907/InstaPrismSourceCode/main/refPhi/RCC_refPhi.RDS -o $R/refphi/RCC_refPhi.RDS
```
---

## Appendix B. Demonstration data and benchmarking

| Arm | Data | Truth set | Comparison |
|---|---|---|---|
| Somatic | SEQC2 HCC1395 (SRR7890883) / HCC1395BL (SRR7890874), SureSelect V6 + UTR; SRA SRP162370 | SEQC2 v1.2.1 high-confidence SNVs/indels and regions (step 9.2) | `rtg vcfeval --sample ALT` on site-only VCFs, SNVs and indels separately; calls at LowConf/Unclassified superset sites and novel calls < 5% VAF excluded |
| Germline | GIAB HG002, Novogene SureSelect V6 exome (github.com/novogene-europe/ngs_demo_datasets) | GIAB v4.2.1 (step 9.3) | genotype-aware `rtg vcfeval` in callable targets ∩ benchmark regions |
| RNA fusion | HCC1395 bulk RNA-seq SRR9134735 (SRA SRP199641) | Fusion genes of the RNA-supported SVs of Talsania et al., 2022 (Additional file 9; 24 entries, 23 distinct events) | `fusion_truth_match.py` (gene pairs in either orientation, HGNC renaming, intergenic gene lists) |

Helper scripts: `make_figures.py` (Figures 2 and 4), `benchmark_plots.py` (Figures 3 and 5), `trace_table.py` (Table 1)
and `fusion_truth_match.py` (step 13); commands in step 14. Table 2 is `somatic_metrics.tsv` (steps 11.5-11.6).

---

### 9. Benchmark truth sets

#### 9.0 Shell variables (set once)
```bash
PIPE=${HOME}/genomics_pipeline
DATA=$PIPE/runs                          # outdirs of the demonstration runs (nextflow/params/*.yaml)
R=$PIPE/references
E=$PIPE/containers
T=$R/truth
REF=$R/genome/GRCh38.d1.vd1.fa
SDF=$T/GRCh38_sdf
SX="singularity exec -B $PIPE"
BCF="$SX $E/bwa_samtools_fastp_fgbio_bcftools.sif bcftools"
RTG="$SX $E/rtg_tools.sif rtg"
BEDTOOLS="$SX $E/gatk.sif bedtools"
```

#### 9.1 RTG reference for vcfeval [Generate]

```bash
$RTG format -o $R/truth/GRCh38_sdf $REF
```
#### 9.2 SEQC2 somatic, HCC1395 [Download + Generate]

```bash
S2=https://ftp-trace.ncbi.nlm.nih.gov/seqc/ftp/Somatic_Mutation_WG/release/v1.2.1
wget -c -P $T $S2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz $S2/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz \
              $S2/High-Confidence_Regions_v1.2.bed $S2/sSNV.MSDUKT.superSet.v1.2.vcf.gz $S2/sINDEL.MDKT.superSet.v1.2.vcf.gz
$SX $E/bcftool_cnvkit_vep_snpeff.sif bash -c "set -euo pipefail; cd $T
  for f in high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz \
           sSNV.MSDUKT.superSet.v1.2.vcf.gz sINDEL.MDKT.superSet.v1.2.vcf.gz; do bcftools index -f -t \$f; done
  bcftools concat -a high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz \
    -Oz -o SEQC2_Truth_v1.2.1.vcf.gz
  bcftools index -f -t SEQC2_Truth_v1.2.1.vcf.gz"
# evaluation region: padded V6 + UTR targets X SEQC2 high-confidence regions
$BEDTOOLS intersect -a $R/capture/S07604624/S07604624_targets_pad100.bed.gz -b $T/High-Confidence_Regions_v1.2.bed \
  | sort -k1,1 -k2,2n | $BEDTOOLS merge -i - > $R/truth/S07604624_targets_pad100_SEQC2_highconf.bed
```

#### 9.3 GIAB HG002 germline, v4.2.1 [Download]

```bash
GB=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
curl -s $GB/ | grep -o 'href="[^"]*benchmark[^"]*"'          # confirm the file names
wget -c -P $R/truth $GB/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz $GB/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
                         $GB/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

### 10. Benchmark settings and scoring helper

```bash
OUT_S=$DATA/somatic_demo/Result     # outdir of the SEQC2 somatic run
OUT_G=$DATA/germline_demo/Result    # outdir of the HG002 germline run

HEADER="callset\ttype\ttruth\tTP\tFP\tFN\trecall\tprecision\tF1"
```
#### 10.1 Counts and metrics per variant type from one vcfeval folder

```bash
count_metrics () {    # usage: count_metrics <vcfeval_dir> <label>
  for t in snps indels; do
    TPb=$($BCF view -H -v $t $1/tp-baseline.vcf.gz | wc -l)   # truth variants found
    FN=$($BCF  view -H -v $t $1/fn.vcf.gz          | wc -l)   # truth variants missed
    TPc=$($BCF view -H -v $t $1/tp.vcf.gz          | wc -l)   # calls that are true
    FP=$($BCF  view -H -v $t $1/fp.vcf.gz          | wc -l)   # calls that are false
    awk -v l=$2 -v t=$t -v a=$TPb -v n=$FN -v c=$TPc -v f=$FP 'BEGIN{
      r=(a+n)?a/(a+n):0; p=(c+f)?c/(c+f):0; F=(p+r)?2*p*r/(p+r):0;
      printf "%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", l, t, a+n, a, f, n, r, p, F }'
  done
}
```

### 11. Somatic benchmark: SEQC2 HCC1395 / HCC1395BL (valid1)

```bash
mkdir -p $OUT_S/benchmark && cd $OUT_S/benchmark
EVAL=$R/truth/S07604624_targets_pad100_SEQC2_highconf.bed
awk '{s+=$3-$2} END{print "somatic evaluation region: " s/1e6 " Mb"}' $EVAL
```

#### 11.1 Truth: HighConf + MedConf, FILTER = PASS, no samples, split alleles

```bash
$BCF view -G $T/SEQC2_Truth_v1.2.1.vcf.gz \
  | awk 'BEGIN{OFS="\t"} /^#/{print; next} {$7="PASS"; print}' \
  | $BCF norm -a --atom-overlaps . -m -any -f $REF | $BCF sort -Oz -o truth_somatic.vcf.gz -
$BCF index -f -t truth_somatic.vcf.gz
```

#### 11.2 Superset and its non-truth labels (FILTER column)

```bash
$BCF concat -a $T/sSNV.MSDUKT.superSet.v1.2.vcf.gz $T/sINDEL.MDKT.superSet.v1.2.vcf.gz \
  | $BCF sort | $BCF norm -a --atom-overlaps . -m -any -f $REF | $BCF sort -Oz -o superset.vcf.gz -

$BCF index -f -t superset.vcf.gz

$BCF query -f '%FILTER\n' superset.vcf.gz | tr ';' '\n' \
  | grep -E '^(HighConf|MedConf|LowConf|Unclassified|REJECT)$' | sort | uniq -c   # label counts

$BCF view -i 'FILTER~"LowConf"'      superset.vcf.gz -Oz -o lowconf.vcf.gz      && $BCF index -f -t lowconf.vcf.gz
$BCF view -i 'FILTER~"Unclassified"' superset.vcf.gz -Oz -o unclassified.vcf.gz && $BCF index -f -t unclassified.vcf.gz

$BCF isec -n=2 -w1 truth_somatic.vcf.gz lowconf.vcf.gz      | grep -vc '^#'    # truth sites labelled LowConf
$BCF isec -n=2 -w1 truth_somatic.vcf.gz unclassified.vcf.gz | grep -vc '^#'    # truth sites labelled Unclassified
$BCF isec -C -w1 lowconf.vcf.gz      truth_somatic.vcf.gz -Oz -o lowconf.x.vcf.gz      && mv lowconf.x.vcf.gz lowconf.vcf.gz           && $BCF index -f -t lowconf.vcf.gz
$BCF isec -C -w1 unclassified.vcf.gz truth_somatic.vcf.gz -Oz -o unclassified.x.vcf.gz && mv unclassified.x.vcf.gz unclassified.vcf.gz && $BCF index -f -t unclassified.vcf.gz
```

#### 11.3 Score one call set

Scoring rules: novel calls (not in the superset) with tumor VAF < 5% are dropped, and calls at LowConf sites
are dropped. The `main` score also drops calls at Unclassified sites; the `strict` score keeps them, so they count
as false positives.

```bash
score_somatic () {    # usage: score_somatic <label> <vcf> [bcftools view selection]
  local L=$1 V=$2; shift 2
  local S=$($BCF query -l $V | grep -iE '(_T|TUMOR)$' | head -1)
  [ -n "$S" ] || { echo "ERROR: no tumor column in $V"; return 1; }
  echo "== $L   tumor column: $S"

  $BCF view -s "$S" "$@" $V | $BCF norm -a --atom-overlaps . -m -any -f $REF | $BCF sort | $BCF norm -d exact -Oz -o $L.calls.vcf.gz -
  $BCF index -f -t $L.calls.vcf.gz

  # (i) novel calls: keep only tumor VAF ≥ 5%
  $BCF isec -C -w1 $L.calls.vcf.gz superset.vcf.gz -Oz -o $L.novel.vcf.gz
  $BCF view -i 'FMT/AF>=0.05' $L.novel.vcf.gz -Oz -o $L.novel_keep.vcf.gz
  $BCF index -f -t $L.novel_keep.vcf.gz

  # calls at superset sites; (ii) drop LowConf
  $BCF isec -n=2 -w1 $L.calls.vcf.gz superset.vcf.gz -Oz -o $L.known.vcf.gz && $BCF index -f -t $L.known.vcf.gz
  $BCF isec -C -w1 $L.known.vcf.gz lowconf.vcf.gz -Oz -o $L.known_strict.vcf.gz && $BCF index -f -t $L.known_strict.vcf.gz
  $BCF isec -C -w1 $L.known_strict.vcf.gz unclassified.vcf.gz -Oz -o $L.known_main.vcf.gz && $BCF index -f -t $L.known_main.vcf.gz

  for s in main strict; do
    $BCF concat -a $L.novel_keep.vcf.gz $L.known_$s.vcf.gz | $BCF view -G | $BCF sort -Oz -o $L.$s.vcf.gz -
    $BCF index -f -t $L.$s.vcf.gz
    rm -rf vcfeval_${L}_$s
    $RTG vcfeval -b truth_somatic.vcf.gz -c $L.$s.vcf.gz -e $EVAL -t $SDF \
                 --sample ALT --all-records -o vcfeval_${L}_$s > vcfeval_${L}_$s.log 2>&1
    count_metrics vcfeval_${L}_$s ${L}_$s | tee -a somatic_metrics.tsv
  done
}

echo -e "$HEADER" > somatic_metrics.tsv
```

#### 11.4 Call sets

```bash
find $OUT_S -name "valid1*.vcf.gz" | sort    # list the VCFs, then set the paths below
CONS=$(find $OUT_S -name "valid1*_consensus.vcf.gz" | head -1)
FILT=$(find $OUT_S -name "valid1*_filtered.vcf.gz"  | head -1)
# per-caller VCFs (optional; set the paths); if not in the outdir, take them from the work directory
M2=$OUT_S/05_variant_calling/somatic/snv/mutect2/valid1_T_mutect2.vcf.gz
STRELKA_SNV=$OUT_S/05_variant_calling/somatic/snv/strelka2/valid1_T_strelka2_snvs.af.vcf.gz
STRELKA_INDEL=$OUT_S/05_variant_calling/somatic/snv/strelka2/valid1_T_strelka2_indels.af.vcf.gz
MUSE=$OUT_S/05_variant_calling/somatic/snv/muse/valid1_T_muse.af.vcf.gz
VARSCAN_SNV=$OUT_S/05_variant_calling/somatic/snv/varscan/valid1_T_varscan_snv.af.vcf.gz
VARSCAN_INDEL=$OUT_S/05_variant_calling/somatic/snv/varscan/valid1_T_varscan_indel.af.vcf.gz
```

#### 11.5 Final call set (Table 2) and consensus

```bash
score_somatic final_pass     $FILT -f PASS     # consensus + VAF/alt reads/depth/normal filters
score_somatic consensus_pass $CONS -f PASS     # ≥ 2 callers
score_somatic consensus_all  $CONS             # all consensus records, LOWCONSENSUS included (plus non-Mutect2 sites with consensus_backbone: union)
```

#### 11.6 Each caller before (raw) and after (pass) its own filter (Figure 3 and Table 2)

```bash
$BCF concat -a $STRELKA_SNV $STRELKA_INDEL | $BCF sort -Oz -o strelka2.vcf.gz - && $BCF index -f -t strelka2.vcf.gz
$BCF concat -a $VARSCAN_SNV $VARSCAN_INDEL | $BCF sort -Oz -o varscan.vcf.gz -  && $BCF index -f -t varscan.vcf.gz
score_somatic mutect2_raw   $M2
score_somatic mutect2_pass  $M2             -f PASS
score_somatic strelka2_raw  strelka2.vcf.gz
score_somatic strelka2_pass strelka2.vcf.gz -f PASS
score_somatic muse_raw      $MUSE
score_somatic muse_pass     $MUSE           -e 'FILTER~"Tier5"'
score_somatic varscan_raw   varscan.vcf.gz
score_somatic varscan_pass  varscan.vcf.gz  -i 'INFO/SS=2'

column -t somatic_metrics.tsv
```

### 12. Germline benchmark: GIAB HG002

```bash
mkdir -p $OUT_G/benchmark && cd $OUT_G/benchmark

echo "$OUT_G"; ls "$OUT_G"
CALLS=$(find $OUT_G -iname "*hg002*_filtered.vcf.gz"  | head -1); echo "calls: $CALLS"
THR=$(find   $OUT_G -iname "*hg002*thresholds.bed.gz" | head -1); echo "thresholds: $THR"
```

#### 12.1 Callable targets: ≥ 90% of bases at ≥ 10X

```bash
zcat $THR | awk 'NR==1{for(i=1;i<=NF;i++) if($i=="10X") c=i; if(!c){print "ERROR: no 10X column" > "/dev/stderr"; exit 1}; next}
                 ($3-$2)>0 && $c/($3-$2)>=0.9 {print $1"\t"$2"\t"$3}' \
  | sort -k1,1 -k2,2n > callable.bed
```

#### 12.2 Evaluation region: callable targets inside GIAB confident regions

```bash
$BEDTOOLS intersect -a callable.bed -b $T/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  | sort -k1,1 -k2,2n > callable_giab.bed
$BEDTOOLS merge -i callable_giab.bed > HG002_eval.bed
awk '{s+=$3-$2} END{print "germline evaluation region: " s/1e6 " Mb"}' HG002_eval.bed
```

#### 12.3 Calls

```bash
$BCF index -f -t $CALLS 2>/dev/null; cp $CALLS hg002_calls.vcf.gz; cp $CALLS.tbi hg002_calls.vcf.gz.tbi
```
#### 12.4 Genotype-aware scoring: before filter (all records) and after filter (PASS)

```bash
echo -e "$HEADER" > germline_metrics.tsv
for m in before_filter after_filter; do
  opt=""; [ $m = before_filter ] && opt="--all-records"
  rm -rf vcfeval_germline_$m
  $RTG vcfeval -b $T/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -c hg002_calls.vcf.gz \
               -e HG002_eval.bed -t $SDF $opt -o vcfeval_germline_$m > vcfeval_germline_$m.log 2>&1
  count_metrics vcfeval_germline_$m germline_$m | tee -a germline_metrics.tsv
done
column -t germline_metrics.tsv
```

### 13. RNA fusion benchmark: HCC1395 RNA-seq vs the Talsania et al. 2022 fusion list

```bash
OUT_R=$DATA/rna_demo/Result         # outdir of the RNA fusion run
SCRIPT=$PIPE/fusion_truth_match.py          # set: location of the script
TRUTH=$T/Talsania_2022_Add_File9.csv        # set: Talsania Additional file 9, saved as CSV
mkdir -p $OUT_R/benchmark && cd $OUT_R/benchmark
```

#### 13.1 Calls and gene-name table

```bash
CALLS=$(find $OUT_R -name "*_fusions_annotated.tsv" | head -1); echo "calls: $CALLS"
# HGNC table to translate renamed genes (downloaded once into the truth folder)
[ -s $T/hgnc_complete_set.txt ] || \
  wget -O $T/hgnc_complete_set.txt https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
head -1 $TRUTH | tr ',' '\n' | grep -n "fusion"    # check: the column "RNA-seq fusion_genes" must be listed
```

#### 13.2 Match calls to the list (gene order, renamed genes, intergenic lists handled)

```bash
python3 $SCRIPT --truth $TRUTH --calls $CALLS --hgnc $T/hgnc_complete_set.txt --out fusion_truth_result.tsv
```

#### 13.3 Recovered fusions by consensus status and tier

```bash
echo "== recovered, by consensus verdict"
awk -F'\t' 'NR>1 && $4=="FOUND"{n[$7]++} END{for(k in n) print k"\t"n[k]}' fusion_truth_result.tsv | sort
echo "== recovered, by tier"
awk -F'\t' 'NR>1 && $4=="FOUND"{n[$8]++} END{for(k in n) print k"\t"n[k]}' fusion_truth_result.tsv | sort
echo "== missed, but at least one partner gene appears in the calls"
awk -F'\t' 'NR>1 && $4=="missed" && $9!="-"' fusion_truth_result.tsv | wc -l

```
#### 13.4 Calls that are not in the reference list

```bash
total=$(grep -v '^#' $CALLS | tail -n +2 | wc -l)
matched=$(awk -F'\t' 'NR>1 && $4=="FOUND"{print $5}' fusion_truth_result.tsv | sort -u | wc -l)
echo "pipeline fusions: $total   matching a listed fusion: $matched   not in the list: $((total - matched))"
```

### 14. Figures and tables

Needs Python 3 with pandas, numpy and matplotlib. Figure 3 uses the rows written by steps 11.5 and 11.6.

```bash
mkdir -p $DATA/figures && cd $DATA/figures
python3 $PIPE/make_figures.py qc  --results $OUT_S --tumor valid1_T --normal valid1_N -o Figure2
python3 $PIPE/make_figures.py cnv --results $OUT_S --tumor valid1_T -o Figure4
python3 $PIPE/benchmark_plots.py somatic  --somatic $OUT_S/benchmark/somatic_metrics.tsv -o Figure3
python3 $PIPE/benchmark_plots.py combined --germline $OUT_G/benchmark/germline_metrics.tsv \
                                          --fusion $OUT_R/benchmark/fusion_truth_result.tsv -o Figure5
python3 $PIPE/trace_table.py $OUT_S/pipeline_trace.txt -o Table1_somatic.tsv    # Table 1, one block per run
python3 $PIPE/trace_table.py $OUT_G/pipeline_trace.txt -o Table1_germline.tsv
python3 $PIPE/trace_table.py $OUT_R/pipeline_trace.txt -o Table1_rna.tsv
```

Table 2 is `$OUT_S/benchmark/somatic_metrics.tsv` (steps 11.5-11.6).

---

## Code availability and citation

Code: https://github.com/chauhan2987/PAPIO (release v1.0.0). Containers: `oras://ghcr.io/chauhan2987/papio_containers`
(tags listed in "Containers"). The test bundle and expected outputs are archived on Zenodo (DOI: [to be added]).

If you use PAPIO, please cite the protocol and the tools listed in the Containers table.

## Contact

Pankaj K. Chauhan - PKChauhan@mdanderson.org