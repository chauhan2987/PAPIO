# PAPIO - Pan-omic Analysis PIpeline for Oncology (v1.0.0)

PAPIO is a Nextflow DSL2 workflow for tumor exome and RNA sequencing. From FASTQ files it produces
germline variants, somatic SNVs and indels (four-caller consensus), copy number and LOH with purity and ploidy,
structural variants, microsatellite instability, tumor mutational burden and RNA gene fusions, with
per-sample QC and one consolidated report per patient. It runs on an HPC cluster (LSF) with Singularity/Apptainer images.

This README covers installation, reference preparation, inputs, running the pipeline and its outputs.
The demonstration data use the Agilent SureSelect V6 + UTR exome kit (SEQC2 HCC1395/HCC1395BL, GIAB HG002);
---

## Contents

1. [Analysis arms](#analysis-arms)
2. [Requirements](#requirements)
3. [Repository layout](#repository-layout)
4. [Containers](#containers)
5. [Reference resources (step by step)](#reference-resources-step-by-step)
6. [Samplesheets](#samplesheets)
7. [Parameter file](#parameter-file)
8. [Running PAPIO](#running-papio)
9. [Processes in each arm](#processes-in-each-arm)
10. [Outputs](#outputs)
11. [Demonstration data and benchmarking](#demonstration-data-and-benchmarking)
12. [Troubleshooting](#troubleshooting)
13. [Citation](#citation)

---

## Analysis arms

| Arm | Set in the parameter file | Input | Main results |
|---|---|---|---|
| Somatic | `assay: wes`, `analysis: somatic` | Matched tumor (DT) and normal (DN) exomes | Somatic SNVs/indels (VCF, MAF), CNV, LOH, purity/ploidy, SVs, MSI, TMB, patient report |
| Germline | `assay: wes`, `analysis: germline` | Normal exomes (DN) | Germline SNVs/indels, SVs, (CNVs for >= 30 normals), germline report |
| Germline within a somatic run | `analysis: somatic`, `germline_on_normals: true` | As somatic | Somatic results plus a germline report for each matched normal |
| RNA fusion | `assay: rna`, `analysis: rna_fusion` | Tumor RNA-seq (RN) | RNA QC, fusions from Arriba and STAR-Fusion, FusionInspector evidence, tiered fusion report |

Tumor-only analysis (`analysis: tumor_only`) is not implemented; the pipeline stops with an error.

---

## Requirements

| Item | Version used | Note |
|---|---|---|
| Nextflow | 26.04.1 or later | Strict-syntax parser; earlier releases fail at parse time |
| Java | 17 or later | For Nextflow |
| Singularity / Apptainer | 3.8+ / 1.2+ | All tools run from images |
| Scheduler | LSF | Queues set with `lsf_queue_*`; other schedulers need a new profile |
| Disk | ~250 GB references; ~150 GB work space per tumor-normal exome pair | Work directory on scratch |
| Memory | up to 93 GB for duplicate marking of a deep normal exome | Other steps <= 64 GB |

Typical run time for one tumor-normal exome pair: about 3.5 h wall clock and about 25 CPU-hours.

---

## Repository layout

```
genomics_pipeline/
|-- nextflow/                 # pipeline code (main.nf, nextflow.config, workflows_*.nf, components_*.nf, utils.nf)
|-- envs/                     # Singularity images (*.sif)
|-- references/               # reference data (section 5)
|-- README.md
```

Run all commands from `genomics_pipeline/`.

```bash
git clone https://github.com/chauhan2987/PAPIO.git genomics_pipeline
cd genomics_pipeline
```

---

## Containers

The pipeline looks for these files in `envs/` (`container_dir`). The names are set in `nextflow/nextflow.config`.

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
| `container_star` | `star_samtools.sif` | STAR 2.7.11b, samtools 1.23.1|
| `container_arriba` | `arriba.sif` | Arriba 2.5.1 |
| `container_starfusion` | `starfusion_fusioninspector_ctat.sif` | STAR-Fusion 1.15.1, FusionInspector 2.10.0, CTAT (GRCh38 GENCODE v44, Oct 2023 plug-n-play)|
| `container_fusioncatcher` | `fusioncatcher.sif` | FusionCatcher 1.33 |
| (benchmarking only) | `rtg_tools.sif` | RTG Tools 3.12.1 (vcfeval) |

---

## Reference resources (step by step)

Every step is marked **[Download]** or **[Generate]**. Output names match the parameter file in the next sections.
Annotation versions: **DNA** - Ensembl 115 (VEP cache, gene BED, MANE CDS); **RNA** - GENCODE v44 (STAR index,
RNA QC files) and CTAT library v44 (STAR-Fusion, FusionInspector), so all RNA tools use the same gene model.

### 0. Set-up

```bash
P=$PWD                                  # pipeline root, e.g. $HOME/genomics_pipeline
R=$P/references
E=$P/envs
REF=$R/GRCh38/GRCh38.d1.vd1.fa
DICT=$R/GRCh38/GRCh38.d1.vd1.dict
U=$R/DBs_Unused_Base                    # raw downloads land here
SX="singularity exec -B $P"
BCF="$SX $E/bwa_samtools_fastp_fgbio_bcftools.sif"
GATK="$SX $E/gatk.sif"
mkdir -p $U $R/{GRCh38,BostonGeneBed,ValidationBed,DBs_SNP_Indels,DBs_Somatic_VCFs,DBs_Germline_VCFs,DBs_Annot} \
         $R/{DBs_Analytic_Metrics,DB_VEP_Cache,DBs_AnnotSV,DBs_RNA,DBs_CTAT,DB_Arriba,DBs_Fusion_Catcher} \
         $R/{RefFlat,RNA_Intervals_Dir,DBs_gCNV/scatter_flat,truth_data} $E
```

Long steps (bwa-mem2 index, STAR index, VEP cache) need an interactive LSF session:

```bash
bsub -Is -q long -W 24:00 -M 64 -R "rusage[mem=64]" -n 16 /bin/bash
```

After opening the session, `cd` to the pipeline root and run the set-up block again (variables do not carry over).

---

### 1. Containers

All images live in `envs/`. Their file names must match the `container_*` parameters in `nextflow/nextflow.config`
(table in "Containers" above). Every image is also published at `oras://ghcr.io/chauhan2987/papio_containers:<tag>`:

```bash
for t in gatk bwa_samtools_fastp_fgbio_bcftools bcftool_cnvkit_vep_snpeff msisensorpro star_samtools; do
  singularity pull $E/$t.sif oras://ghcr.io/chauhan2987/papio_containers:$t
done
```

The images below can also be rebuilt from public sources.

#### 1.1 Single-tool images [Download]

```bash
cd $E
apptainer pull strelka2.sif      docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1
apptainer pull annotsv.sif       docker://quay.io/biocontainers/annotsv:3.5.3--py313hdfd78af_0
apptainer pull starfusion.sif    docker://trinityctat/starfusion:1.15.1
apptainer pull fusioncatcher.sif docker://clinicalgenomics/fusioncatcher:1.33
apptainer pull rtg_tools.sif     docker://biocontainers/rtg-tools:v3.12.1dfsg-1-deb_cv1
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
    vcf2maf=1.6.21 samtools=1.23.1 htslib=1.23.1 \
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
tar -xzf $U/GRCh38.d1.vd1.fa.tar.gz -C $R/GRCh38
ls $REF
```

#### 2.2 FASTA index, bwa-mem2 index, sequence dictionary [Generate]

```bash
$BCF bash -c "samtools faidx $REF && bwa-mem2 index $REF"
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

#### 3.2 Somalier sites [Download]

```bash
wget -c -O $R/DBs_Somatic_VCFs/somalier_hg38.vcf.gz https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
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
reconcile $U/Homo_sapiens_assembly38.dbsnp138.vcf.gz          $R/DBs_SNP_Indels/Homo_sapiens_assembly38.dbsnp138.d1vd1.vcf.gz
reconcile $U/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz $R/DBs_SNP_Indels/Mills_and_1000G_gold_standard.indels.hg38.d1vd1.vcf.gz
reconcile $U/Homo_sapiens_assembly38.known_indels.vcf.gz      $R/DBs_SNP_Indels/Homo_sapiens_assembly38.known_indels.d1vd1.vcf.gz
reconcile $U/af-only-gnomad.hg38.vcf.gz                       $R/DBs_Somatic_VCFs/af-only-gnomad.d1vd1.vcf.gz
reconcile $U/small_exac_common_3.hg38.vcf.gz                  $R/DBs_Somatic_VCFs/small_exac_common_3.d1vd1.vcf.gz
reconcile $U/1000g_pon.hg38.vcf.gz                            $R/DBs_Somatic_VCFs/1000g_pon.d1vd1.vcf.gz
reconcile $U/hapmap_3.3.hg38.vcf.gz                           $R/DBs_Germline_VCFs/hapmap_3.3.d1vd1.vcf.gz
reconcile $U/1000G_omni2.5.hg38.vcf.gz                        $R/DBs_Germline_VCFs/1000g_omni2.5.d1vd1.vcf.gz
reconcile $U/1000G_phase1.snps.high_confidence.hg38.vcf.gz    $R/DBs_Germline_VCFs/1000g_phase1_snps_high_conf.d1vd1.vcf.gz
'
```

Check one file: the record count drops only slightly, and no contig outside the dictionary remains.

```bash
F=$R/DBs_SNP_Indels/Mills_and_1000G_gold_standard.indels.hg38.d1vd1.vcf.gz
zcat $U/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | grep -vc '^#'
zcat $F | grep -vc '^#'
comm -23 <(zcat $F | grep -v '^#' | cut -f1 | sort -u) <(grep '^@SQ' $DICT | sed 's/.*SN:\([^\t]*\).*/\1/' | sort)   # must print nothing
$GATK gatk ValidateVariants -V $F -R $REF --validation-type-to-exclude ALL
```

---

### 4. Capture intervals

Design BEDs come from Agilent SureDesign (login required): https://earray.chem.agilent.com/suredesign/

| Kit | Design | Files | Folder |
|---|---|---|---|
| Demonstration (SEQC2, HG002): SureSelect Human All Exon V6 + UTR | S07604624 | `S07604624_Regions.bed`, `S07604624_Covered.bed` | `references/ValidationBed/` |
| Own data (RCC cohort): SureSelect Human All Exon V7 | S31285117 | `S31285117_Regions.bed`, `S31285117_MergedProbes.bed` | `references/BostonGeneBed/` |

Check that a BED is on hg38 (EGFR exons lie between chr7:55.0 and 55.2 Mb):

```bash
cut -f1 $R/ValidationBed/S07604624_Regions.bed | sort -u | head
grep -w chr7 $R/ValidationBed/S07604624_Regions.bed | awk '$2>55019000 && $3<55211000' | head -3
```

Sections 4.1-4.4 are run once per kit. The demonstration (SEQC2, HG002) uses the V6 + UTR kit; for your own
data, comment the V6 block and uncomment the V7 block, then run 4.1-4.4 again.

```bash
# --- demonstration kit: SureSelect V6 + UTR (S07604624) ---
K=$R/ValidationBed;  ID=S07604624; PROBES=$K/S07604624_Covered.bed
BAITS=$K/Covered_Baits.interval_list; TARGETS=$K/Targets_Regions.interval_list
SOM_SHARDS=$R/Validation_Scatter_Intervals_Dir_V6_pad100; GERM_SHARDS=$R/Validation_Germline_Scatter_150_V6_pad100

# --- own data: SureSelect V7 (S31285117) ---
# K=$R/BostonGeneBed;  ID=S31285117; PROBES=$K/S31285117_MergedProbes.bed
# BAITS=$K/MergedProbes.interval_list;  TARGETS=$K/Regions.interval_list
# SOM_SHARDS=$R/Scatter_Intervals_Dir_V7_pad100;            GERM_SHARDS=$R/Germline_Scatter_150_V7_pad100
```

#### 4.1 Bait and target interval lists [Generate]

```bash
$GATK gatk BedToIntervalList -I $PROBES             -SD $DICT -O $BAITS
$GATK gatk BedToIntervalList -I $K/${ID}_Regions.bed -SD $DICT -O $TARGETS
```

#### 4.2 Clean, sorted target BED [Generate]

Vendor BEDs start with `track`/`browser` lines; this copy is used for evaluation regions.

```bash
awk 'BEGIN{OFS="\t"} !/^#/ && !/^track/ && !/^browser/ && NF>=3 {print $1,$2,$3}' $K/${ID}_Regions.bed \
  | LC_ALL=C sort -k1,1 -k2,2n > $K/${ID}_Regions.sorted.bed
head -3 $K/${ID}_Regions.sorted.bed
```

#### 4.3 Padded calling shards: somatic (24) and germline (150) [Generate]

Targets +- 100 bp; overlapping padded targets are merged; no target is split between shards.

```bash
$GATK gatk SplitIntervals -R $REF -L $K/${ID}_Regions.bed --interval-padding 100 --interval-merging-rule ALL \
  --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION --scatter-count 24  -O $SOM_SHARDS
$GATK gatk SplitIntervals -R $REF -L $K/${ID}_Regions.bed --interval-padding 100 --interval-merging-rule ALL \
  --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION --scatter-count 150 -O $GERM_SHARDS
```

#### 4.4 Padded call regions for Strelka2 and Manta (`call_regions_gz`) [Generate]

Built from the somatic shards, so all callers see the same territory.

```bash
grep -hv '^@' $SOM_SHARDS/*.interval_list \
  | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' | LC_ALL=C sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"} {if($1==c && $2<=e){if($3>e)e=$3} else {if(c!="")print c,s,e; c=$1;s=$2;e=$3}} END{if(c!="")print c,s,e}' \
  > $K/${ID}_Regions.pad100.bed
$BCF bash -c "bgzip -f $K/${ID}_Regions.pad100.bed && tabix -f -p bed $K/${ID}_Regions.pad100.bed.gz"
```

Params YAML entries for the two kits:

| Parameter | Demonstration (V6 + UTR) | Own data (V7) |
|---|---|---|
| `interval_bed` | `ValidationBed/S07604624_Regions.bed` | `BostonGeneBed/S31285117_Regions.bed` |
| `call_regions_gz` | `ValidationBed/S07604624_Regions.pad100.bed.gz` | `BostonGeneBed/S31285117_Regions.pad100.bed.gz` |
| `bait_interval_list` | `ValidationBed/Covered_Baits.interval_list` | `BostonGeneBed/MergedProbes.interval_list` |
| `target_interval_list` | `ValidationBed/Targets_Regions.interval_list` | `BostonGeneBed/Regions.interval_list` |
| `scatter_intervals` | `Validation_Scatter_Intervals_Dir_V6_pad100` | `Scatter_Intervals_Dir_V7_pad100` |
| `germline_scatter_intervals` | `Validation_Germline_Scatter_150_V6_pad100` | `Germline_Scatter_150_pad100` |

---

### 5. Copy-number and MSI resources

#### 5.1 Sequenza and MSI list [Generate]

```bash
$SX $E/sequenza.sif sequenza-utils gc_wiggle --fasta $REF -w 50 -o $R/DBs_Somatic_VCFs/sequenza_d1.vd1_gc50.wig.gz
$SX $E/msisensorpro.sif msisensor-pro scan -d $REF -o $R/DBs_Somatic_VCFs/msisensor_scan.list
```

#### 5.2 Germline copy-number resources (only with `run_gcnv: true` and ≥ 30 normals) [Generate]

```bash
G=$R/DBs_gCNV
GCNV_BED=$R/ValidationBed/S07604624_Regions.bed          # demonstration kit
# GCNV_BED=$R/BostonGeneBed/S31285117_Regions.bed        # own data (V7)
$GATK gatk PreprocessIntervals -R $REF -L $GCNV_BED \
  --bin-length 0 --padding 250 --interval-merging-rule OVERLAPPING_ONLY -O $G/preprocessed_intervals.interval_list
$GATK gatk AnnotateIntervals -R $REF -L $G/preprocessed_intervals.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY -O $G/annotated_intervals.tsv
$GATK gatk IntervalListTools --INPUT $G/preprocessed_intervals.interval_list \
  --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_COUNT 24 --OUTPUT $G/scatter
i=0; for d in $G/scatter/temp_*; do i=$((i+1)); printf -v n '%04d' $i; cp "$d/scattered.interval_list" "$G/scatter_flat/shard_${n}.interval_list"; done
awk 'BEGIN{OFS="\t"; print "CONTIG_NAME","PLOIDY_PRIOR_0","PLOIDY_PRIOR_1","PLOIDY_PRIOR_2","PLOIDY_PRIOR_3"
  for (i=1;i<=22;i++) print "chr"i,"0.01","0.01","0.97","0.01"
  print "chrX","0.01","0.49","0.49","0.01"; print "chrY","0.495","0.495","0.01","0.01"}' > $G/contig_ploidy_priors.tsv
```

CNVkit runs against the matched normal (`cnvkit_pon: null`); a panel of normals can be added later.

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
  | LC_ALL=C sort -k1,1 -k2,2n > $R/GRCh38/GRCh38_ensembl_115_genes.bed

# merge overlapping intervals of a sorted BED
merge3() { awk 'BEGIN{OFS="\t"} {if($1==c && $2<=e){if($3>e)e=$3} else {if(c!="")print c,s,e; c=$1;s=$2;e=$3}} END{if(c!="")print c,s,e}'; }

# tmb_coding_bed: MANE Select coding exons, merged
zcat $GTF | awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS" && $9 ~ /tag "MANE_Select"/ {
    c=$1; if (c=="MT") c="M"; if (c !~ /^([0-9]+|X|Y|M)$/) next
    print "chr" c, $4-1, $5 }' | LC_ALL=C sort -k1,1 -k2,2n | merge3 > $R/GRCh38/GRCh38_ens115_MANE_CDS.bed

# germline_report_genes: your own list, one symbol per line
cat > $R/DBs_Analytic_Metrics/reportable_genes.txt <<'EOF'
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
zcat $GTF | awk -F'\t' -v list=$R/DBs_Analytic_Metrics/reportable_genes.txt 'BEGIN{OFS="\t"; while ((getline g < list) > 0) G[toupper(g)]=1}
  $3=="CDS" && $9 ~ /tag "MANE_Select"/ && match($9, /gene_name "[^"]+"/) {
    g=substr($9,RSTART+11,RLENGTH-12); c=$1; if (c=="MT") c="M"
    if (c ~ /^([0-9]+|X|Y|M)$/ && (toupper(g) in G)) print "chr" c, $4-1, $5, g }' \
  | LC_ALL=C sort -k1,1 -k2,2n -u > $R/DBs_Analytic_Metrics/reportable_genes_MANE_CDS.bed
```

#### 6.3 Chromosome-arm BED [Download + Generate]

```bash
wget -c -O $R/DBs_Analytic_Metrics/cytoBand.hg38.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
zcat $R/DBs_Analytic_Metrics/cytoBand.hg38.txt.gz | awk -F'\t' 'BEGIN{OFS="\t"}
  $1 ~ /^chr([0-9]+|X|Y)$/ { k=$1 SUBSEP substr($4,1,1)
    if (!(k in s) || $2<s[k]) s[k]=$2
    if (!(k in e) || $3>e[k]) e[k]=$3 }
  END { for (k in s) { split(k,a,SUBSEP); arm=substr(a[1],4) a[2]
          if (arm ~ /^(13|14|15|21|22)p$/) continue      # acrocentric p arms: no assembled sequence
          print a[1], s[k], e[k], arm } }' \
  | LC_ALL=C sort -k1,1 -k2,2n > $R/DBs_Analytic_Metrics/GRCh38_arms.bed
```

---

### 7. Variant annotation databases

#### 7.1 VEP cache [Download]

The cache release must match the VEP version in the container (115).

```bash
$SX $E/bcftool_cnvkit_vep_snpeff.sif vep --help 2>&1 | grep -i "ensembl-vep" | head -2
curl -L -o $U/homo_sapiens_vep_115_GRCh38.tar.gz \
  https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar -xzf $U/homo_sapiens_vep_115_GRCh38.tar.gz -C $R/DB_VEP_Cache
```

#### 7.2 Chromosome-name map for ClinVar and COSMIC [Generate]

```bash
for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' $c $c; done > $R/DBs_Annot/chr_map.txt
printf 'MT\tchrM\n' >> $R/DBs_Annot/chr_map.txt
```

#### 7.3 ClinVar [Download + Generate]

```bash
wget -c -P $R/DBs_Annot https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
$BCF bash -c "set -euo pipefail; cd $R/DBs_Annot
  tabix -f -p vcf clinvar.vcf.gz
  bcftools annotate --rename-chrs chr_map.txt clinvar.vcf.gz -Oz -o clinvar_chr.vcf.gz
  tabix -f -p vcf clinvar_chr.vcf.gz"
```

#### 7.4 COSMIC v104 (licence required) [Download + Generate]

COSMIC's scripted download takes two calls: the first, with your COSMIC login, returns a temporary link;
the second downloads the file. File paths are listed at https://cancer.sanger.ac.uk/cosmic/download (GRCh38, v104).

```bash
cd $R/DBs_Annot
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
tar -xf Cosmic_CancerGeneCensus_Tsv_v104_GRCh38.tar           # -> Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz
tar -xf Cosmic_GenomeScreensMutant_Normal_Vcf_v104_GRCh38.tar  # -> Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.vcf.gz
# chr names, sorted, transcript-level duplicate rows collapsed to one record per variant
$BCF bash -c "set -euo pipefail
  bcftools annotate --rename-chrs chr_map.txt Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.vcf.gz -Ou \
    | bcftools sort -Ou | bcftools norm -d exact -Oz -o cosmic_genom_mut_norm_v104_chr_dedup.vcf.gz
  tabix -f -p vcf cosmic_genom_mut_norm_v104_chr_dedup.vcf.gz"
cd $P
```

#### 7.5 AnnotSV annotations [Download]

```bash
wget -c -P $U https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz
tar -xzf $U/Annotations_Human_3.5.tar.gz -C $R/DBs_AnnotSV
ls $R/DBs_AnnotSV
```

---

### 8. RNA fusion resources

#### 8.1 GENCODE v44 annotation [Download]

v44 is the release the CTAT library (8.2) is built on, so STAR, Arriba, RNA QC and STAR-Fusion share one gene model.

```bash
wget -c -P $R/DBs_RNA https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -kf $R/DBs_RNA/gencode.v44.annotation.gtf.gz      # the params YAML uses the uncompressed .gtf
RGTF=$R/DBs_RNA/gencode.v44.annotation.gtf
```

#### 8.2 CTAT genome library, GENCODE v44 [Download]

Check the file name in https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
```bash
wget -c -P $U https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
tar -xzf $U/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz -C $R/DBs_CTAT
ls $R/DBs_CTAT/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir | head
```

#### 8.3 Arriba database files [Download + Generate]

```bash
wget -c -P $U https://github.com/suhrig/arriba/releases/download/v2.5.1/arriba_v2.5.1.tar.gz
tar -xzf $U/arriba_v2.5.1.tar.gz -C $U
mkdir -p $R/DB_Arriba/arriba_v2.5.1
cp $U/arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz \
   $U/arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
   $U/arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3 $R/DB_Arriba/arriba_v2.5.1/
# known_fusions_tsv: uncompressed copy for the fusion annotation step
zcat $R/DB_Arriba/arriba_v2.5.1/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz > $R/DBs_RNA/known_fusions_hg38_GRCh38_v2.5.1.tsv
```

#### 8.4 FusionCatcher database (only with `use_fusioncatcher: true`) [Download]

```bash
cd $R/DBs_Fusion_Catcher
for p in aa ab ac ad; do wget -c http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.$p; done
cat human_v102.tar.gz.* > human_v102.tar.gz && rm -f human_v102.tar.gz.*
tar -xzf human_v102.tar.gz && ls human_v102 | head
cd $P
```

#### 8.5 STAR index on the DNA reference genome [Generate]

Built on the same FASTA as the DNA arm, with GENCODE v44. Use a new folder so an older index is not overwritten.
`--sjdbOverhang` is read length - 1: 149 for 2 x 150 bp reads, 99 for 2 x 100 bp; 100 works well for both.

```bash
$SX $E/star_samtools.sif STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir $R/StarIndex_gencode_v44 --genomeFastaFiles $REF --sjdbGTFfile $RGTF --sjdbOverhang 100 \
  --limitGenomeGenerateRAM 60000000000
```

#### 8.6 refFlat and rRNA interval list for Picard CollectRnaSeqMetrics [Generate]

```bash
cat > $P/make_refflat.py <<'EOF'
# GENCODE GTF -> Picard refFlat; exons written in ascending order (required by Picard, also for minus-strand genes)
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
$SX $E/bcftool_cnvkit_vep_snpeff.sif python3 $P/make_refflat.py $RGTF $R/RefFlat/gencode.v44.refFlat

# rRNA interval list = sequence dictionary header + rRNA / Mt_rRNA transcripts
{ grep '^@' $DICT
  awk -F'\t' '$3=="transcript" && $9 ~ /gene_type "(Mt_)?rRNA(_pseudogene)?"/ {
      match($9, /transcript_id "[^"]+"/); print $1"\t"$4"\t"$5"\t"$7"\t"substr($9,RSTART+15,RLENGTH-16) }' $RGTF
} > $R/RNA_Intervals_Dir/gencode.v44.rrna.interval_list
```


---

### 9. Library check before the first run: inline UMIs

Look at the base composition of the first 15 cycles of R1 and R2, for tumor and normal, in several patients.
A random block followed by a fixed block means an inline UMI. In the RCC cohort: 3 random bases, 2 fixed, then genome,
which is `umi_read_structure: 3M2S+T` (3 bases UMI, skip 2, rest template).

```bash
FQ=$HOME/test_data/fastq_hg002_new/HG002_novogene_wes_18g_rep1_R2.fq.gz      # any R1 or R2 file
zcat $FQ | head -40000 | awk 'NR%4==2 {n++; for (i=1;i<=15;i++) c[i, substr($0,i,1)]++}
  END {for (i=1;i<=15;i++) printf "pos%2d  A:%5.1f  C:%5.1f  G:%5.1f  T:%5.1f\n", i,
       100*c[i,"A"]/n, 100*c[i,"C"]/n, 100*c[i,"G"]/n, 100*c[i,"T"]/n}'
```

---


## Samplesheets

One CSV for DNA (`dna_samplesheet`) and one for RNA (`rna_samplesheet`), with the same columns.
One row per FASTQ pair (lane or chunk). Rows with the same `patient_id`, `sample_type` and `batch` are merged
into one sample.

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier; links tumor, normal and RNA of one patient |
| `sex` | `0` female, `1` male, `NA` unknown (NA switches off the sex check) |
| `tumor_type` | Free-text label, e.g. `BRCA`, `RCC`; `NA` for normal-only samples |
| `batch` | Library or batch identifier (`LB` in the read group) |
| `sample_type` | `DT` tumor DNA, `DN` normal DNA, `RN` tumor RNA |
| `run_num`, `flowcell`, `lane_num` | Read-group fields (`PU`) |
| `chunk` | `001`, `002`, … for lanes or replicate libraries of one sample |
| `read1`, `read2` | FASTQ paths (absolute, or relative to `fastq_root`) |

Sample names are derived as `<patient_id>_T`, `<patient_id>_N` and `<patient_id>_R`.

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

Use the same `patient_id` for a patient's DNA and RNA so the results can be joined.

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
pipeline_root:          /home/<user>/genomics_pipeline
container_dir:          /home/<user>/genomics_pipeline/envs
fastq_root:             /scratch/<user>/fastq
work_root:              /scratch/<user>/work_demo
outdir:                 /home/<user>/Result_demo
dna_samplesheet:        /home/<user>/samples_dna_demo.csv
rna_samplesheet:        /home/<user>/samples_rna_demo.csv

# directories the containers may read/write
bind_paths:             [/home/<user>, /scratch/<user>]

# ---------------- Reference genome ----------------
# (GDC GRCh38.d1.vd1)
ref_index_dir:          ${projectDir}/../references/GRCh38
ref_fasta:              ${projectDir}/../references/GRCh38/GRCh38.d1.vd1.fa
ref_fai:                ${projectDir}/../references/GRCh38/GRCh38.d1.vd1.fa.fai
ref_dict:               ${projectDir}/../references/GRCh38/GRCh38.d1.vd1.dict

# DNA arms - shared (somatic and germline)
# ---------------- Library ----------------
dedup_mode:             markdup         # umi_markdup| umi_consensus | markdup (demonstration: no UMIs)
umi_read_structure:     null
# dedup_mode:           umi_markdup     # own data: inline UMIs
# umi_read_structure:   '3M2S+T' # per-mate: 3bp UMI, 2bp linker skip, rest template
umi_read_structure_r2:  null     # R2 structure if it differs from R1 (e.g. '+T' for UMI on R1 only); null = same as R1
skip_trimming:          false
bqsr:                   false
make_cram:              false     # BAM->CRAM conversion (storage only). false = keep BAM.

# ---------------- Capture: demonstration = SureSelect V6 + UTR (S07604624) ----------------
interval_bed:           ${projectDir}/../references/ValidationBed/S07604624_Regions.bed
call_regions_gz:        ${projectDir}/../references/ValidationBed/S07604624_Regions.pad100.bed.gz
bait_interval_list:     ${projectDir}/../references/ValidationBed/Covered_Baits.interval_list
target_interval_list:   ${projectDir}/../references/ValidationBed/Targets_Regions.interval_list


# ---------------- Capture: own data = SureSelect V7 (S31285117) ----------------
# (reconciled to d1.vd1; need .tbi)
# interval_bed:           ${projectDir}/../references/BostonGeneBed/S31285117_Regions.bed
# call_regions_gz:        ${projectDir}/../references/BostonGeneBed/S31285117_Regions.pad100.bed.gz
# bait_interval_list:     ${projectDir}/../references/BostonGeneBed/MergedProbes.interval_list
# target_interval_list:   ${projectDir}/../references/BostonGeneBed/Regions.interval_list

# ---------------- Known sites (BQSR, germline) ----------------
# (reconciled to d1.vd1; need .tbi)
dbsnp:                  ${projectDir}/../references/DBs_SNP_Indels/Homo_sapiens_assembly38.dbsnp138.d1vd1.vcf.gz
mills:                  ${projectDir}/../references/DBs_SNP_Indels/Mills_and_1000G_gold_standard.indels.hg38.d1vd1.vcf.gz
known_indels:           ${projectDir}/../references/DBs_SNP_Indels/Homo_sapiens_assembly38.known_indels.d1vd1.vcf.gz

# ---------------- Sample QC ----------------
somalier_sites:         ${projectDir}/../references/DBs_Somatic_VCFs/somalier_hg38.vcf.gz
verifybamid_db:         null
coverage_report_bed:    ${projectDir}/../references/DBs_Analytic_Metrics/reportable_genes_MANE_CDS.bed


# ---------------- Annotation (both DNA arms) ----------------
vep_cache_dir:          ${projectDir}/../references/DB_VEP_Cache
vep_plugins_dir:        null
clinvar_vcf:            ${projectDir}/../references/DBs_Annot/clinvar_chr.vcf.gz
annotsv_dir:            ${projectDir}/../references/DBs_AnnotSV
cosmic_cgc_tsv:         ${projectDir}/../references/DBs_Annot/Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz


# ------------ Somatic arm (analysis: somatic) ----------------
scatter_intervals:      ${projectDir}/../references/Validation_Scatter_Intervals_Dir_V6_pad100
germline_resource:      ${projectDir}/../references/DBs_Somatic_VCFs/af-only-gnomad.d1vd1.vcf.gz
# scatter_intervals:      ${projectDir}/../references/Scatter_Intervals_Dir_V7_pad100
pileup_variants:        ${projectDir}/../references/DBs_Somatic_VCFs/small_exac_common_3.d1vd1.vcf.gz
# Choose either 1000G PON or your custom PON for RCC
panel_of_normals:       ${projectDir}/../references/DBs_Somatic_VCFs/1000g_pon.d1vd1.vcf.gz
# panel_of_normals:     ${projectDir}/../references/DBs_Germline_VCFs/Kidney_Can_PON.vcf.gz   # own data: in-house PON
orientation_bias:       true
cnvkit_pon:             null
msisensor_list:         ${projectDir}/../references/DBs_Somatic_VCFs/msisensor_scan.list
sequenza_gc_wiggle:     ${projectDir}/../references/DBs_Somatic_VCFs/sequenza_d1.vd1_gc50.wig.gz
gene_bed:               ${projectDir}/../references/GRCh38/GRCh38_ensembl_115_genes.bed
arm_bed:                ${projectDir}/../references/DBs_Analytic_Metrics/GRCh38_arms.bed
cosmic_vcf:             ${projectDir}/../references/DBs_Annot/cosmic_genom_mut_norm_v104_chr_dedup.vcf.gz
tmb_coding_bed:         ${projectDir}/../references/GRCh38/GRCh38_ens115_MANE_CDS.bed
germline_on_normals:    false     # true = germline report for each matched normal 

# ------------ Germline arm (analysis: germline) ----------------
germline_scatter_intervals: ${projectDir}/../references/Validation_Germline_Scatter_150_V6_pad100
# germline_scatter_intervals: ${projectDir}/../references/Germline_Scatter_150_pad100
germline_caller:        haplotypecaller   # haplotypecaller | deepvariant
germline_filter:        nv                # nv (any sample number NVScoreVariants) | vqsr (>= 30 samples)
vqsr_min_samples:       30  # minimum samples before VQSR is allowed (else use 'nv')
run_gcnv:               false             # true only with >= 30 normals of one kit
gcnv_min_samples:       30
germline_ba1_af:         0.05
hapmap:                 ${projectDir}/../references/DBs_Germline_VCFs/hapmap_3.3.d1vd1.vcf.gz
omni:                   ${projectDir}/../references/DBs_Germline_VCFs/1000g_omni2.5.d1vd1.vcf.gz
thousandg_snps:         ${projectDir}/../references/DBs_Germline_VCFs/1000g_phase1_snps_high_conf.d1vd1.vcf.gz
germline_report_genes:  ${projectDir}/../references/DBs_Analytic_Metrics/reportable_genes.txt
preprocessed_intervals:   ${projectDir}/../references/DBs_gCNV/preprocessed_intervals.interval_list
contig_ploidy_priors:     ${projectDir}/../references/DBs_gCNV/contig_ploidy_priors.tsv
gcnv_scatter_intervals:   ${projectDir}/../references/DBs_gCNV/scatter_flat
gcnv_annotated_intervals: ${projectDir}/../references/DBs_gCNV/annotated_intervals.tsv

# ---------------- RNA fusion arm(analysis: rna_fusion) ----------------
rna_use_ctat_reference: false
gtf:                    ${projectDir}/../references/DBs_RNA/gencode.v44.annotation.gtf
star_index:             ${projectDir}/../references/StarIndex_gencode_v44
refflat:                ${projectDir}/../references/RefFlat/gencode.v44.refFlat
rrna_interval_list:     ${projectDir}/../references/RNA_Intervals_Dir/gencode.v44.rrna.interval_list
rna_strandedness:       auto
ctat_lib:               ${projectDir}/../references/DBs_CTAT/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir
arriba_blacklist:       ${projectDir}/../references/DB_Arriba/arriba_v2.5.1/blacklist_hg38_GRCh38_v2.5.1.tsv.gz
arriba_known_fusions:   ${projectDir}/../references/DB_Arriba/arriba_v2.5.1/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz
arriba_protein_domains: ${projectDir}/../references/DB_Arriba/arriba_v2.5.1/protein_domains_hg38_GRCh38_v2.5.1.gff3
known_fusions_tsv:      ${projectDir}/../references/DBs_RNA/known_fusions_hg38_GRCh38_v2.5.1.tsv
use_fusioncatcher:      false
fusioncatcher_db:       ${projectDir}/../references/DBs_Fusion_Catcher/human_v102

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

## Running PAPIO

Launch from a login node or a long-running interactive job; Nextflow submits every task to LSF.

```bash
nextflow run ${HOME}/genomics_pipeline/nextflow/main.nf -c ${HOME}/genomics_pipeline/nextflow/nextflow.config \
  -params-file params_demo.yaml -profile production -resume
```

- `-resume` reuses finished tasks from `work_root`; always add it when rerunning.
- **Partial runs:** set `stop_at` (e.g. `variant_calling`). To restart later from published files, set `start_at`
  and the matching input folder (`input_bam_dir` from `markdup` on, `input_vcf_dir` for `annotation`,
  `input_annot_dir` for `classification`). The consolidated reports need upstream results of the same run,
  so they are not made when starting at `classification`; use `-resume` instead.
- **Process-specific resources:** put overrides in a small extra config and add it with a second `-c`:

```groovy
// more_time.config
process {
    withName: 'PREPROCESS_AND_ALIGN_DNA' { time = { 12.h * task.attempt }; queue = params.lsf_queue_medium }
    withName: 'METRICS'                  { time = { 5.h  * task.attempt }; queue = params.lsf_queue_medium }
}
```

A run writes `pipeline_trace.txt`, `report.html` and `timeline.html` to `outdir`.

---

## Processes in each arm

Kind: **required** - always runs; **toggle** - switched by the parameter shown; **resource** - runs only when that
reference file is set. Stages run in the order listed.

### DNA stages shared by the somatic and germline arms

| Stage | Processes | Kind |
|---|---|---|
| 0 | `PREFLIGHT_REFERENCES`; `MANIFEST_ROW` → `WRITE_MANIFEST` | toggle (`preflight_references`, `write_manifest`) |
| 1 mapping | `PREPROCESS_AND_ALIGN_DNA` (fastp → bwa-mem2 → fixmate) or `ALIGN_DNA_UMI` (fgbio); `MULTIQC_TRIMMING` | required (one, by `dedup_mode`) |
| 2 markdup | `MERGE_AND_MARK_DUPLICATES` / `MERGE_AND_MARKDUP_UMI` / `MERGE_AND_CONSENSUS`; `MULTIQC_MARKDUP` | required (one, by `dedup_mode`) |
| | `CRAM_CONVERT` | toggle (`make_cram`) |
| 3 metrics | `MOSDEPTH` → `INFER_SEX`; `SAMTOOLS_STATS`; `METRICS`; `MULTIQC_METRICS` | required |
| | `VERIFYBAMID2` | resource (`verifybamid_db`) |
| | `QC_GATE` → `QC_GATE_SUMMARY` | toggle (`qc_gate`) |
| | `SOMALIER_EXTRACT` → `SOMALIER_RELATE` → `CHECK_SAMPLE_SWAPS`; `MULTIQC_SOMALIER` | resource (`somalier_sites`) |
| 4 bqsr | `BASE_RECALIBRATOR` → `APPLY_BQSR` | toggle (`bqsr`) |

### Somatic arm

| Stage | Processes | Kind |
|---|---|---|
| 5 variant calling | `MANTA` | required |
| | `MUTECT2` (per shard) → `LEARN_READ_ORIENTATION_MODEL` | required; model = toggle (`orientation_bias`) |
| | `GETPILEUPSUMMARIES` → `CALCULATE_CONTAMINATION` | resource (`pileup_variants`) |
| | `MERGE_MUTECT2` (merge + FilterMutectCalls) | required |
| | `STRELKA2` → `TAG_STRELKA_AF`; `MUSE_CALL` → `MUSE_SUMP` → `TAG_MUSE_AF`; `VARSCAN` → `TAG_VARSCAN_AF` | toggle (`use_*_in_consensus`, all on) |
| | `ENSEMBLE_CONSENSUS` → `FILTER_SOMATIC_VCF` | toggle (`somatic_consensus`, `somatic_variant_filter`) |
| | `MUTATIONAL_SIGNATURES` | toggle (`mutational_signatures`, off) |
| | `CNVKIT` | required |
| | `SEQUENZA` → `SEQUENZA_GENE_LOH` | resource (`sequenza_gc_wiggle`; gene/arm tables need `gene_bed` / `arm_bed`) |
| | `MSISENSOR_PRO` → `CLASSIFY_MSI`; `MSISENSOR_GENE_MAP` | required; toggle (`classify_msi`); resource (`gene_bed`) |
| 6 annotation | `ANNOTATE_VCF` (VEP + ClinVar + COSMIC); `COMPUTE_TMB`; `MULTIQC_ANNOT` | required; TMB = toggle (`compute_tmb`) |
| 7 classification | `CLASSIFY_VARIANTS`; `VCF2MAF` | toggle (`classify_variants`, `generate_maf`) |
| | `ANNOTATE_SV_SOMATIC` → `CONSOLIDATE_PATIENT_REPORT` | toggle (`consolidate_somatic_reports`); SV annotation needs `annotsv_dir` |

### Germline arm

| Stage | Processes | Kind |
|---|---|---|
| 5 variant calling | `HAPLOTYPECALLER` (per shard) → `MERGE_GVCFS` → `GENOTYPE_GVCFS` → `NV_SCORE_VARIANTS` → `FILTER_VARIANT_TRANCHES` | required (default) |
| | Cohort path: `REBLOCK_GVCF` → `GENOMICSDB_IMPORT` → `GENOTYPE_GVCFS_JOINT` or `GNARLY_GENOTYPER` → `MERGE_COHORT_VCF` → VQSR (SNP, indel) → `SPLIT_COHORT_VCF` | toggle (`germline_filter: vqsr`, ≥ 30 samples) |
| | `DEEPVARIANT` instead of the above | toggle (`germline_caller: deepvariant`) |
| | `MANTA_GERMLINE` | required |
| | `GATK_COLLECT_READ_COUNTS` → `GATK_PLOIDY_MODEL` → `GATK_GCNV_CALLER` → `GATK_GCNV_POSTPROCESS` | toggle (`run_gcnv`, ≥ 30 normals) |
| | `BCFTOOLS_ROH` | toggle (`run_roh`; blocked for WES, needs `roh_af_file`) |
| 6 annotation | `ANNOTATE_VCF` (VEP + ClinVar; no COSMIC); `MULTIQC_ANNOT` | required |
| 7 classification | `CLASSIFY_VARIANTS` (germline labels); `VCF2MAF` | toggle |
| | `ANNOTATE_SV_MANTA`, `ANNOTATE_SV_GCNV` → `CONSOLIDATE_GERMLINE_REPORT` | toggle (`consolidate_germline_reports`) + `annotsv_dir` |

### RNA fusion arm

| Stage | Processes | Kind |
|---|---|---|
| 1 mapping | `STAR_ALIGN` (two-pass, chimeric output); `CRAM_CONVERT` | required; CRAM = toggle (`rna_cram_convert`) |
| 3 metrics | `MARKDUP_RNA` (QC only), `RNASEQ_METRICS`, `SAMTOOLS_STATS`, strandedness inference; `MULTIQC` | required |
| 5 variant calling | `ARRIBA`; `STARFUSION` (kickstart from STAR chimeric junctions); `FUSIONCATCHER` | toggle (`use_arriba` on, `use_starfusion` on, `use_fusioncatcher` off) |
| | `FUSION_CONSENSUS` → `FILTER_FUSION_CANDIDATES` → `FUSIONINSPECTOR` | consensus required; inspector = toggle (`use_fusioninspector`) |
| 6 annotation | `ANNOTATE_FUSIONS` (known fusions + COSMIC CGC → tiers) | required |
| 7 classification | `CONSOLIDATE_FUSION_REPORT` | required |

### Which process assigns which label

| Labels | Process | Arm |
|---|---|---|
| ONCOGENIC_STRONG, ONCOGENIC_MODERATE, VUS_CANCER_GENE, VUS, LIKELY_BENIGN | `CLASSIFY_VARIANTS` (somatic mode) | Somatic SNVs/indels |
| PATHOGENIC_OR_LIKELY, BENIGN_OR_LIKELY, LOF_CANCER_GENE_REVIEW, VUS | `CLASSIFY_VARIANTS` (germline mode) | Germline SNVs/indels |
| ACMG classes 1-5 | AnnotSV (`ANNOTATE_SV_*`) | Germline SV/CNV, somatic SVs |
| KNOWN_ONCOGENIC, KNOWN_LOW_SUPPORT, LIKELY_ONCOGENIC, VUS, LIKELY_ARTIFACT, LOW_CONFIDENCE | `ANNOTATE_FUSIONS` | RNA fusions |

All labels are computational pre-classifications for review, not AMP/ASCO/CAP or ACMG/AMP clinical tiers.

---

## Outputs

Results are grouped by stage under `outdir`, and from variant calling on by analysis type
(`somatic/`, `germline/`, `rna_fusion/`). Main files:

| Folder | Files |
|---|---|
| `01_trimming/` | fastp JSON/HTML per lane, MultiQC |
| `02_markdup/` | `<sample>_dup_metrics.json`, MultiQC |
| metrics folder | Picard (`*.hs_metrics`, `*.insert_size_metrics`, alignment), mosdepth (`*.mosdepth.summary.txt`, `*.thresholds.bed.gz`), somalier (`somalier.pairs.tsv`, `somalier.samples.tsv`), QC gate, MultiQC |
| `05_variant_calling/somatic/` | `snv/` per caller and `callers_consensus/<tumor>_filtered.vcf.gz`; `cnv/cnvkit/`, `cnv/sequenza/` (`_segments.txt`, `_alternative_solutions.txt`, `_gene_loh.tsv`, `_arm_cn.tsv`); `sv/`; `msi/` |
| `05_variant_calling/germline/` | `snv/<normal>_filtered.vcf.gz`, Manta SVs, gCNV |
| `06_annotation/` | VEP VCF, TMB; `rna_fusion/<sample>_fusions_annotated.tsv` |
| `07_classification/` | variant classification tables, MAF, consolidated patient / germline / fusion reports |
| `outdir/` | `pipeline_trace.txt`, `report.html`, `timeline.html`, run manifest |

Check the exact top-level folder names with `ls <outdir>` after the first run.

---

## Demonstration data and benchmarking

| Arm | Data | Truth set | Comparison |
|---|---|---|---|
| Somatic | SEQC2 HCC1395 (SRR7890883) / HCC1395BL (SRR7890874), SureSelect V6 + UTR; SRA SRP162370 | SEQC2 v1.2.1 high-confidence SNVs/indels and regions (step 10.2) | `rtg vcfeval --sample ALT` on site-only VCFs, SNVs and indels separately; LowConf/Unclassified/REJECT sites and novel calls < 5% VAF excluded |
| Germline | GIAB HG002, Novogene SureSelect V6 exome (github.com/novogene-europe/ngs_demo_datasets) | GIAB v4.2.1 (step 10.2) | genotype-aware `rtg vcfeval` in callable targets ∩ benchmark regions |
| RNA fusion | HCC1395 bulk RNA-seq SRR9134735, SRR9134709, SRR9134757 (SRA SRP199641) | 138 validated HCC1395 fusions (Zhang et al., 2016) and RNA-supported SVs (Talsania et al., 2022) | `fusion_truth_match.py` (gene pairs in either orientation, HGNC renaming, intergenic gene lists) |

Helper scripts: `make_figures.py`, `benchmark_plots.py` (Figures 2-5, Tables 1-2), `fusion_truth_match.py`.

---

### 10. Benchmark truth sets

#### Set once
```bash
PIPE=${HOME}/genomics_pipeline
DATA=${HOME}/test_data
R=$PIPE/references
E=$PIPE/envs
T=$R/truth_data
REF=$R/GRCh38/GRCh38.d1.vd1.fa
SDF=$T/GRCh38_sdf
SX="singularity exec -B $PIPE -B $DATA"
BCF="$SX $E/bwa_samtools_fastp_fgbio_bcftools.sif bcftools"
RTG="$SX $E/rtg_tools.sif rtg"
BEDTOOLS="$SX $E/gatk.sif bedtools"
```

#### 10.1 RTG reference for vcfeval [Generate]

```bash
$RTG format -o $R/truth_data/GRCh38_sdf $REF
```
#### 10.2 SEQC2 somatic, HCC1395 [Download + Generate]

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
# evaluation region: padded V6 + UTR targets x SEQC2 high-confidence regions
$BEDTOOLS intersect -a $R/ValidationBed/S07604624_Regions.pad100.bed.gz -b $T/High-Confidence_Regions_v1.2.bed \
  > $R/ValidationBed/S07604624_Regions_pad100_intersection_SEQC2.bed
```

#### 10.3 GIAB HG002 germline, v4.2.1 [Download]

```bash
GB=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
curl -s $GB/ | grep -o 'href="[^"]*benchmark[^"]*"'          # confirm the file names
wget -c -P $R/truth_data $GB/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz $GB/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
                         $GB/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```


#### =================  Benchmark your calls =================

```bash
OUT_S=$DATA/Result_valid_8_Oct_26_som          # <<< outdir of the SEQC2 somatic run
OUT_G=$DATA/Result_valid_8_Oct_26_germ                 # <<< outdir of the HG002 germline run

HEADER="callset\ttype\ttruth\tTP\tFP\tFN\trecall\tprecision\tF1"
```
#### ---------- Counts and metrics per variant type from one vcfeval folder ----------

```bash
count_metrics () {               # usage: count_metrics <vcfeval_dir> <label>
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

####  A. SOMATIC - SEQC2 HCC1395 / HCC1395BL (valid1)

```bash
mkdir -p $OUT_S/benchmark && cd $OUT_S/benchmark
EVAL=$R/ValidationBed/S07604624_Regions_pad100_intersection_SEQC2.bed
awk '{s+=$3-$2} END{print "somatic evaluation region: " s/1e6 " Mb"}' $EVAL
```

#### ---------- A1. Truth: HighConf + MedConf, FILTER = PASS, no samples, split alleles ----------

```bash
$BCF view -G $T/SEQC2_Truth_v1.2.1.vcf.gz \
  | awk 'BEGIN{OFS="\t"} /^#/{print; next} {$7="PASS"; print}' \
  | $BCF norm -a --atom-overlaps . -m -any -f $REF | $BCF sort -Oz -o truth_somatic.vcf.gz -
$BCF index -f -t truth_somatic.vcf.gz
```


#### ---------- A2. Superset and its non-truth labels (FILTER column) ----------

```bash
$BCF concat -a $T/sSNV.MSDUKT.superSet.v1.2.vcf.gz $T/sINDEL.MDKT.superSet.v1.2.vcf.gz \
  | $BCF sort | $BCF norm -a --atom-overlaps . -m -any -f $REF | $BCF sort -Oz -o superset.vcf.gz -

$BCF index -f -t superset.vcf.gz

$BCF query -f '%FILTER\n' superset.vcf.gz | tr ';' '\n' \
  | grep -E '^(HighConf|MedConf|LowConf|Unclassified|REJECT)$' | sort | uniq -c   # label counts

$BCF view -i 'FILTER~"LowConf"'      superset.vcf.gz -Oz -o lowconf.vcf.gz      && $BCF index -f -t lowconf.vcf.gz
$BCF view -i 'FILTER~"Unclassified"' superset.vcf.gz -Oz -o unclassified.vcf.gz && $BCF index -f -t unclassified.vcf.gz
```

#### ---------- A3. Score one call set ----------

#####  Rules: (i) novel calls (not in superset) with tumor VAF < 5% are dropped
#####         (ii) calls at LowConf sites are dropped
#####  main  : calls at Unclassified sites also dropped
#####  strict: calls at Unclassified sites kept -> counted as false positives

```bash
score_somatic () {               # usage: score_somatic <label> <vcf> [bcftools view selection]
  local L=$1 V=$2; shift 2
  local S=$($BCF query -l $V | grep -iE '(_T|TUMOR)$' | head -1)
  [ -n "$S" ] || { echo "ERROR: no tumor column in $V"; return 1; }
  echo "== $L   tumor column: $S"

  $BCF view -s "$S" "$@" $V | $BCF norm -a --atom-overlaps . -m -any -f $REF | $BCF sort -Oz -o $L.calls.vcf.gz -
  $BCF index -f -t $L.calls.vcf.gz

  # (i) novel calls: keep only tumor VAF >= 5%
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

#### ---------- A4. Call sets ----------

```bash
find $OUT_S -name "valid1*.vcf.gz" | sort          # list, then set the paths below
CONS=$(find $OUT_S -name "valid1*_consensus.vcf.gz" | head -1)
FILT=$(find $OUT_S -name "valid1*_filtered.vcf.gz"  | head -1)
#   Per-caller VCFs (optional). If not in the outdir, take them from the work directory.
M2=$OUT_S/05_variant_calling/somatic/snv/mutect2/valid1_T_mutect2.vcf.gz                                  # <<<
STRELKA_SNV=$OUT_S/05_variant_calling/somatic/snv/strelka2/valid1_T_strelka2_snvs.af.vcf.gz                                   # <<<
STRELKA_INDEL=$OUT_S/05_variant_calling/somatic/snv/strelka2/valid1_T_strelka2_indels.af.vcf.gz                              # <<<
MUSE=$OUT_S/05_variant_calling/somatic/snv/muse/valid1_T_muse.af.vcf.gz                                                 # <<<
VARSCAN_SNV=$OUT_S/05_variant_calling/somatic/snv/varscan/valid1_T_varscan_snv.af.vcf.gz                                  # <<<
VARSCAN_INDEL=$OUT_S/05_variant_calling/somatic/snv/varscan/valid1_T_varscan_indel.af.vcf.gz                              # <<<
```

#### ---------- A5. Final call set (Table 1) and consensus ----------

```bash
score_somatic final_pass     $FILT -f PASS     # consensus + VAF/alt reads/depth/normal filters
score_somatic consensus_pass $CONS -f PASS     # >= 2 callers
score_somatic consensus_all  $CONS             # all Mutect2 PASS calls (LOWCONSENSUS included)
```

#### ---------- A6. Each caller before (raw) and after (pass) its own filter (optional) ----------

```bash
$BCF concat -a $STRELKA_SNV $STRELKA_INDEL | $BCF sort -Oz -o strelka2.vcf.gz - && $BCF index -f -t strelka2.vcf.gz
$BCF concat -a $VARSCAN_SNV $VARSCAN_INDEL | $BCF sort -Oz -o varscan.vcf.gz -  && $BCF index -f -t varscan.vcf.gz
score_somatic mutect2_raw   $M2
score_somatic mutect2_pass  $M2             -f PASS
score_somatic strelka2_raw  strelka2.vcf.gz
score_somatic strelka2_pass strelka2.vcf.gz -f PASS
score_somatic muse_raw      $MUSE
score_somatic muse_pass     $MUSE           -e 'FILTER~"Tier5"'
score_somatic varscan_raw   varscan.vcf.gz  -i 'INFO/SS=2'
score_somatic varscan_pass  varscan.vcf.gz  -f PASS -i 'INFO/SS=2'

column -t somatic_metrics.tsv
```


#### ---------------- B. GERMLINE (GIAB HG002) ----------------

```bash
mkdir -p $OUT_G/benchmark && cd $OUT_G/benchmark

echo "$OUT_G"; ls "$OUT_G" 
CALLS=$(find $OUT_G -iname "*hg002*_filtered.vcf.gz"  | head -1); echo "calls: $CALLS"
THR=$(find   $OUT_G -iname "*hg002*thresholds.bed.gz" | head -1); echo "thresholds: $THR"
```

#### ---------- B1. Callable targets: >= 90% of bases at >= 10X ----------

```bash
zcat $THR | awk 'NR==1{for(i=1;i<=NF;i++) if($i=="10X") c=i; if(!c){print "ERROR: no 10X column" > "/dev/stderr"; exit 1}; next}
                 ($3-$2)>0 && $c/($3-$2)>=0.9 {print $1"\t"$2"\t"$3}' \
  | sort -k1,1 -k2,2n > callable.bed
```

#### ---------- B2. Evaluation region = callable targets inside GIAB confident regions ----------

```bash
$BEDTOOLS intersect -a callable.bed -b $T/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  | sort -k1,1 -k2,2n > callable_giab.bed
$BEDTOOLS merge -i callable_giab.bed > HG002_eval.bed
awk '{s+=$3-$2} END{print "germline evaluation region: " s/1e6 " Mb"}' HG002_eval.bed
```

#### ---------- B3. Calls:  ------------------------------------

```bash
$BCF index -f -t $CALLS 2>/dev/null; cp $CALLS hg002_calls.vcf.gz; cp $CALLS.tbi hg002_calls.vcf.gz.tbi
```
#### ---------- B4. Genotype-aware scoring: before filter (all records) and after filter (PASS) ----------

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

#### ---------------- C. RNA FUSION - HCC1395 RNA-seq vs Talsania et al. 2022 fusion list ----------------

```bash
OUT_R=$DATA/Result_valid_8_Oct_26_som                                    # <<< outdir of the RNA fusion run
SCRIPT=$PIPE/fusion_truth_match.py                       # <<< where you keep the script
TRUTH=$T/Talsania_2022_Add_File9.csv                      # <<< Talsania Additional file 9, saved as CSV
mkdir -p $OUT_R/benchmark && cd $OUT_R/benchmark
```

#### ---------- C1. Calls and gene-name table ----------

```bash
CALLS=$(find $OUT_R -name "*_fusions_annotated.tsv" | head -1); echo "calls: $CALLS"
# HGNC table, to translate renamed genes (download once into the truth folder)
[ -s $T/hgnc_complete_set.txt ] || \
  wget -O $T/hgnc_complete_set.txt https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
head -1 $TRUTH | tr ',' '\n' | grep -n "fusion"                  # check: the column "RNA-seq fusion_genes" must be listed
```

#### ---------- C2. Match calls to the list (gene order, renamed genes, intergenic lists handled) ----------

```bash
python3 $SCRIPT --truth $TRUTH --calls $CALLS --hgnc $T/hgnc_complete_set.txt --out fusion_truth_result.tsv
```

#### ---------- C3. Recovered fusions by consensus status and tier ----------

```bash
echo "== recovered, by consensus verdict"
awk -F'\t' 'NR>1 && $4=="FOUND"{n[$7]++} END{for(k in n) print k"\t"n[k]}' fusion_truth_result.tsv | sort
echo "== recovered, by tier"
awk -F'\t' 'NR>1 && $4=="FOUND"{n[$8]++} END{for(k in n) print k"\t"n[k]}' fusion_truth_result.tsv | sort
echo "== missed, but at least one partner gene appears in the calls"
awk -F'\t' 'NR>1 && $4=="missed" && $9!="-"' fusion_truth_result.tsv | wc -l

```
#### ---------- C4. Calls that are not in the reference list ----------

```bash
total=$(grep -v '^#' $CALLS | tail -n +2 | wc -l)
matched=$(awk -F'\t' 'NR>1 && $4=="FOUND"{print $5}' fusion_truth_result.tsv | sort -u | wc -l)
echo "pipeline fusions: $total   matching a listed fusion: $matched   not in the list: $((total - matched))"
```
---

## Troubleshooting

| Symptom | Cause | Solution |
|---|---|---|
| Exit 140 | LSF run-time limit | Raise `time` for that process in an extra config (see Running PAPIO) |
| Exit 255 in `PREPROCESS_AND_ALIGN_DNA` | fastp stopped on a damaged FASTQ (bwa-mem2 log ends normally) | `gzip -t` both FASTQs; download again (`prefetch` → `vdb-validate` → `fasterq-dump`) |
| `0 total baseline variants` in vcfeval | Truth FILTER is not PASS and `--all-records` was not used | Set FILTER to PASS in the truth VCF (step 10.2 and step 10.3) |
| Many LOWSUPPORT calls | Tumor depth below `filter_min_tumor_dp` (20), or `call_regions_gz` narrower than the Mutect2 shards | Check coverage; use the padded `call_regions_gz` (step 4.4) |
| Params ignored | Values in a `-c` config, or `${HOME}` in the params file | Put them in the params YAML with absolute paths |
| Germline recall far below 0.9 | `germline_scatter_intervals` from another capture kit | Use the shards built for the sample's kit (step 4.3) |

---

## Citation

If you use PAPIO, please cite the protocol and the tools listed in the Containers table. Code and versioned containers: Zenodo DOI [10.5281/zenodo.23176505].

## Contact

Pankaj K. Chauhan - PKChauhan@mdanderson.org

