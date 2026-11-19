# PAPIO reference resources

Which resources each analysis needs, what happens without the optional ones,
and which require registration. "Required" means `main.nf` refuses to start
without it for that mode; "optional" means the dependent step is skipped with a
warning. Run `scripts/check_installation.sh <params.yaml>` to check a setup, and
`scripts/make_reference_manifest.sh <params.yaml>` to record versions and checksums.

All DNA resources are GRCh38 (GDC `GRCh38.d1.vd1`). The RNA arm aligns to the same genome with GENCODE v49
(the gene set of Ensembl 115, which the DNA annotation uses); STAR-Fusion and FusionInspector use the CTAT
library (GENCODE v44). The test bundle uses a mini CTAT library for every RNA step.

## Core (every DNA run)

| Parameter | Resource | Required for | Source | Registration |
|---|---|---|---|---|
| `ref_fasta`, `ref_fai`, `ref_dict` | GRCh38.d1.vd1 FASTA + index + sequence dictionary | all DNA modes | NCI GDC reference files | none |
| `ref_index_dir` | bwa-mem2 index of the same FASTA | mapping | built once: `bwa-mem2 index` | none |
| `interval_bed` | capture-kit target BED | WES (all modes) | capture-kit vendor (e.g. Agilent SureDesign) | vendor account |
| `scatter_intervals` | directory of interval shards (padded targets) | variant calling | built once from the target BED (GATK `SplitIntervals`) | none |
| `vep_cache_dir` (+ `vep_cache_version`) | Ensembl VEP offline cache, release 115 | annotation | Ensembl FTP | none |

## Somatic tumour-normal

| Parameter | Resource | Status | Without it | Source | Registration |
|---|---|---|---|---|---|
| `germline_resource` | gnomAD allele-frequency-only VCF | **required** | – | GATK somatic best-practices bundle (`af-only-gnomad.hg38.vcf.gz`) | none |
| `panel_of_normals` | Mutect2 panel of normals | optional | more artefactual calls | build from your normals (`components_pon.nf`), or GATK `1000g_pon.hg38.vcf.gz` | none |
| `pileup_variants` | common biallelic SNPs | optional | contamination filtering skipped | GATK bundle (`small_exac_common_3.hg38.vcf.gz`) | none |
| `cnvkit_pon` | CNVkit reference `.cnn` | optional | reference built from the matched normal only (noisier copy-number profile) | build from normals sequenced with the **same capture kit** | none |
| `msisensor_list` | microsatellite site list | optional | MSI skipped | built once: `msisensor-pro scan` on the FASTA | none |
| `sequenza_gc_wiggle` | GC-content wiggle (50 bp) | optional | purity/ploidy/LOH skipped | built once: `sequenza-utils gc_wiggle -w 50` | none |
| `call_regions_gz` | bgzipped + tabixed padded targets | optional (WES) | Manta/Strelka2 scan the whole genome | built once from the target BED | none |

## Germline

| Parameter | Resource | Status | Without it | Source | Registration |
|---|---|---|---|---|---|
| `dbsnp` | dbSNP VCF | **required** (HaplotypeCaller; also MuSE in somatic runs) | – | GATK hg38 bundle | none |
| `hapmap`, `omni`, `thousandg_snps`, `mills` | training/truth sites | **required** (HaplotypeCaller; run stops at start) | – | GATK hg38 bundle | none |
| `preprocessed_intervals`, `contig_ploidy_priors`, `gcnv_scatter_intervals` | gCNV cohort model inputs | **required if** `run_gcnv: true` | – | built once (GATK `PreprocessIntervals` etc.); needs ≥ 30 normals | none |
| `roh_af_file` | population AF file (tabixed) | optional (WGS only) | ROH skipped | derived from gnomAD | none |

## Shared QC

| Parameter | Resource | Status | Without it | Source | Registration |
|---|---|---|---|---|---|
| `bait_interval_list`, `target_interval_list` | Picard interval lists | optional | CollectHsMetrics skipped | built once: Picard `BedToIntervalList` | none |
| `verifybamid_db` | VerifyBamID2 SVD resource | optional | contamination QC skipped | shipped in the VerifyBamID2 repository (`resource/`) | none |
| `somalier_sites` | somalier sites VCF | optional | sample-swap/relatedness check skipped | somalier GitHub releases | none |
| `dbsnp`, `mills`, `known_indels` | BQSR known sites | **required if** `bqsr: true` | – | GATK hg38 bundle | none |

## Annotation and classification (DNA)

| Parameter | Resource | Status | Without it | Source | Registration |
|---|---|---|---|---|---|
| `tmb_coding_bed` | MANE Select coding exons (merged BED, Ensembl 115) | **required** for WES TMB (`compute_tmb: true`) | – | built once from the Ensembl GTF (README Appendix A, step 6) | none |
| `clinvar_vcf` | ClinVar VCF (dated release) | optional | no ClinVar fields | NCBI ClinVar FTP | none |
| `cosmic_vcf` | COSMIC Genome Screens Mutant VCF | optional | no COSMIC recurrence | COSMIC | **COSMIC licence** (free for academic use) |
| `gnomad_faf_vcf` | gnomAD v4.1 joint grpmax filtering AF (`fafmax_faf95_max_joint`), MANE regions | recommended | classification uses VEP MAX_AF | Zenodo 10.5281/zenodo.23531816 (CC BY 4.0), or `scripts/prepare_gnomad_faf.sh` for all sites | none (open data) |
| `germline_af_thresholds` | Gene-specific germline BA1/BS1 (TSV: GENE, BA1, BS1) | optional | global 0.05 / 0.01 | ClinGen Criteria Specification Registry (cspec.genome.network) | none |
| `vep_plugins_dir/NMD.pm` | Ensembl VEP NMD plugin | recommended | last-exon check only | github.com/Ensembl/VEP_plugins (match the VEP release) | Apache-2.0 |
| `cosmic_cgc_tsv` | Cancer Gene Census | optional | no gene-role labels in classification | COSMIC | **COSMIC licence** |
| `annotsv_dir` | AnnotSV annotation bundle | optional | SV/CNV annotation and germline SV/CNV report sections skipped | AnnotSV installer | none (GeneHancer annotations, if added, need their own licence) |
| `vep_plugins_dir` | VEP plugin data | optional | plugin fields absent | per plugin | varies by plugin |

## RNA fusion

| Parameter | Resource | Status | Without it | Source | Registration |
|---|---|---|---|---|---|
| `ctat_lib` | CTAT genome library, GENCODE v44 (plug-n-play) | **required if** STAR-Fusion, FusionInspector or `rna_use_ctat_reference: true` is used, or `star_index` is not set | – | Broad CTAT resource library (`GRCh38_gencode_v44_CTAT_lib_Oct292023`) | none |
| `ref_fasta`, `ref_fai`, `gtf` | `GRCh38.d1.vd1` and the GENCODE v49 GTF the STAR index was built from (taken from the CTAT library when `rna_use_ctat_reference: true`) | **required** | – | GDC genome and GENCODE release 49 (`setup.sh --full`) | none |
| `star_index` | STAR index of `GRCh38.d1.vd1` with GENCODE v49 | optional | the CTAT library's own index is used | built by `setup.sh --full` (about 40 GB RAM) or README Appendix A, step 8.5 | none |
| `fusioncatcher_db` | FusionCatcher human database (`human_v102`) | **required if** `use_fusioncatcher: true` | – | FusionCatcher SourceForge | none |
| `arriba_blacklist`, `arriba_known_fusions`, `arriba_protein_domains` | Arriba v2.5.1 database files (hg38) | optional (blacklist strongly recommended) | Arriba false positives increase sharply without the blacklist | Arriba release tarball (`database/`) | none |
| `known_fusions_tsv` | known oncogenic fusion pairs | optional | no single-caller known-fusion rescue | Mitelman / COSMIC Fusions / FusionGDB | COSMIC Fusions needs the COSMIC licence |
| `cosmic_cgc_tsv` | Cancer Gene Census | optional | no gene-role fusion tiering | COSMIC | **COSMIC licence** |
| `refflat`, `rrna_interval_list` | RNA QC files from the same GTF (GENCODE v49) | optional | RNA-seq metrics / %rRNA skipped | built by `setup.sh --full` or README Appendix A, step 8.6 | none |

## Minimum sets

- **Somatic WES:** core + `germline_resource` + `dbsnp` (MuSE, on by default) + `tmb_coding_bed` (TMB, on by default). Everything else adds a module.
- **Germline WES:** core + `dbsnp`, HapMap, Omni, 1000G SNPs, Mills (HaplotypeCaller filtering).
- **RNA fusion:** `ctat_lib` (+ `fusioncatcher_db` unless `use_fusioncatcher: false`).
- **Licensed:** only COSMIC files (`cosmic_vcf`, `cosmic_cgc_tsv`, COSMIC Fusions). The pipeline runs without them; classification and fusion tiering lose gene-role information.
