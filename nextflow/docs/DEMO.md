# Testing PAPIO

Two separate paths:

| | **Test run** | **Full benchmark** |
|---|---|---|
| Answers | Is PAPIO installed correctly, and does every mode run end to end? | Does PAPIO reproduce the manuscript's accuracy? |
| Data | demo reads (truth-variant and microsatellite targets of one chromosome for DNA; fusion partner genes for RNA) | full SEQC2 WES pair + HCC1395 RNA-seq |
| References | test bundle (one chromosome cut to the demo windows, mini CTAT library), tens of MB; setup builds its bwa-mem2 and STAR indexes in a few minutes | full GRCh38 resources |
| Hardware | one Linux server, 16 CPUs, 64 GB RAM, no scheduler | HPC (LSF / SLURM) or a large server |
| Setup | `bash setup.sh --test <dir>` | `bash setup.sh --full <dir>` |
| Check | `scripts/compare_to_expected.sh` | `rtg vcfeval` against SEQC2 truth |

The test run checks **execution, not accuracy**: its calls are not
comparable with the manuscript's benchmark numbers.

---

## 1. Test run (users and reviewers)

```bash
git clone https://github.com/chauhan2987/PAPIO.git genomics_pipeline
cd genomics_pipeline/nextflow                 # every command in this file runs from here
bash setup.sh --test ~/papio_test              # containers + test bundle + params files

bash scripts/check_installation.sh ~/papio_test/params_test_somatic.yaml

nextflow run main.nf -params-file ~/papio_test/params_test_somatic.yaml  -profile local
nextflow run main.nf -params-file ~/papio_test/params_test_germline.yaml -profile local
nextflow run main.nf -params-file ~/papio_test/params_test_rna.yaml      -profile local

bash scripts/compare_to_expected.sh check ~/papio_test/runs/test_somatic/Result  tests/expected/somatic
bash scripts/compare_to_expected.sh check ~/papio_test/runs/test_germline/Result tests/expected/germline
bash scripts/compare_to_expected.sh check ~/papio_test/runs/test_rna/Result      tests/expected/rna
```

`-profile slurm` on a SLURM cluster; no profile = LSF. Edit the CPU/memory
numbers in `conf/local.config` to your machine.

Switched off in the test run, and why:

| Step | Reason |
|---|---|
| Sequenza, TMB, mutational signatures | need the whole exome; one chromosome gives wrong answers |
| (CNVkit runs) | it uses the matched normal; on one chromosome its copy-number profile is not meaningful |
| COSMIC annotation, CGC gene roles | COSMIC is licensed and cannot be redistributed |
| AnnotSV SV/CNV annotation | its annotation database is a separate large install |
| FusionCatcher | its database (~19 GB) has no small version |
| Arriba blacklist / known fusions / protein domains | GRCh38 coordinates do not fit the mini RNA genome |

Expected warnings: QC_GATE coverage/on-target warnings (DNA), and the
`check_installation.sh` WARN lines for the resources above.

`compare_to_expected.sh` passes when every expected output file exists,
variant/row counts are within 5% (or 2 records), and the set of PASS
fusions is identical.

---

## 2. Full benchmark

```bash
bash setup.sh --full ~/papio_full                     # public GRCh38 resources + indexes + params files
bash scripts/download_demo_fastq.sh ~/papio_full/demo  # SEQC2 FASTQs (MD5-verified) + samplesheets (path used by the params files)
# add the capture-kit BED/interval lists and (optionally) COSMIC, then:
bash scripts/check_installation.sh ~/papio_full/params/params_somatic.yaml
nextflow run main.nf -params-file ~/papio_full/params/params_somatic.yaml
```

---

## 3. Building the test bundle (developer, once per release)

```bash
# 1. full benchmark run finished (section 2)
# 2. demo reads: picks the smallest chromosome with ≥ 50 truth SNVs and ≥ 5 truth indels, keeps
#    the capture targets holding a truth variant + 100 microsatellite targets (demo windows) and
#    caps the RNA at 100,000 non-chimeric pairs (all chimeric pairs kept)
bash scripts/make_demo_subset.sh
# 3. test bundle: cuts resources, VEP cache and MSI list to the demo windows, builds targets +
#    scatter + mini CTAT library; ships no index (setup.sh --test builds them); fails if any
#    file is > 95 MB or the tarball > 150 MB (GitHub / Zenodo friendly)
bash scripts/make_test_bundle.sh
# 4. try it first:  bash setup.sh --test ~/papio_test --bundle <papio_test_bundle_*.tar.gz>
#    then upload the tarball to Zenodo; set TEST_BUNDLE_URL and TEST_BUNDLE_MD5 in setup.sh
# 5. run section 1 yourself, then record the expected outputs:
bash scripts/compare_to_expected.sh record ~/papio_test/runs/test_somatic/Result  tests/expected/somatic
bash scripts/compare_to_expected.sh record ~/papio_test/runs/test_germline/Result tests/expected/germline
bash scripts/compare_to_expected.sh record ~/papio_test/runs/test_rna/Result      tests/expected/rna
# 6. commit tests/expected/ and the three pipeline_trace.txt files (demo compute requirements)
```
