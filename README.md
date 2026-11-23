# local-gwas
A basic setup for GWAS workflows

---

## Starting Notes

#### Data

**Genotypes:** 1000 Genomes Project Phase 3
  - 2,504 Samples
  - Reference: GRCh37
  - Chromosome 22 (10 GB)

**Phenotypes:** simulated


#### Software

**Tools:**
- PLINK 2.0
- PLINK 1.9
- bcftools
- Python
- R
- Nextflow
- qqman

**Environment:**
- MacOS (Apple silicon, arm64)
- Conda (Miniforge)


#### Phases Overview
1. Setup
2. Prep Data
3. Quality Control
4. Population Structure
5. Phenotype Simulation
6. Association Testing
7. Stress Test
8. Pipeline Setup

----

## Project layout

```
gwas-project/
├── data/
│   ├── raw/          # downloaded files
│   └── processed/    # PLINK files after processing
├── results/          # plots, association output
├── scripts/          # standalone scripts for manual workflow
├── workflow/         # Nextflow pipeline
│   ├── main.nf
│   ├── nextflow.config
│   └── modules/
├── envs/gwas.yml     # conda environment
├── logs/
└── README.md
```

## Raw data (`data/raw/`)

Source: 1000 Genomes Phase 3, GRCh37 (`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`)

Retrieved: 2026-11-20

| Title | File name | Description | Size |
|---|---|---|---|
| chr22 genotypes | `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | Phased genotypes for all 2,504 samples on chromosome 22 as compressed VCF | 197 MB |
| chr22 genotype index | `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi` | Tabix index for bcftools efficiency | 36 KB |
| Sample panel | `integrated_call_samples_v3.20130502.ALL.panel` | Maps each sample to its population, super-population, and sex | 55 KB |
| Related individuals | `20140625_related_individuals.txt` | Samples flagged as related to others in the release | 1.2 KB |

## Processed data (`data/processed/`)

Input: raw chr22 VCF 

Generated: 2026-11-20 

Tool: PLINK v2.0.0-b.1.4

All filtered datasets contain **2,504 samples** and **1,055,454 variants** (biallelic SNPs on chr22).

### 1. PLINK 2 conversion (`chr22_snps`)

Created with:
```bash
plink2 --vcf $VCF --snps-only just-acgt --min-alleles 2 --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' --make-pgen --out data/processed/chr22_snps
```
Of 1,103,547 variants in the VCF, 48,093 were excluded (multiallelic sites, indels, structural variants, and non-ACGT alleles). Variant IDs follow the format `chr:pos:ref:alt`.

| Title | File name | Description | Size |
|---|---|---|---|
| Genotypes | `chr22_snps.pgen` | PLINK 2 binary genotype matrix (samples × variants) | 76 MB |
| Variant info | `chr22_snps.pvar` | Text table of variants: chromosome, position, ID, REF, ALT, INFO | 156 MB |
| Sample info | `chr22_snps.psam` | Text table of 2,504 sample IDs | 27 KB |
| Conversion log | `chr22_snps.log` | PLINK's record of the run | 1.3 KB |

### 2. PLINK 1 conversion (`chr22_snps_bed`)

Created with:
```bash
plink2 --pfile data/processed/chr22_snps --make-bed --out data/processed/chr22_snps_bed
```
Same data as above in the older `.bed/.bim/.fam` format, kept for comparison and for tools that require it.

| Title | File name | Description | Size |
|---|---|---|---|
| Genotypes | `chr22_snps_bed.bed` | PLINK 1 binary genotypes, fixed 2 bits per genotype (SNP-major) | 631 MB |
| Variant info | `chr22_snps_bed.bim` | Text table of variants: chromosome, ID, genetic distance, position, allele 1, allele 2 | 35 MB |
| Sample info | `chr22_snps_bed.fam` | Text table of samples: family ID, individual ID, parental IDs, sex, phenotype | 47 KB |
| Conversion log | `chr22_snps_bed.log` | PLINK's record of the format conversion | 822 B |

### 3. Variant QC (`chr22_varqc`)

Created with:
```bash
plink2 --pfile data/processed/chr22_snps --geno 0.02 --maf 0.01 \
  --make-pgen --out data/processed/chr22_varqc
```

| Step | Filter | Removed | Remaining |
|---|---|---|---|
| Start | biallelic SNPs (`chr22_snps`) | n/a | 1,055,454 |
| Missingness | `--geno 0.02` | 0 | 1,055,454 |
| Minor allele frequency | `--maf 0.01` | 884,505 | 170,949 |

| Title | File name | Description | Size |
|---|---|---|---|
| Genotypes | `chr22_varqc.pgen` | PLINK 2 genotypes after variant QC | 58 MB |
| Variant info | `chr22_varqc.pvar` | Variant table for the retained variants | 28 MB |
| Sample info | `chr22_varqc.psam` | Same 2,504 samples as before | 27 KB |
| QC log | `chr22_varqc.log` | PLINK's record of the filter run | 1.1 KB |

**QC decisions**
- **Missingness:** `--geno 0.02` removed nothing. The data are imputed and complete, so this filter is a safeguard here.
- **MAF ≥ 1%:** rare variants have almost no power and unstable estimates, and a variant with one carrier is driven by one person. 83.8% of variants were rare.
- **HWE not applied to the pooled sample.** 24.5% of retained variants fail at p < 1e-6 in the pooled data versus 0.9% in the European-only subset. The difference is mostly population structure (Wahlund effect, a heterozygote deficit from mixing populations), not genotyping error.

### 4. Kinship files (`chr22_unrelated`)

Created with:
```bash
plink2 --pfile data/processed/chr22_varqc \
  --king-cutoff 0.0884 \
  --out data/processed/chr22_unrelated
```

| Title | File name | Description | Size |
|---|---|---|---|
| Kept samples | `chr22_unrelated.king.cutoff.in.id` | Samples kept after removing relatives at kinship ≥ 0.0884 | 19 KB |
| Removed samples | `chr22_unrelated.king.cutoff.out.id` | Samples removed for relatedness | 661 B |
| Kinship log | `chr22_unrelated.log` | PLINK's record of the `--king-cutoff` run | 1.1 KB |

**Kinship filtering:** 82 removed samples

### 5. Sample-QC dataset (`chr22_sampleqc`)

Created with:
```bash
plink2 --pfile data/processed/chr22_varqc \
  --keep data/processed/chr22_unrelated.king.cutoff.in.id \
  --make-pgen --out data/processed/chr22_sampleqc
```

| Title | File name | Description | Size |
|---|---|---|---|
| Genotypes | `chr22_sampleqc.pgen` | PLINK 2 genotypes after variant and sample QC | 56 MB |
| Variant info | `chr22_sampleqc.pvar` | Variant table (unchanged from `chr22_varqc`, 170,949 variants) | 28 MB |
| Sample info | `chr22_sampleqc.psam` | Sample table for the 2,422 retained samples | 27 KB |
| QC log | `chr22_sampleqc.log` | PLINK's record of the run | 1013 B |

### 6. LD pruning (`chr22_pruned`)

Created with:
```bash
plink2 --pfile data/processed/chr22_sampleqc \
  --indep-pairwise 1000 50 0.2 \
  --out data/processed/chr22_pruned
```

| Title | File name | Description | Size |
|---|---|---|---|
| Kept variants | `chr22_pruned.prune.in` | IDs of the 14,744 variants retained after LD pruning | 231 KB |
| Removed variants | `chr22_pruned.prune.out` | IDs of the 156,205 variants removed for high LD | 2.4 MB |
| Pruning log | `chr22_pruned.log` | PLINK's record of the pruning run | 908 B |

**Why prune:** nearby variants are correlated (linkage disequilibrium), so a block of correlated variants would cause LD structure to dominate PCA instead of ancestry. Pruning keeps one representative per correlated cluster.

**Result:** 14,744 variants kept (8.6%) from 170,949. Pruning writes variant lists only; it does not create a new genotype dataset.


## Phases

### Phase 0: Setup
- [x] Install Miniforge, create the conda environment from `envs/gwas.yml`
- [x] Verify `plink2`, `bcftools`, `nextflow` all run
- [x] Initialize git repo and project directories

### Phase 1: Look at the data
- [x] Download chr22 from 1000 Genomes Phase 3
- [x] Inspect the VCF (header, records, genotype encoding)
- [x] Convert to PLINK binary format; understand `.bed/.bim/.fam`

### Phase 2: Variant and sample QC
- [x] Variant filters: call rate, minor allele frequency, Hardy-Weinberg
- [x] Sample filters: missingness, heterozygosity, relatedness
- [x] Record how many variants/samples each filter removes, and why

### Phase 3: Population structure
- [x] LD pruning
- [x] PCA
- [x] Plot PC1 vs PC2 colored by population

### Phase 4: Simulate the phenotype
- [ ] Write a Python script: choose causal SNPs, effect sizes, noise
- [ ] Save the true causal SNPs separately as the answer key

### Phase 5: Association testing
- [ ] Run `plink2 --glm` with and without PC covariates
- [ ] Manhattan plot, QQ plot, genomic inflation factor (λ)
- [ ] Check whether the planted signals are recovered

### Phase 6: Break it on purpose
- [ ] Build a phenotype confounded with ancestry
- [ ] Run association without PCs; observe λ inflation and false positives
- [ ] Add PCs back; observe the correction

### Phase 7: Port to Nextflow
Convert the working manual steps into a pipeline, one process at a time, checking
each against the manual results.
- [ ] 7a. Nextflow basics: processes, channels, `workflow` block (toy example first)
- [ ] 7b. Port conversion + QC as processes
- [ ] 7c. Port LD pruning + PCA
- [ ] 7d. Port association + plotting
- [ ] 7e. Parameterize in `nextflow.config` (input paths, thresholds, covariate on/off)
- [ ] 7f. Verify the pipeline reproduces the manual results; try `-resume`

### Phase 8 (optional): Extensions
- LD clumping and gene annotation of top hits
- Mixed-model association (REGENIE)
- Heritability (LDSC), polygenic scores
- Real data: public summary statistics (GWAS Catalog, Pan-UK Biobank)

---

## Additional Notes

#### Data Options

- 1000 Genomes
- Arabidopsis thaliana

