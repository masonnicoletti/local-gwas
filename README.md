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

### 2. PLINK 1 format conversion (`chr22_snps_bed`)

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
- [ ] Variant filters: call rate, minor allele frequency, Hardy-Weinberg
- [ ] Sample filters: missingness, heterozygosity, relatedness
- [ ] Record how many variants/samples each filter removes, and why

### Phase 3: Population structure
- [ ] LD pruning
- [ ] PCA
- [ ] Plot PC1 vs PC2 colored by population

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

