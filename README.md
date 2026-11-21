# local-gwas
A basic setup for GWAS workflows

## Starting Notes

### Data

**Genotypes:** 1000 Genomes Project Phase 3
  - 2,504 Samples
  - Reference: GRCh37
  - Chromosome 22 (10 GB)

**Phenotypes:** simulated


### Software

**Tools:**
- PLINK 2.0
- PLINK 1.9
- bcftools
- Python
- R
- Nextflow
- qqman

**Environment:**
- macOS (Apple silicon, arm64)
- conda (Miniforge)


### Phases Overview
1. Setup
2. Prep Data
3. Quality Control
4. Population Structure
5. Phenotype Simulation
6. Association Testing
7. Stress Test
8. Pipeline Setup


### Project layout

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

## Phases

### Phase 0: Setup
- [ ] Install Miniforge, create the conda environment from `envs/gwas.yml`
- [ ] Verify `plink2`, `bcftools`, `nextflow` all run
- [ ] Initialize git repo and project directories

### Phase 1: Look at the data
- [ ] Download chr22 from 1000 Genomes Phase 3
- [ ] Inspect the VCF (header, records, genotype encoding)
- [ ] Convert to PLINK binary format; understand `.bed/.bim/.fam`

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


## Additional Notes

#### Data Options

- 1000 Genomes
- Arabidopsis thaliana

