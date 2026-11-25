process VERSIONS {
    output:
    path "versions.txt", emit: versions

    script:
    """
    {
      plink2 --version
      plink --version | head -1
      python3 --version
      python3 -c "import pandas, numpy, scipy, matplotlib; print('pandas', pandas.__version__, 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'matplotlib', matplotlib.__version__)"
    } > versions.txt
    """
}

process VCF_TO_PGEN {
    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_snps.{pgen,pvar,psam}"), emit: fileset
    path "${id}_snps.log", emit: log

    script:
    """
    plink2 \\
      --vcf ${vcf} \\
      --snps-only just-acgt \\
      --min-alleles 2 \\
      --max-alleles 2 \\
      --set-all-var-ids '@:#:\$r:\$a' \\
      --make-pgen \\
      --out ${id}_snps
    """
}

process VARIANT_QC {
    input:
    tuple val(id), path(snps)

    output:
    tuple val(id), path("${id}_varqc.{pgen,pvar,psam}"), emit: fileset
    path "${id}_varqc.log", emit: log

    script:
    """
    plink2 \\
      --pfile ${id}_snps \\
      --geno ${params.geno} \\
      --maf ${params.maf} \\
      --make-pgen \\
      --out ${id}_varqc
    """
}

process SAMPLE_QC {
    input:
    tuple val(id), path(varqc)

    output:
    tuple val(id), path("${id}_sampleqc.{pgen,pvar,psam}"), emit: fileset
    path "${id}_unrelated.king.cutoff.{in,out}.id", emit: kinship
    path "*.log", emit: log

    script:
    """
    plink2 \\
      --pfile ${id}_varqc \\
      --king-cutoff ${params.king} \\
      --out ${id}_unrelated

    plink2 \\
      --pfile ${id}_varqc \\
      --keep ${id}_unrelated.king.cutoff.in.id \\
      --make-pgen \\
      --out ${id}_sampleqc
    """
}

process LD_PRUNE {
    input:
    tuple val(id), path(sampleqc)

    output:
    tuple val(id), path("${id}_pruned.prune.{in,out}"), emit: lists
    path "${id}_pruned.log", emit: log

    script:
    """
    plink2 \\
      --pfile ${id}_sampleqc \\
      --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
      --out ${id}_pruned
    """
}

process PCA {
    input:
    tuple val(id), path(sampleqc), path(pruned)

    output:
    tuple val(id), path("${id}_pca.eigenvec"), path("${id}_pca.eigenval"), emit: results
    path "${id}_pca.log", emit: log

    script:
    """
    plink2 \\
      --pfile ${id}_sampleqc \\
      --extract ${id}_pruned.prune.in \\
      --pca ${params.n_pcs} \\
      --out ${id}_pca
    """
}

process ASSOCIATION {
    input:
    tuple val(id), path(sampleqc), path(eigenvec), val(mode)
    path pheno

    output:
    tuple val(id), val(mode), path("${id}_${mode}.*.glm.linear"), emit: results
    path "${id}_${mode}.log", emit: log

    script:
    def model = mode == 'pcs' \
        ? "--covar ${eigenvec} --covar-name PC1-PC${params.n_pcs} --glm hide-covar" \
        : "--glm allow-no-covars"
    """
    plink2 \\
      --pfile ${id}_sampleqc \\
      --pheno ${pheno} \\
      ${model} \\
      --out ${id}_${mode}
    """
}

process ASSOC_REPORT {
    input:
    tuple val(id), path(nocov), path(pcs)
    path key

    output:
    path "${id}_summary.tsv", emit: summary
    path "${id}_causal.tsv", optional: true, emit: causal
    path "${id}_*.png", emit: plots

    script:
    def key_arg = key.name == 'NO_KEY' ? '' : "--key ${key}"
    """
    assoc_report.py \\
      --nocov ${nocov} \\
      --pcs ${pcs} \\
      ${key_arg} \\
      --prefix ${id}
    """
}

workflow {
    main:
    if( !params.pheno )
        error "Please provide a phenotype file: --pheno <file> (or pheno: in a params file)"

    vcf_ch   = channel.fromPath(params.vcf).map { vcf -> tuple(params.id, vcf) }
    pheno_ch = channel.fromPath(params.pheno).first()
    key_path = params.key ?: "${projectDir}/assets/NO_KEY"
    key_ch   = channel.fromPath(key_path).first()
    modes_ch = channel.of('nocov', 'pcs')

    ver     = VERSIONS()
    snps    = VCF_TO_PGEN(vcf_ch)
    varqc   = VARIANT_QC(snps.fileset)
    samples = SAMPLE_QC(varqc.fileset)
    prune   = LD_PRUNE(samples.fileset)
    pca     = PCA(samples.fileset.join(prune.lists))

    assoc_in = samples.fileset
        .join(pca.results)
        .map { id, files, vec, vals -> tuple(id, files, vec) }
        .combine(modes_ch)
    assoc = ASSOCIATION(assoc_in, pheno_ch)

    nocov_ch = assoc.results.filter { id, mode, f -> mode == 'nocov' }.map { id, mode, f -> tuple(id, f) }
    pcs_ch   = assoc.results.filter { id, mode, f -> mode == 'pcs' }.map { id, mode, f -> tuple(id, f) }
    report   = ASSOC_REPORT(nocov_ch.join(pcs_ch), key_ch)

    publish:
    snps_files     = snps.fileset
    varqc_files    = varqc.fileset
    sampleqc_files = samples.fileset
    kinship_lists  = samples.kinship
    prune_lists    = prune.lists
    pca_results    = pca.results
    assoc_results  = assoc.results
    report_tables  = report.summary.mix(report.causal)
    report_plots   = report.plots
    tool_versions  = ver.versions
    plink_logs     = snps.log.mix(varqc.log, samples.log, prune.log, pca.log, assoc.log)
}

output {
    snps_files {
        path 'plink'
    }
    varqc_files {
        path 'plink'
    }
    sampleqc_files {
        path 'plink'
    }
    kinship_lists {
        path 'kinship'
    }
    prune_lists {
        path 'pruning'
    }
    pca_results {
        path 'pca'
    }
    assoc_results {
        path 'assoc'
    }
    report_tables {
        path 'report'
    }
    report_plots {
        path 'report'
    }
    tool_versions {
        path 'report'
    }
    plink_logs {
        path 'logs'
    }
}