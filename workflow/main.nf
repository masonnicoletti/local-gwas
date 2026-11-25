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

workflow {
    main:
    vcf_ch  = channel.fromPath(params.vcf).map { vcf -> tuple(params.id, vcf) }

    snps    = VCF_TO_PGEN(vcf_ch)
    varqc   = VARIANT_QC(snps.fileset)
    samples = SAMPLE_QC(varqc.fileset)

    publish:
    snps_files     = snps.fileset
    varqc_files    = varqc.fileset
    sampleqc_files = samples.fileset
    kinship_lists  = samples.kinship
    plink_logs     = snps.log.mix(varqc.log, samples.log)
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
    plink_logs {
        path 'logs'
    }
}