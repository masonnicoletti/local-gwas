process SPLIT_PANEL {
    input:
    path panel

    output:
    path "*.txt"

    script:
    """
    awk 'NR>1 && NF>=3 {print \$1 > (\$3 ".txt")}' ${panel}
    """
}

process COUNT_SAMPLES {
    input:
    path ids

    output:
    path "${ids.baseName}.count"

    script:
    """
    awk -v g=${ids.baseName} 'END {print g "\\t" NR}' ${ids} > ${ids.baseName}.count
    """
}

process SUMMARIZE {
    input:
    path counts

    output:
    path "summary.tsv"

    script:
    """
    cat ${counts} | sort > summary.tsv
    """
}

workflow {
    main:
    panel_ch   = channel.fromPath(params.panel)
    groups_ch  = SPLIT_PANEL(panel_ch).flatten()
    counts_ch  = COUNT_SAMPLES(groups_ch)
    summary_ch = SUMMARIZE(counts_ch.collect())

    publish:
    summary = summary_ch
}

output {
    summary {
        path 'summary'
    }
}