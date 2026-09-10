#!/usr/bin/env nextflow

process report_stats {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path adata
    val descriptor
    val batch

    output:
    path "report_stats_${descriptor}.csv",       emit: stats
    path "report_stats_${descriptor}_batch.csv", emit: batch_stats, optional: true

    script:
    def batch_arg = batch ? "--batch ${batch}" : ''
    """
    report_stats.py --adata ${adata} --descriptor ${descriptor} ${batch_arg}
    """
}
