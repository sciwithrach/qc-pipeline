#!/usr/bin/env nextflow

process qc_plots {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path adata
    val batch
    val descriptor

    output:
    path "qcplot_*${descriptor}.png"

    script:
    """
    qc_plots.py --adata ${adata} --batch ${batch} --descriptor ${descriptor}
    """
}
