#!/usr/bin/env nextflow

process plt_proportion {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path adata
    val batch
    val descriptor

    output:
    path "plt_proportion_*${descriptor}.png"

    script:
    """
    plt_proportion.py --adata ${adata} --batch ${batch} --descriptor ${descriptor}
    """

}
