#!/usr/bin/env nextflow

process outlier_filter {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path adata
    val clusters
    val doublets
    val descriptor

    output:
    path "plt_mads*${descriptor}.png", emit : madplots
    path "outliers_${descriptor}.csv", emit : csv
    path "adata_final_filter_${descriptor}.h5ad", emit : adata
    path "outliers_${descriptor}.log", emit : log

    script:
    """
    outlier_filter.py --adata ${adata} --clusters ${clusters} --doublets ${doublets} --descriptor ${descriptor}
    """
}
