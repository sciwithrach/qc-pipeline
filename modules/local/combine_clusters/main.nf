#!/usr/bin/env nextflow

process combine_clusters {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path adata
    val scSHC_csv
    path celltype_csv

    output:
    path 'cluster_comparison_on_embedding.png', emit : plot
    path "adata_clusters.h5ad", emit : adata

    script:
    """
    combine_clusters.py --adata ${adata} --celltype-csv ${celltype_csv} --scshc-csv ${scSHC_csv}
    """
 }
