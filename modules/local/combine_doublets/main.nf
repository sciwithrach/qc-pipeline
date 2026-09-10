#!/usr/bin/env nextflow

/*
 * combine results
 */

 process combine_doublets {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path adata
    path vaeda_csv
    path scdblfinder_csv
    val batch_id

    output:
    path 'adata_doublets.h5ad', emit : adata
    path 'doublets_comparison.png', emit : plot

    script:
    """
    combine_doublets.py --adata ${adata} --vaeda-csv ${vaeda_csv} --scdblfinder-csv ${scdblfinder_csv} --batch ${batch_id}
    """

    stub:
    """
    touch adata_doublets.h5ad
    touch doublets_comparison.png
    """
 }
