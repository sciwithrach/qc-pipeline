#!/usr/bin/env nextflow

/*
 * Get unique batch names from adata
 */
process get_batches {

    container "${projectDir}/containers/utricle-qc.sif"

    input :
    path raw_adata
    val batch_col

    output :
    stdout

    script :
    """
    get_batches.py --adata ${raw_adata} --batch-col ${batch_col}
    """
}
