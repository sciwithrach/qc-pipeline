#!/usr/bin/env nextflow

/*
 * minimal qc filter
 */

 process minimal_filter {

    container "${projectDir}/containers/utricle-qc.sif"

    input :
    path adata
    val descriptor

    output:
    path "adata_minimal_filter_${descriptor}.h5ad", emit : adata

    script:
    """
    minimal_filter.py --adata ${adata} --descriptor ${descriptor}
    """
 }
