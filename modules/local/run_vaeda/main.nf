#!/usr/bin/env nextflow

process run_vaeda {

    container "${projectDir}/containers/vaeda.sif"

    input:
    path adata
    val batch

    output:
    path 'doublets_vaeda.csv', emit : csv

    script:
    """
    run_vaeda.py --adata ${adata} --batch ${batch}
    """
 }
