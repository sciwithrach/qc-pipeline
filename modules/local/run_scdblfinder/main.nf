#!/usr/bin/env nextflow

/*
 * scDblFinder
 */

 process run_scDblFinder {

    container "${projectDir}/containers/r-clustering.sif"

    input:
    path adata
    val batch

    output:
    path 'doublets_scdblfinder.csv', emit : csv
    path 'sessionInfo_doublets.txt', emit : session

    script:
    """
    run_scDblFinder.R ${adata} ${batch} ${task.cpus}
    """

    stub:
    """
    touch doublets_scdblfinder.csv
    touch sessionInfo_doublets.txt
    """
 }
