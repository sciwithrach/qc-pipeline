#!/usr/bin/env nextflow

process scSHC {

    container "${projectDir}/containers/r-clustering.sif"

    input:
    path adata

    output:
    path "scSHC_clusters.csv", emit: csv
    path "scSHC_clusters.png", emit: plot
    path "scSHC_results.rds", emit: results_object

    script:
    """
    scSHC.R ${adata} ${task.cpus}
    """

    stub:
    """
    touch scSHC_clusters.csv
    touch scSHC_clusters.png
    touch scSHC_results.rds
    """
 }
