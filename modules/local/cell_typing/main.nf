#!/usr/bin/env nextflow

/*
 * simple classifier - adapted from topometry tutorial
 */

 process cell_typing {

    container "${projectDir}/containers/utricle-qc.sif"

    input :
    path adata
    path type_csv
    val descriptor

    output:
    path "simple_cell_types_${descriptor}.csv", emit : csv

    script:
    """
    cell_typing.py --adata ${adata} --type-csv ${type_csv} --descriptor ${descriptor}
    """
 }
