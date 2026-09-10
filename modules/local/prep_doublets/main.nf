#!/usr/bin/env nextflow

/*
 * vaeda
 */

 process prep_doublets {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path adata

    output:
    path "adata_scDblFinder_in.h5ad", emit : scDblFinder
    path "adata_vaeda_in.h5ad", emit : vaeda

    script:
    """
    prep_doublets.py --adata ${adata}
    """
 }
