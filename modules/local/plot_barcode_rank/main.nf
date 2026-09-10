#!/usr/bin/env nextflow

/*
 * Barcode rank plot
 */
process plot_barcode_rank {

    container "${projectDir}/containers/utricle-qc.sif"

    input :
    path raw_adata
    val batch

    output :
    path "plot_barcode_rank_${batch}.png", emit : plot

    script :
    """
    plot_barcode_rank.py --adata ${raw_adata} --batch ${batch}
    """
}
