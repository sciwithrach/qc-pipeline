#!/usr/bin/env nextflow

/*
 * Run cellbender
 */
process run_cellbender {

    container "${projectDir}/containers/cellbender.sif"

    input :
    tuple val(batch), path(split_path), val(total_droplets_included), val(expected_cells)

    output :
    tuple val(batch), path("cellbender_${batch}.h5")              , emit : h5
    tuple val(batch), path("cellbender_${batch}_report.html")     , emit : report
    tuple val(batch), path("cellbender_${batch}.log")             , emit : log
    tuple val(batch), path("cellbender_${batch}.pdf")             , emit : pdf
    tuple val(batch), path("cellbender_${batch}_metrics.csv")     , emit : metrics

    script:
    def droplets_arg = total_droplets_included ? "--total-droplets-included '${total_droplets_included}'" : ""
    def cells_arg    = expected_cells          ? "--expected-cells '${expected_cells}'"                   : ""
    """
    cellbender remove-background \
            --input '${split_path}' \
            --output 'cellbender_${batch}.h5' \
            --cuda \
            --learning-rate 0.00005 \
            ${droplets_arg} \
            ${cells_arg}
    """
}
