#!/usr/bin/env nextflow

process filtering_summary {

    container "${projectDir}/containers/utricle-qc.sif"

    input:
    path csv_list
    path batch_csv_list

    output:
    path "plt_filter*.png"

    script:
    def batch_arg = batch_csv_list ? "--batch_csvs ${batch_csv_list}" : ''
    """
    filtering_summary.py ${csv_list} ${batch_arg}
    """
}
