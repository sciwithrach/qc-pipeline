#!/usr/bin/env nextflow

/*
* utricleqc — primary QC pipeline entry point.
*
* A second, independent workflow (outlier removal) lives at
* workflows/outlier_removal.nf and is run directly, e.g.:
*   nextflow run workflows/outlier_removal.nf -profile alldata,apptainer
* It is not wired through this entry point — see README.md.
*/

include { QC_PIPELINE } from './workflows/qc_pipeline.nf'

workflow {

    main:
    QC_PIPELINE()

    publish:

    /*
    * Cell calling
    */

    // barcode rank plot
    plot_barcode_rank    = QC_PIPELINE.out.plot_barcode_rank

    // adatas
    adata_joined_batches = QC_PIPELINE.out.adata_joined_batches

    // all cellbender outputs
    cellbender_h5       = QC_PIPELINE.out.cellbender_h5
    cellbender_report   = QC_PIPELINE.out.cellbender_report
    cellbender_log      = QC_PIPELINE.out.cellbender_log
    cellbender_pdf      = QC_PIPELINE.out.cellbender_pdf
    cellbender_metrics  = QC_PIPELINE.out.cellbender_metrics

    /*
    * Doublets
    */

    // adatas
    adata_doublets = QC_PIPELINE.out.adata_doublets

    // csvs
    csv_vaeda = QC_PIPELINE.out.csv_vaeda
    csv_scdblfinder = QC_PIPELINE.out.csv_scdblfinder

    // plots
    plot_doublets_comparison = QC_PIPELINE.out.plot_doublets_comparison

    // session info
    sessionInfo_scdblfinder = QC_PIPELINE.out.sessionInfo_scdblfinder

    /*
    * Dimensionality reduction
    */
    plot_topo_pca = QC_PIPELINE.out.plot_topo_pca
    plot_topo_geom = QC_PIPELINE.out.plot_topo_geom
    plot_topo_embeddings = QC_PIPELINE.out.plot_topo_embeddings

    /*
    * Clustering
    */

    // clustree
    plot_clustree = QC_PIPELINE.out.plot_clustree

    // scSHC
    csv_scSHC  = QC_PIPELINE.out.csv_scSHC
    plot_scSHC = QC_PIPELINE.out.plot_scSHC
    rds_scSHC  = QC_PIPELINE.out.rds_scSHC

    // combine clusters
    plot_cluster_comparison = QC_PIPELINE.out.plot_cluster_comparison

    /*
    * QC reporting
    */

    // plots
    qcplots_unfiltered = QC_PIPELINE.out.qcplots_unfiltered
    qcplots_cellcalling = QC_PIPELINE.out.qcplots_cellcalling
    qcplots_minFilter = QC_PIPELINE.out.qcplots_minFilter

    // stats
    stats_unfiltered        = QC_PIPELINE.out.stats_unfiltered
    stats_cellCalling       = QC_PIPELINE.out.stats_cellCalling
    stats_minFilter         = QC_PIPELINE.out.stats_minFilter
    batch_stats_unfiltered  = QC_PIPELINE.out.batch_stats_unfiltered
    batch_stats_cellCalling = QC_PIPELINE.out.batch_stats_cellCalling
    batch_stats_minFilter   = QC_PIPELINE.out.batch_stats_minFilter

}

output {

    /*
    * NB: Final adata for outliers etc is combine_doublets.out.adata
    */

    // adatas
    // raw - cell calling labelled
    adata_joined_batches {
        path 'anndatas'
        mode 'copy'
    }
    // final output
    adata_doublets {
        path 'anndatas'
        mode 'copy'
    }

    // plots
    plot_barcode_rank {
        path 'plots/cellcalling'
        mode 'copy'
    }
    plot_doublets_comparison {
        path 'plots/doublets'
        mode 'copy'
    }
    plot_cluster_comparison {
        path 'plots/clusters'
        mode 'copy'
    }
    plot_clustree {
        path 'plots/clusters'
        mode 'copy'
    }
    plot_scSHC {
        path 'plots/clusters'
        mode 'copy'
    }

    // cellbender
    cellbender_h5 {
        path 'anndatas/cellbender'
        mode 'copy'
    }
    cellbender_report {
        path 'reports/cellbender'
        mode 'copy'
    }
    cellbender_log {
        path 'reports/cellbender'
        mode 'copy'
    }
    cellbender_pdf {
        path 'reports/cellbender'
        mode 'copy'
    }
    cellbender_metrics {
        path 'reports/cellbender'
        mode 'copy'
    }

    // filter reports
    stats_unfiltered {
        path 'reports/stats'
        mode 'copy'
    }
    stats_cellCalling {
        path 'reports/stats'
        mode 'copy'
    }
    stats_minFilter {
        path 'reports/stats'
        mode 'copy'
    }
    batch_stats_unfiltered {
        path 'reports/stats'
        mode 'copy'
    }
    batch_stats_cellCalling {
        path 'reports/stats'
        mode 'copy'
    }
    batch_stats_minFilter {
        path 'reports/stats'
        mode 'copy'
    }

    // other reports
    sessionInfo_scdblfinder {
        path 'reports/doublets'
        mode 'copy'
    }
    rds_scSHC {
        path 'reports/clusters'
        mode 'copy'
    }

    // qc plots
    qcplots_unfiltered {
        path 'plots/filtering'
        mode 'copy'
    }
    qcplots_cellcalling {
        path 'plots/filtering'
        mode 'copy'
    }
    qcplots_minFilter {
        path 'plots/filtering'
        mode 'copy'
    }

    // topometry
    plot_topo_pca {
        path 'plots/topometry'
        mode 'copy'
    }
    plot_topo_geom {
        path 'plots/topometry'
        mode 'copy'
    }
    plot_topo_embeddings {
        path 'plots/topometry'
        mode 'copy'
    }

    // csvs
    csv_vaeda {
        path 'csv'
        mode 'copy'
    }
    csv_scdblfinder {
        path 'csv'
        mode 'copy'
    }
    csv_scSHC {
        path 'csv'
        mode 'copy'
    }

}
