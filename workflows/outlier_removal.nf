/*
* Modules
*
* NB: this is a standalone entry script (its own `workflow {}` / `output {}`),
* not a named sub-workflow included by main.nf — run it directly with:
*   nextflow run workflows/outlier_removal.nf -profile <profile>,apptainer
*/
include { plt_proportion }    from '../modules/local/plt_proportion/main.nf'
include { outlier_filter }    from '../modules/local/outlier_filter/main.nf'
include { filtering_summary } from '../modules/local/filtering_summary/main.nf'
include { report_stats }  from '../modules/local/report_stats/main.nf'
include { qc_plots }      from '../modules/local/qc_plots/main.nf'
include { run_topometry } from '../modules/local/run_topometry/main.nf'
include { cell_typing }   from '../modules/local/cell_typing/main.nf'
include { clustree }         from '../modules/local/clustree/main.nf'
include { scSHC }            from '../modules/local/scshc/main.nf'
include { combine_clusters } from '../modules/local/combine_clusters/main.nf'

/*
* Workflow
*/
workflow {

    main:

    //////////////////////////////////////////////////
    // Outlier removal & reporting
    //////////////////////////////////////////////////

    // outlier filter (returns data in adata.raw if supplied)
    outlier_filter(
        params.adata,
        params.clusters,
        params.doublets,
        params.descriptor
    )

    // qc plots
    qc_plots(
        outlier_filter.out.adata,
        params.batch,
        'outlierFilter'
    )

    // report stats
    report_stats(
        outlier_filter.out.adata,
        'outlierFilter',
        params.batch
    )

    // gather csv paths — mix prior stats CSVs with the new one, then collect into a list
    ch_csvs = channel.fromPath(params.csv_paths)
                .mix(report_stats.out.stats)
                .collect()

    // gather batch csv paths — mix prior batch stats CSVs with the new one
    // (batch_csv_paths is optional; skip fromPath entirely when not supplied)
    ch_prior_batch_csvs = params.batch_csv_paths ? channel.fromPath(params.batch_csv_paths) : channel.empty()
    ch_batch_csvs = ch_prior_batch_csvs
                        .mix(report_stats.out.batch_stats)
                        .collect()
                        .ifEmpty([])

    // filtering_summary
    filtering_summary(
        ch_csvs,
        ch_batch_csvs
    )

    // plot final proportions per batch and cell type
    plt_proportion(
        outlier_filter.out.adata,
        params.batch,
        params.descriptor
    )

    //////////////////////////////////////////////////
    // Generate the final dataset
    //////////////////////////////////////////////////

    // re-run topometry on filtered data
    run_topometry(
        outlier_filter.out.adata,
        'plate',
        'postqc'
    )

    // scSHC clusters (optional — skip with --skip_scSHC)
    if (!params.skip_scSHC) {
        scSHC(outlier_filter.out.adata)
    }
    scSHC_csv = params.skip_scSHC ? channel.value('NO_FILE') : scSHC.out.csv

    // simple cell typing
    cell_typing(
        outlier_filter.out.adata,
        params.celltype_csv,
        'postqc'
    )

    // clustree on topometry
    clustree(
        run_topometry.out.adata,
        cell_typing.out.csv
    )

    // cluster comparison plot
    combine_clusters(
        run_topometry.out.adata,
        scSHC_csv,
        cell_typing.out.csv
    )

    publish:
    // outliers
    madplots                = outlier_filter.out.madplots
    csv_outliers            = outlier_filter.out.csv
    adata_outliers          = outlier_filter.out.adata
    log_outliers            = outlier_filter.out.log
    // qc reporting
    plot_qc                 = qc_plots.out
    plot_stats              = report_stats.out.stats
    plot_stats_batch        = report_stats.out.batch_stats
    // summary plots
    plot_filter_summary     = filtering_summary.out
    plot_proportions        = plt_proportion.out
    // topometry
    plot_topo_pca           = run_topometry.out.plt_pca
    plot_topo_geom           = run_topometry.out.plt_geometry
    plot_topo_embeddings    = run_topometry.out.plt_comparison
    topo_object             = run_topometry.out.tg
    // clustree
    plot_clustree           = clustree.out.plot
    // scSHC
    csv_scSHC  = params.skip_scSHC ? channel.empty() : scSHC.out.csv
    plot_scSHC = params.skip_scSHC ? channel.empty() : scSHC.out.plot
    rds_scSHC  = params.skip_scSHC ? channel.empty() : scSHC.out.results_object
    // combine clusters
    plot_cluster_comparison = combine_clusters.out.plot
    adata_cluster           = combine_clusters.out.adata
}

output {
    // reports
    csv_outliers {
        path 'reports'
    }
    log_outliers {
        path 'reports'
    }
    csv_scSHC {
        path 'reports'
    }
    rds_scSHC {
        path 'reports'
    }

    // objects
    adata_outliers {
        path 'anndatas'
    }
    adata_cluster {
        path 'anndatas'
    }
    topo_object {
        path 'anndatas'
    }
    // plots - qc
    plot_qc {
        path 'plots/qc'
    }
    plot_stats {
        path 'plots/qc'
    }
    plot_stats_batch {
        path 'plots/qc'
    }
    plot_proportions {
        path 'plots/qc'
    }
    madplots {
        path 'plots/qc'
    }
    plot_filter_summary {
        path 'plots/qc'
    }
    // plots - topometry
    plot_topo_pca {
        path 'plots/topo'
    }
    plot_topo_geom {
        path 'plots/topo'
    }
    plot_topo_embeddings {
        path 'plots/topo'
    }
    // plots - clusters
    plot_scSHC {
        path 'plots/clusters'
    }
    plot_clustree {
        path 'plots/clusters'
    }
    plot_cluster_comparison {
        path 'plots/clusters'
    }

}
