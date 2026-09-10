/*
* Modules
*/
include { plot_barcode_rank } from '../modules/local/plot_barcode_rank/main.nf'
include { get_batches }       from '../modules/local/get_batches/main.nf'
include { split_adata }       from '../modules/local/split_adata/main.nf'
include { run_cellbender }    from '../modules/local/run_cellbender/main.nf'
include { join_adatas }       from '../modules/local/join_adatas/main.nf'

include { prep_doublets }    from '../modules/local/prep_doublets/main.nf'
include { run_vaeda }        from '../modules/local/run_vaeda/main.nf'
include { run_scDblFinder }  from '../modules/local/run_scdblfinder/main.nf'
include { combine_doublets } from '../modules/local/combine_doublets/main.nf'

include { qc_plots as qc_unfiltered; qc_plots as qc_cellCalling; qc_plots as qc_minFilter } from '../modules/local/qc_plots/main.nf'
include { report_stats as report_unfiltered; report_stats as report_cellCalling; report_stats as report_minFilter } from '../modules/local/report_stats/main.nf'

include { cell_typing }    from '../modules/local/cell_typing/main.nf'
include { minimal_filter } from '../modules/local/minimal_filter/main.nf'
include { run_topometry }  from '../modules/local/run_topometry/main.nf'

include { clustree }         from '../modules/local/clustree/main.nf'
include { scSHC }            from '../modules/local/scshc/main.nf'
include { combine_clusters } from '../modules/local/combine_clusters/main.nf'

/*
* Workflow
*/
workflow QC_PIPELINE {

    main:

    /*
    * Cell calling
    */

    // qcplots - unfiltered
    if (!params.skip_cellbender) {
        qc_unfiltered(
            params.raw_adata,
            params.batch,
            'unfiltered'
        )
    }

    // report stats - unfiltered
    if (!params.skip_cellbender) {
        report_unfiltered(
            params.raw_adata,
            'unfiltered',
            params.batch
        )
    }

    // barcode_rank plot
    if (!params.skip_cellbender) {
        plot_barcode_rank(
            params.raw_adata,
            params.batch
        )
    }

    // get batch names as a channel
    if (!params.skip_cellbender) {
        ch_batches = get_batches(params.raw_adata, params.batch)
                        .splitText()
                        .map { it.trim() }
    }

    // split adata by batch (one process per batch)
    if (!params.skip_cellbender) {
        split_adata(
            params.raw_adata,
            ch_batches,
            params.batch
        )
    }

    // load per-batch cellbender settings, or default to null for all batches
    if (!params.skip_cellbender) {
        ch_cellbender_params = params.cellbender_csv
            ? channel.fromPath(params.cellbender_csv)
                .splitCsv(header: true)
                .map { row -> tuple(row.batch, row.total_droplets_included ?: null, row.expected_cells ?: null) }
            : split_adata.out.split_path.map { batch, path -> tuple(batch, null, null) }
    }

    // run cellbender — (batch, split.h5ad, total_droplets, expected_cells)
    if (!params.skip_cellbender) {
        run_cellbender(split_adata.out.split_path.join(ch_cellbender_params))
    }

    // re-join split paths with cellbender outputs by batch key
    // -> (batch, split.h5ad, cellbender.h5)
    if (!params.skip_cellbender) {
        ch_for_join = split_adata.out.split_path
                                        .join(run_cellbender.out.h5)
    }

    // join raw adata with cellbender metadata — collect all batches into lists
    if (!params.skip_cellbender) {
        join_adatas(
            params.raw_adata,
            ch_for_join.map { batch, split, cb -> cb }.collect()
        )
    }

    // supply either the adata from params or the outputs from the workflow
    ch_cell_calling_adata = params.skip_cellbender ? params.raw_adata : join_adatas.out.filtered

    // qcplots - post cell calling
    qc_cellCalling(
        ch_cell_calling_adata,
        params.batch,
        'cellCalling'
    )

    // report stats - post cell calling
    report_cellCalling(
        ch_cell_calling_adata,
        'cellCalling',
        params.batch
    )

    // simple cell typing
    cell_typing(
        ch_cell_calling_adata,
        params.celltype_csv,
        'qc'
    )

    /*
    * QC minimal filter
    */
    minimal_filter(
        ch_cell_calling_adata,
        'qc'
    )

    // qcplots - post minimal filter
    qc_minFilter(
        minimal_filter.out,
        params.batch,
        'minimalFilter'
    )

    // report stats - post minimal filter
    report_minFilter(
        minimal_filter.out,
        'minimalFilter',
        params.batch
    )

    /*
    * Dimensionality reduction
    */
    run_topometry(
        minimal_filter.out,
        params.metadata,
        'qc'
    )

    /*
    * Clustering & cluster stability
    */
    // clustree on topometry
    clustree(
        run_topometry.out.adata,
        cell_typing.out.csv
    )

    // scSHC clusters (optional — skip with --skip_scSHC)
    if (!params.skip_scSHC) {
        scSHC(minimal_filter.out)
    }
    scSHC_csv = params.skip_scSHC ? channel.value('NO_FILE') : scSHC.out.csv

    // cluster comparison plot
    combine_clusters(
        run_topometry.out.adata,
        scSHC_csv,
        cell_typing.out.csv
    )

    /*
    * Doublet analysis
    */

    // prep vaeda by removing adata.uns
    prep_doublets(run_topometry.out.adata)

    // run vaeda analysis
    run_vaeda(prep_doublets.out.vaeda, params.batch)

    // run scDblFinder
    run_scDblFinder(prep_doublets.out.scDblFinder, params.batch)

    // combine doublet results
    combine_doublets(
        combine_clusters.out.adata,
        run_vaeda.out.csv,
        run_scDblFinder.out.csv,
        params.batch
    )

    emit:
    /*
    * Cell calling
    */

    // barcode rank plot
    plot_barcode_rank    = params.skip_cellbender ? channel.empty() : plot_barcode_rank.out.plot

    // adatas
    adata_joined_batches = params.skip_cellbender ? channel.empty() : join_adatas.out.all

    // all cellbender outputs
    cellbender_h5       = params.skip_cellbender ? channel.empty() : run_cellbender.out.h5.map      { _batch, path -> path }
    cellbender_report   = params.skip_cellbender ? channel.empty() : run_cellbender.out.report.map  { _batch, path -> path }
    cellbender_log      = params.skip_cellbender ? channel.empty() : run_cellbender.out.log.map     { _batch, path -> path }
    cellbender_pdf      = params.skip_cellbender ? channel.empty() : run_cellbender.out.pdf.map     { _batch, path -> path }
    cellbender_metrics  = params.skip_cellbender ? channel.empty() : run_cellbender.out.metrics.map { _batch, path -> path }

    /*
    * Doublets
    */

    // adatas
    adata_doublets = combine_doublets.out.adata

    // csvs
    csv_vaeda = run_vaeda.out.csv
    csv_scdblfinder = run_scDblFinder.out.csv

    // plots
    plot_doublets_comparison = combine_doublets.out.plot

    // session info
    sessionInfo_scdblfinder = run_scDblFinder.out.session

    /*
    * Dimensionality reduction
    */
    plot_topo_pca = run_topometry.out.plt_pca
    plot_topo_geom = run_topometry.out.plt_geometry
    plot_topo_embeddings = run_topometry.out.plt_comparison

    /*
    * Clustering
    */

    // clustree
    plot_clustree = clustree.out.plot

    // scSHC
    csv_scSHC  = params.skip_scSHC ? channel.empty() : scSHC.out.csv
    plot_scSHC = params.skip_scSHC ? channel.empty() : scSHC.out.plot
    rds_scSHC  = params.skip_scSHC ? channel.empty() : scSHC.out.results_object

    // combine clusters
    plot_cluster_comparison = combine_clusters.out.plot

    /*
    * QC reporting
    */

    // plots
    qcplots_unfiltered = params.skip_cellbender ? channel.empty() : qc_unfiltered.out
    qcplots_cellcalling = qc_cellCalling.out
    qcplots_minFilter = qc_minFilter.out

    // stats
    stats_unfiltered        = params.skip_cellbender ? channel.empty() : report_unfiltered.out.stats
    stats_cellCalling       = report_cellCalling.out.stats
    stats_minFilter         = report_minFilter.out.stats
    batch_stats_unfiltered  = params.skip_cellbender ? channel.empty() : report_unfiltered.out.batch_stats
    batch_stats_cellCalling = report_cellCalling.out.batch_stats
    batch_stats_minFilter   = report_minFilter.out.batch_stats
}
