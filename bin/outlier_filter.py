#!/usr/bin/env python3
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sb
from scipy.stats import median_abs_deviation
import warnings
import re
import logging

###################
# OUTLIER BOOL ####
###################

def is_outlier(df, metadata, nmads, which):
    M = df[metadata]
    upper = np.median(M) + nmads * median_abs_deviation(M) < M
    lower = M < np.median(M) - nmads * median_abs_deviation(M)
    if which == 'upper':
        return upper
    elif which == 'lower':
        return lower
    elif which == 'both':
        outlier = upper | lower
        return outlier

###################
#### BOOL PLOT ####
###################

def plt_clusterwise_bool(adata, clusters, outlier_cols):
    n_outliers = len(outlier_cols)
    fig, axs = plt.subplots(n_outliers, figsize=(12, 3.5 * n_outliers))

    for i, outlier in enumerate(outlier_cols):
        sb.histplot(
            data=adata.obs,
            x=clusters,
            hue=outlier,
            multiple='fill',
            ax=axs[i],
            palette={False: 'lightgray', True: 'tab:red'},
            edgecolor='white',
            legend=False
        )
        label = outlier.replace("n_", "").replace("outlier_", "").replace("mread", "read").replace("pct_counts_mito", "mito%").replace("_union", "s").replace("_", " ")
        axs[i].set_title(label, fontsize=15, loc='left')
        axs[i].set_xlabel('')
        axs[i].set_ylabel('Proportion', fontsize=8)
        axs[i].bar_label(axs[i].containers[0], fmt='{:.2f}', fontsize=8)

        if re.search('count$', outlier) and len(axs[i].containers) > 1:
            with warnings.catch_warnings():
                warnings.simplefilter(action='ignore', category=FutureWarning)
                thresholds = [f'{x:,.0f}' for x in adata[adata.obs[outlier]].obs.groupby(clusters)[outlier.replace("outlier_", "")].min()]
            axs[i].bar_label(
                axs[i].containers[1],
                #label_type='center',
                labels = thresholds,
                fontsize=8, padding = -40, rotation=90
                #weight = 'bold'
                )

    axs[i].spines[['top', 'right']].set_visible(False)
    plt.suptitle(f'Proportion of outliers by cluster ({clusters})', fontsize=20)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    return fig

###################
### PLOT BOXEN ####
###################

def plt_clusterwise_boxen(adata, clusters, outlier_cols):
    mito_col = 'outlier_pct_counts_mito'
    n_extra = 1 if mito_col in outlier_cols else 0
    n_rows = len(outlier_cols) + n_extra

    fig, axs = plt.subplots(n_rows, figsize=(15, 3 * n_rows))

    ax_i = 0
    for outlier in outlier_cols:
        metric = outlier.replace('outlier_', '')
        sb.boxenplot(
            data=adata.obs,
            x=clusters,
            y=metric,
            hue=outlier,
            palette={False: 'gray', True: 'tab:red'},
            legend=False,
            ax=axs[ax_i]
        )
        metric = metric.replace("pct_counts_mito", "mito%").replace("_", " ").replace("tscp", "transcript").replace("mread", "read")
        axs[ax_i].set_title(metric, loc='left', fontsize=15)
        axs[ax_i].set_xlabel('')
        axs[ax_i].set_ylabel('')
        if metric in ('transcript count', 'read count'):
            axs[ax_i].set_yscale('log', base=np.exp(1))
            axs[ax_i].set_yticks([np.exp(x) for x in range(2,18,2)], [f'{np.exp(x):,.1e}' for x in range(2,18,2)])
            axs[ax_i].set_ylim(ymin=10)
            axs[ax_i].set_title(f'log1p {metric}', loc='left', fontsize=15)
        axs[ax_i].spines[['top', 'right']].set_visible(False)
        
        ax_i += 1
    
    if mito_col in outlier_cols:
        # plot zoomed in mito% up to 1%
        sb.boxenplot(
            data=adata.obs,
            x=clusters,
            y="pct_counts_mito",
            hue=mito_col,
            palette=['gray', 'tab:red'],
            legend=False,
            ax=axs[ax_i]
        )
        axs[ax_i].set_title("mito% (up to 1%)", loc='left', fontsize=15)
        axs[ax_i].set_xlabel('')
        axs[ax_i].set_ylabel('')
        axs[ax_i].set_ylim(0, 1)
        axs[ax_i].spines[['top', 'right']].set_visible(False)

        fig.suptitle(f'Per cell outlier values by cluster: {clusters}', fontsize=20)
        fig.tight_layout(rect=[0, 0, 1, 0.98])

    return fig

########################
# PLOT LIBRARY SIZE ####
########################

def plt_library_size(adata, metric_list):
    
    # setup
    n_rows = len(metric_list)
    fig, axs = plt.subplots(n_rows, figsize=(12, 6 * n_rows))

    # plot
    for i, metric in enumerate(metric_list):
        outlier_metric = f'outlier_{metric}'
        sb.scatterplot(
            data = adata.obs, 
            x='log1p_total_counts', 
            y=metric, 
            hue=outlier_metric, 
            palette={'gray', 'tab:red'},
            legend=False,
            ax=axs[i]
            )
        axs[i].spines[['top', 'right']].set_visible(False)
        
        # set log axes where applicable
        metric = metric.replace("pct_counts_mito", "mito%").replace("_", " ").replace("tscp", "transcript").replace("mread", "read")
        if metric in ('transcript count', 'read count'):
            axs[i].set_yscale('log', base=np.exp(1))
            axs[i].set_yticks([np.exp(x) for x in range(2,18,2)], [f'{np.exp(x):,.1e}' for x in range(2,18,2)])
            axs[i].set_ylim(ymin=10)
            axs[i].set_title(f'log1p {metric}', loc='left', fontsize=15)
        else:
            axs[i].set_title(metric, loc='left', fontsize=15)
    
    # layout
    fig.suptitle('Outliers by library size', fontsize=20)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    return fig

###################
# CALCULATE QC ####
###################

def calculate_qc(adata):
    adata.var['mito'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.str.startswith(('Rrp', 'Rp'))
    adata.var['hemo'] = adata.var_names.str.startswith('Hb')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito', 'ribo', 'hemo'], percent_top=[20], inplace=True)

###################
# GET OUTLIERS ####
###################

def get_outliers(adata, filter_list, clusters, doublet_column, descriptor):

    # calculate qc metrics if needed
    if 'pct_counts_mito' not in adata.obs.columns:
        calculate_qc(adata)

    # initial report
    logging.warning('OUTLIER REPORT')
    logging.warning(f'\nTotal starting reads: {adata.obs.mread_count.sum():,.2e}')
    logging.warning(f'Total starting cells: {adata.obs.shape[0]:,.0f}')
    logging.warning(f'Cell count by cluster: {[f'C{i}: {x:,.0f}' for i, x in enumerate(adata.obs.groupby(clusters).size())]}')
    logging.warning('\n')

    # doublets are excluded first, ahead of MAD calculation, so per-cluster
    # medians/MADs used below reflect the post-doublet-removal population
    adata.obs['outlier_doublet'] = (adata.obs[doublet_column] == 'doublet')
    logging.warning(f'total doublets: {adata.obs.outlier_doublet.sum():,.0f}')

    # reference population for MAD-based outlier detection below: doublets
    # already removed, so thresholds aren't skewed by them
    survivors = adata.obs.loc[~adata.obs.outlier_doublet]

    # determine outliers
    
    for metric in filter_list:
        # 3 mads lower filter for all metrics
        which = 'lower' 
        # apply both 3 mads upper AND lower filter for gene count
        if metric =='gene_count': 
            which = 'both'
        # initialise obs column
        adata.obs[f'outlier_{metric}'] = False
        # identify outliers in singlets only
        # flag cells with >75 000 reads and/or 3 mads above cluster median (note: 75k threshold applied independently to median calculation)
        if metric == 'mread_count':
            flagged = survivors.groupby(clusters, group_keys=False).apply(is_outlier, metric, 3, which) | (survivors.mread_count > 75000)
        # flag outliers for all other metrics
        else:
            flagged = survivors.groupby(clusters, group_keys=False).apply(is_outlier, metric, 3, which)
        # mark flagged outliers in obs column
        adata.obs.loc[flagged.index, f'outlier_{metric}'] = flagged
        # report total outliers
        logging.warning(f'total {metric} outliers: {adata.obs[f"outlier_{metric}"].sum():,.0f}')


    # determine mito% outliers
    
    # initialise obs column
    adata.obs['outlier_pct_counts_mito'] = False
    # identify outliers in singlets only - 5 mads upper bound and >1%, or >5% regardless of mads
    mito_flagged = survivors.groupby(clusters, group_keys=False).apply(is_outlier, 'pct_counts_mito', 5, 'upper') & (survivors.pct_counts_mito > 1) | (survivors.pct_counts_mito > 5)
    # mark flagged outliers in obs column
    adata.obs.loc[mito_flagged.index, 'outlier_pct_counts_mito'] = mito_flagged
    # report total outliers
    logging.warning(f'total mito outliers: {adata.obs["outlier_pct_counts_mito"].sum():,.0f}')
    # select columns denoting outliers from adata.obs
    outlier_cols = [x for x in adata.obs.columns if x.startswith('outlier')]

    # create qc fail column for boolean plots

    # identify cells failing any metric
    adata.obs['outlier_qc_fail'] = adata.obs[outlier_cols].astype(bool).any(axis=1)
    # log total outliers across all metrics
    logging.warning(f'total qc fail: {adata.obs.outlier_qc_fail.sum():,.0f}')
    # report outlier columns
    logging.warning(f'\noutlier columns: {outlier_cols}')

    # intialise outlier df
    outlier_df = pd.DataFrame(index=adata.obs[clusters].cat.categories)

    # report outlier metrics, where qc_fail gives the total stat
    for outlier_metric in outlier_cols + ['outlier_qc_fail']:
        logging.warning(f'\n{outlier_metric}')
        
        # column name without outlier prefix
        metric = outlier_metric.replace('outlier_', '')
        
        # total reads in outliers per cluster
        outlier_df[f'n_reads_removed_{outlier_metric}'] = adata[adata.obs[outlier_metric]].obs.groupby(clusters).mread_count.sum().astype(int)
        logging.warning(f'total {metric} reads removed: {outlier_df[f'n_reads_removed_{outlier_metric}'].sum():,.2e}')
       
        # number of outliers per cluster
        outlier_df[f'n_outliers_removed_{outlier_metric}'] = adata.obs.groupby(clusters)[outlier_metric].sum().astype(int)
        logging.warning(f'total {metric} outliers removed: {outlier_df[f'n_outliers_removed_{outlier_metric}'].sum():,.0f}')
        
        # for quantifiable outliers (tscps, reads, genes - could change this to assert numerical type instead)
        if outlier_metric.replace('outlier_', '') in filter_list:
            
            # median value amongst outliers
            outlier_df[f'median_of_{metric}_outliers'] = adata[adata.obs[outlier_metric]].obs.groupby(clusters)[metric].median()
            logging.warning(f'{metric} threshold by cluster:\n {[f"C{i}: {x:,.0f}" for i,x in enumerate(adata[adata.obs[outlier_metric]].obs.groupby(clusters)[metric].min())]}')
            logging.warning(f'median {metric} in REMOVED cells by cluster:\n {[f"cluster {i}: {x:,.0f}, " for i,x in enumerate(outlier_df[f'median_of_{metric}_outliers'])]}')
            logging.warning(f'median {metric} in RETAINED cells by cluster:\n {[f"cluster {i}: {x:,.0f}, " for i,x in enumerate(adata[~adata.obs[outlier_metric]].obs.groupby(clusters)[metric].median())]}\n')

    # save outlier dataframe
    outlier_df.to_csv(f'outliers_{descriptor}.csv')

    # return outlier columns for downstream analysis
    return outlier_cols

###################
###### RUN CLI ####
###################

def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--clusters', required=True)
    p.add_argument('--doublets', required=True)
    p.add_argument('--descriptor', required=True)
    return p


def main():
    # load arguments
    args = build_parser().parse_args()

    # logging setup
    logging.basicConfig(
        filename = f"outliers_{args.descriptor}.log",
        encoding="utf-8",
        filemode="a",
        level=logging.WARNING,
        format='{asctime} {levelname} {message}',
        style='{'
    )

    # read adata
    adata = ad.read_h5ad(args.adata)

    # convert to raw if available
    if adata.raw:
        adata = adata.raw.to_adata()

    # metadata to filter by (except mito%)
    filter_list = ['tscp_count', 'mread_count', 'gene_count']

    # setup outliers
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)
        outlier_cols = get_outliers(adata, filter_list, args.clusters, args.doublets, args.descriptor)   

    # plot boolean summary
    fig = plt_clusterwise_bool(adata, args.clusters, outlier_cols)
    fig.savefig(f'plt_mads_boolean_{args.descriptor}.png', transparent=True, dpi=300)
    plt.close(fig)

    # plot boxen summary (continuous metrics only, not doublet/small-cluster exclusions)
    boxen_cols = [f'outlier_{m}' for m in filter_list] + ['outlier_pct_counts_mito']
    fig = plt_clusterwise_boxen(adata, args.clusters, boxen_cols)
    fig.savefig(f'plt_mads_boxen_{args.descriptor}.png', transparent=True, dpi=300)
    plt.close(fig)

    # plot against library size 
    fig = plt_library_size(adata, filter_list + ['pct_counts_mito'])
    fig.savefig(f'plt_mads_outliers_by_counts_{args.descriptor}.png', transparent=True, dpi=300)
    plt.close(fig)

    # filter adata
    adata = adata[~adata.obs.outlier_qc_fail].copy()
    
    # report after filtering
    logging.warning(f'\nRemaining cells: {adata.shape[0]:,.0f}')
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)
        logging.warning(f'Remaining cell count by cluster: {[f'C{i}: {x:,.0f}' for i, x in enumerate(adata.obs.groupby(args.clusters).size())]}')
    
    logging.warning(f'\nRemaining reads: {adata.obs.mread_count.sum():,.2e}\n')
    if 'plate' in adata.obs.columns:
        for plate in adata.obs.plate.cat.categories:
            logging.warning(f'Remaining reads for plate {plate}: {adata[adata.obs.plate == plate].obs.mread_count.sum():,.2e}')
    
    # save filtered adata
    adata.write_h5ad(f'adata_final_filter_{args.descriptor}.h5ad')


if __name__ == '__main__':
    main()
