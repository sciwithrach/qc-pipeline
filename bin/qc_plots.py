#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sb


sc.settings.verbosity = 0
sc.set_figure_params(dpi_save=300, transparent=True, frameon=False, format="png")


def calculate_qc(adata):
    adata.var['mito'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.str.startswith(('Rrp', 'Rp'))
    adata.var['hemo'] = adata.var_names.str.startswith('Hb')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito', 'ribo', 'hemo'], percent_top=[20], inplace=True)


def pltLogCounts(adata, batch_id, suffix):
    g = sb.pairplot(
        adata.obs,
        vars=["tscp_count", "mread_count", 'gene_count'],
        hue=batch_id,
        corner=True,
        diag_kws=dict(log_scale=True, common_norm=False)
    )
    g.fig.subplots_adjust(top=0.93)
    g.fig.suptitle(f'Pairwise relationships between number of transcripts, reads, and genes ({adata.shape[0]} barcodes)')
    plt.tight_layout()
    g.savefig(f"qcplot_logCounts_{suffix}.png")
    plt.close()


def pltCounts(adata, batch_id, suffix):
    g = sb.pairplot(
        adata.obs,
        vars=["tscp_count", "mread_count", 'gene_count'],
        hue=batch_id,
        corner=True,
        diag_kws=dict(common_norm=False)
    )
    g.fig.subplots_adjust(top=0.93)
    g.fig.suptitle(f'Pairwise relationships between number of transcripts, reads, and genes ({adata.shape[0]} barcodes)')
    plt.tight_layout()
    g.savefig(f"qcplot_rawCounts_{suffix}.png")
    plt.close()


def pltGenesLibrarySize(adata, suffix):
    fig = plt.figure(figsize=(8*3, 6*1))
    ax = fig.add_subplot(1, 3, 1)
    ax.hist(adata.obs['log1p_total_counts'], bins=100)              # Total number of counts for a cell
    ax.set_xlabel('Log1p library size', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.set_title('Histogram of log1p library size', fontsize=14)

    ax = fig.add_subplot(1, 3, 2)
    ax.hist(adata.obs['log1p_n_genes_by_counts'], bins=100)         # The number of genes with at least 1 count in a cell
    ax.set_xlabel('Log of num. genes per cell', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.set_title('Histogram of number of genes per barcode', fontsize=14)

    ax = fig.add_subplot(1, 3, 3)
    x = adata.obs['log1p_total_counts']
    y = adata.obs['log1p_n_genes_by_counts']
    ax.scatter(x, y, s=5)
    ax.set_ylabel('Log of num. genes per cell', fontsize=14)
    ax.set_xlabel('Log library size', fontsize=14)
    corr_coef = np.corrcoef(x, y)[0, 1]
    ax.set_title('Correlation = ' + str(round(corr_coef, 3)), fontsize=14)

    plt.tight_layout()
    fig.savefig(f'qcplot_geneLibrarySize_{suffix}.png')
    plt.close()


def pltBarcodesPerGene(adata, suffix):
    fig = plt.figure(figsize=(8*3, 6*1))
    ax = fig.add_subplot(1, 3, 1)
    ax.hist(adata.var['n_cells_by_counts'], bins=100)               # Number of cells this expression is measured in
    ax.set_xlabel('Number of cells a gene is expressed in', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.set_title('Histogram of number of cells each gene is expressed in', fontsize=14)

    ax = fig.add_subplot(1, 3, 2)
    ax.hist(np.log(adata.var['n_cells_by_counts'] + 1), bins=100)
    ax.set_xlabel('Log - Number of cells a gene is expressed in', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.set_title('Histogram of log of number of cells each gene is expressed in', fontsize=14)

    ax = fig.add_subplot(1, 3, 3)
    ax.hist(np.log(adata.var['n_cells_by_counts'] + 1), bins=100)
    ax.set_xlabel('Log - Number of cells a gene is expressed in', fontsize=14)
    ax.set_ylabel('Frequency, ylim trimmed', fontsize=14)
    ax.set_title('Histogram of log of number of cells each gene is expressed in', fontsize=14)
    ax.set_ylim([0, 1000])

    plt.tight_layout()
    fig.savefig(f'qcplot_barcodesPerGene_{suffix}.png')
    plt.close()


def pltViolinQC(adata, suffix):
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15, 4), gridspec_kw={"wspace": 0.5})
    sc.pl.violin(adata, ['pct_counts_hemo'], jitter=0.4, show=False, ax=ax1)
    sc.pl.violin(adata, ['pct_counts_mito'], jitter=0.4, show=False, ax=ax2)
    sc.pl.violin(adata, ['pct_counts_ribo'], jitter=0.4, show=False, ax=ax3)
    sc.pl.violin(adata, ["log1p_n_genes_by_counts"], jitter=0.4, show=False, ax=ax4)
    fig.savefig(f'qcplot_metadataViolins_{suffix}.png')
    plt.close()


def pltHistQC(adata, suffix):
    metadata_list = ['mread_count', 'tscp_count', 'gene_count']
    n_metadata = len(metadata_list)
    figure_ids = zip(*(iter(np.arange(n_metadata*3) + 1),) * 3)
    fig_dict = dict(zip(np.arange(n_metadata), figure_ids))
    fig = plt.figure(figsize=(8*3, 6*n_metadata))
    for i, metadata in enumerate(metadata_list):
        ids = fig_dict[i]
        ax = fig.add_subplot(n_metadata, 3, ids[0])
        ax.hist(adata.obs[metadata], 100)
        ax.set_xlabel(metadata, fontsize=14)
        ax.set_ylabel(metadata)
        ax = fig.add_subplot(n_metadata, 3, ids[1])
        ax.scatter(adata.obs['log1p_total_counts'], adata.obs[metadata])
        ax.set_xlabel('Log library size', fontsize=14)
        ax.set_ylabel(metadata, fontsize=14)
        ax = fig.add_subplot(n_metadata, 3, ids[2])
        ax.scatter(adata.obs['log1p_n_genes_by_counts'], adata.obs[metadata])
        ax.set_xlabel('Log num. genes per cell', fontsize=14)
        ax.set_ylabel(metadata, fontsize=14)
    plt.tight_layout()
    fig.savefig(f'qcplot_metadataHistograms_{suffix}.png')
    plt.close()


def pltLogReads(adata, suffix):
    metadata_list = ['mread_count', 'tscp_count', 'gene_count']
    n_metadata = len(metadata_list)
    fig, axs = plt.subplots(1, n_metadata, figsize=(8*n_metadata, 6))
    for i, metadata in enumerate(metadata_list):
        axs[i].scatter([np.log1p(x) for x in adata.obs['mread_count']], adata.obs[metadata])
        axs[i].set_xlabel('Log read count', fontsize=14)
        axs[i].set_ylabel(metadata, fontsize=14)
    fig.savefig(f'qcplot_metadataLogReads_{suffix}.png')
    plt.close()


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--batch', required=True)
    p.add_argument('--descriptor', required=True)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)

    if 'log1p_total_counts' not in adata.obs.columns:
        if adata.raw:
            adata.layers['counts'] = adata.raw.to_adata()
        calculate_qc(adata)

    pltLogCounts(adata, suffix=args.descriptor, batch_id=args.batch)
    pltCounts(adata, suffix=args.descriptor, batch_id=args.batch)
    pltGenesLibrarySize(adata, suffix=args.descriptor)
    pltBarcodesPerGene(adata, suffix=args.descriptor)
    pltViolinQC(adata, suffix=args.descriptor)
    pltHistQC(adata, suffix=args.descriptor)
    pltLogReads(adata, suffix=args.descriptor)


if __name__ == '__main__':
    main()
