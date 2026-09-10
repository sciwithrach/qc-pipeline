#!/usr/bin/env python3
import argparse
import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--celltype-csv', required=True)
    p.add_argument('--scshc-csv', required=True)
    p.add_argument('--descriptor', required=False)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)
    ctype = pd.read_csv(args.celltype_csv, index_col=0, header=0)

    # rename old clusters & cell typing if they exist
    tag = 'qc' if args.descriptor == 'postqc' else f'before_{args.descriptor}'
    for key in ['shc_clusters', 'predicted_celltype']:
        if key in adata.obs:
            adata.obs[f'{key}_{tag}'] = adata.obs[key]
            del adata.obs[key]

    if os.path.exists(args.scshc_csv):
        shc = pd.read_csv(args.scshc_csv, index_col=0, header=0)
        adata.obs = adata.obs.join([ctype.predicted_celltype, shc.shc_clusters])
        adata.obs['shc_clusters'] = adata.obs['shc_clusters'].astype(int).astype('category')
    else:
        adata.obs = adata.obs.join(ctype.predicted_celltype)
    adata.obs['predicted_celltype'] = adata.obs['predicted_celltype'].astype('category')

    if 'ms_spectral' in adata.uns['basis']:
        cluster_cols = [col for col in adata.obs.columns if col.startswith('topo_clusters_ms_')]
    elif 'spectral' in adata.uns['basis']:
        cluster_cols = [col for col in adata.obs.columns if col.startswith('topo_clusters_res_')]
    else:
        cluster_cols = [col for col in adata.obs.columns if col.startswith('pca_leiden_res')]

    sc.pl.embedding(
        adata,
        basis='projection',
        color=['predicted_celltype'] + (['shc_clusters'] if 'shc_clusters' in adata.obs.columns else []) + cluster_cols,
        legend_loc='on data',
        legend_fontsize='xx-small',
        frameon=False,
        show=False,
    )
    plt.savefig('cluster_comparison_on_embedding.png')
    plt.close()

    adata.write('adata_clusters.h5ad')


if __name__ == '__main__':
    main()
