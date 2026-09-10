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


sc.settings.verbosity = 0
sc.settings.figdir = "."


def plt_proportion(adata, batch, descriptor):
    batches = adata.obs[batch].unique()
    total_cells = adata.shape[0]

    proportions = [adata[adata.obs[batch] == x].shape[0] / total_cells for x in batches]
    prop_dict = dict(zip(batches, proportions))
    prop_df = pd.Series(dict(sorted(prop_dict.items(), key=lambda item: item[1], reverse=False)))

    fig, ax = plt.subplots(figsize=(8, 3))
    p = ax.barh(y=prop_df.index, width=prop_df.values)
    ax.bar_label(p, label_type='edge', fmt='{:.2f}', fontsize=5, padding=5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(which='both', labelsize=5)
    fig.suptitle(f"Proportion of data per {batch}", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    return fig


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--batch', required=True)
    p.add_argument('--descriptor', required=True)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)

    fig = plt_proportion(adata, args.batch, args.descriptor)
    fig.savefig(f'plt_proportion_{args.batch}_{args.descriptor}.png', bbox_inches='tight', transparent=True, dpi=300)
    plt.close(fig)

    fig = plt_proportion(adata, "predicted_celltype", args.descriptor)
    fig.savefig(f'plt_proportion_predicted_celltype_{args.descriptor}.png', bbox_inches='tight', transparent=True, dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
