#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def add_doublet_results(adata, results_csv):
    results = pd.read_csv(results_csv, index_col=0)
    colname = results.columns[0]
    results = pd.Categorical(results[colname], ordered=True, categories=['singlet', 'doublet'])
    adata.obs[colname] = results


def _to_bool(series):
    if series.dtype.name == 'category' or series.dtype == object:
        return series == 'doublet'
    return series.astype(bool)


def find_doublet_union(adata, doublet_columns):
    bool_df = adata.obs[doublet_columns].apply(_to_bool)
    adata.obs['doublet_union'] = (bool_df.sum(axis=1) >= 1).map({True: 'doublet', False: 'singlet'}).astype('category')
    adata.obs.doublet_union.cat.reorder_categories(['singlet', 'doublet'])


def plot_doublets(adata, doublet_columns, basis='projection'):
    sc.pl.embedding(
        adata,
        basis=basis,
        color=doublet_columns,
        palette={'singlet': 'gray', 'doublet': 'red'},
        alpha=0.6,
        frameon=False,
        show=False,
    )
    plt.savefig('doublets_comparison.png')
    plt.close()


def report_doublets(adata, doublet_columns, batch=None):
    for column in doublet_columns:
        is_doublet = _to_bool(adata.obs[column])
        print(f'Total doublets - {column} : {is_doublet.sum():,}')
        if batch:
            print(f'Doublets per {batch} - {column} :')
            for group in adata.obs[batch].unique():
                print(f'  {group} : {is_doublet[adata.obs[batch] == group].sum():,}')


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--vaeda-csv', required=False)
    p.add_argument('--scdblfinder-csv', required=False)
    p.add_argument('--batch', required=True)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)
    arg_csvs = [csv for csv in [args.vaeda_csv, args.scdblfinder_csv] if csv is not None]

    for csv in arg_csvs:
        add_doublet_results(adata, csv)

    doublet_columns = [col for col in adata.obs if col in ['vaeda_calls', 'scDblFinder_calls']]
    find_doublet_union(adata, doublet_columns)
    doublet_columns = doublet_columns + ['doublet_union']

    plot_doublets(adata, doublet_columns=doublet_columns)
    report_doublets(adata, doublet_columns, batch=args.batch)
    adata.write('adata_doublets.h5ad')


if __name__ == '__main__':
    main()
