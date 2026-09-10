#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd


def compute_stats(obs, name):
    metadata_list = ['tscp_count', 'gene_count', 'mread_count']
    df = pd.DataFrame(index=[name])
    for metadata in metadata_list:
        df[f'total_{metadata}'] = obs[metadata].sum()
        df[f'median_{metadata}'] = obs[metadata].median()
        df[f'mean_{metadata}'] = obs[metadata].mean()
        df[f'min_{metadata}'] = obs[metadata].min()
        df[f'max_{metadata}'] = obs[metadata].max()
    df['total_nuclei'] = len(obs)
    return df


def report_stats(adata, colname, batch_col=None):
    compute_stats(adata.obs, colname).T.to_csv(f'report_stats_{colname}.csv')
    if batch_col is not None:
        rows = [compute_stats(adata[adata.obs[batch_col] == val].obs, val)
                for val in sorted(adata.obs[batch_col].unique())]
        pd.concat(rows).T.to_csv(f'report_stats_{colname}_batch.csv')


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--descriptor', required=True)
    p.add_argument('--batch', default=None)
    return p


def main():
    args = build_parser().parse_args()
    adata = ad.read_h5ad(args.adata)
    report_stats(adata, args.descriptor, batch_col=args.batch)


if __name__ == '__main__':
    main()
