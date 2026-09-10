#!/usr/bin/env python3
import argparse
import anndata as ad
import vaeda
import pandas as pd


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--batch', required=True)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)

    adatas = []
    for batch in adata.obs[args.batch].cat.categories:
        subset = adata[adata.obs[args.batch] == batch].copy()
        subset = vaeda.vaeda(subset, filter_genes=False, seed=42)
        adatas.append(subset)
    adata = ad.concat(adatas)

    results = adata.obs[['vaeda_calls', 'vaeda_scores']]
    results.to_csv('doublets_vaeda.csv')


if __name__ == '__main__':
    main()
