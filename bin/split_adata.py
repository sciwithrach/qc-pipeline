#!/usr/bin/env python3
import argparse
import anndata as ad


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--batch', required=True)
    p.add_argument('--batch-col', required=True)
    return p


def main():
    args = build_parser().parse_args()
    adata = ad.read_h5ad(args.adata)
    if adata.raw:
        adata = adata.raw.to_adata().copy()
    bdata = adata[adata.obs[args.batch_col] == args.batch].copy()
    bdata.write_h5ad(f'adata_split_{args.batch}.h5ad')


if __name__ == '__main__':
    main()
