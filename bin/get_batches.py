#!/usr/bin/env python3
import argparse
import anndata as ad


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--batch-col', required=True)
    return p


def main():
    args = build_parser().parse_args()
    adata = ad.read_h5ad(args.adata)
    if adata.raw:
        adata = adata.raw.to_adata()
    for b in adata.obs[args.batch_col].unique():
        print(b)


if __name__ == '__main__':
    main()
