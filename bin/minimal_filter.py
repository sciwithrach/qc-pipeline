#!/usr/bin/env python3
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--descriptor', required=True)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.write(f'adata_minimal_filter_{args.descriptor}.h5ad')


if __name__ == '__main__':
    main()
