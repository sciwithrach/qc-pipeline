#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd
import cellbender.remove_background.downstream as cb


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--raw-adata', required=True)
    p.add_argument('cellbender_paths', nargs='+')
    return p


def main():
    args = build_parser().parse_args()

    cellbender_paths = sorted(args.cellbender_paths)

    cell_probs_list = []
    for cellbender in cellbender_paths:
        adata_cb = cb.anndata_from_h5(cellbender)
        column = adata_cb.obs.cell_probability
        cell_probs_list.append(column)
    cell_probs_col = pd.concat(cell_probs_list)

    adata = ad.read_h5ad(args.raw_adata)
    adata.obs = adata.obs.join(cell_probs_col)
    adata.write_h5ad('adata_cellbender_unfiltered.h5ad')

    adata = adata[adata.obs.cell_probability > 0.5].copy()
    adata.write_h5ad('adata_cellbender_filtered.h5ad')


if __name__ == '__main__':
    main()
