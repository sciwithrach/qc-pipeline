#!/usr/bin/env python

import anndata as ad
import scanpy as sc
import argparse

# parse arguments
p = argparse.ArgumentParser()
p.add_argument(
    "-a", "--adata", 
    help = "adata with raw counts in X", 
    type = str
    )
p.add_argument(
    "-b", "--batch", 
    default = None, 
    help = "batch variable", 
    type = str)
args = p.parse_args()

# load data
raw = ad.read_h5ad(args.adata)

# set universal conditions
# high umi cells
cond1 = raw.obs.tscp_count.between(500,1000)
# low umi cells
cond2 = raw.obs.tscp_count.between(10,50)

# get high and low umi count barcodes

# per batch
if args.batch:
    # initialise list
    adata_list = []
    for batch in raw.obs[args.batch].unique():
        # cells with the correct plate
        cond3 = raw.obs[args.batch]==batch
        new = ad.concat([raw[(cond1 & cond3)][:500,], raw[(cond2 & cond3)][:50000,]])
        adata_list.append(new)
    # join batches
    test = ad.concat(adata_list)
# for the whole dataset
else:
    test = ad.concat([raw[cond1][:500,], raw[cond2][:50000,]])

# save test (before cell calling)
test.write_h5ad('test/adata_test.h5ad')

# save test (after cell calling)
tiny_test = test[test.obs.tscp_count > 500].copy()
tiny_test.write_h5ad('test/adata_test_tiny.h5ad')