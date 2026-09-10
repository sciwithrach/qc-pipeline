#!/usr/bin/env Rscript

set.seed(42)

library(scSHC)
library(anndataR)
library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)

args       <- commandArgs(trailingOnly = TRUE)
adata_path <- args[1]
nr_cores   <- as.integer(args[2])

adata     <- read_h5ad(adata_path)
sce       <- adata$as_SingleCellExperiment()
count_mat <- assay(sce, assayNames(sce)[1])

results <- scSHC(count_mat, cores = nr_cores)

results_csv <- data.frame(
    shc_clusters = results[[1]],
    row.names    = colnames(sce)
)

write.csv(results_csv, "scSHC_clusters.csv")
saveRDS(results, "scSHC_results.rds")

p <- ggplot(results_csv, aes(x = shc_clusters)) +
    geom_bar() +
    theme_minimal() +
    labs(x = "Cluster", y = "Cell count", title = "scSHC cluster sizes")
ggsave("scSHC_clusters.png", p, width = 8, height = 6)
