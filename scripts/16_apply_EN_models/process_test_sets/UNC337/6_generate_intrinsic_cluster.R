# Generate intrinsic cluster for UNC337 (n239 recombined)

library(tidyverse)
library(matrixStats)

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/intrinsic_cluster/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

expr <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_primaries_combined_gene_symbols.tsv", 
                  sep = "\t", header = T, row.names = 1, check.names = F)

## Filter expression data to include only intrinsic genes 
intrinsic_list <- unlist(read.table("source_files/Human intrinisic wMouseID.txt", 
                                    sep = "\t", header = T)[[3]])

intrinsic_expr <- as.matrix(expr[row.names(expr) %in% intrinsic_list,])
intrinsic_expr <- as.data.frame(intrinsic_expr - rowMedians(intrinsic_expr))
intrinsic_expr <- tibble::rownames_to_column(intrinsic_expr, var = "Gene")

write.table(intrinsic_expr,
            paste0(out_dir, "/UNC337_n239_primaries_symbols_intrinsic_genes.txt"), 
            sep = "\t", row.names = F, col.names = T)

## Load expression data externally into Cluster 3.0 for hierarchical clustering:
# Cluster genes and arrays (correlation, centered), using centroid linkage as 
# the clustering method.

