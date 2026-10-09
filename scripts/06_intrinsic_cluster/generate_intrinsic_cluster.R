# Generate intrinsic cluster + cluster 3.0 files for PAM50 subtyping

library(tidyverse)
library(matrixStats)

dwd <- read.table("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols.tsv", 
                   sep = "\t", header = T, row.names = 1, check.names = F)

## Filter expression data to include only intrinsic genes 
intrinsic_list <- unlist(read.table("source_files/Human intrinisic wMouseID.txt", 
                             sep = "\t", header = T)[[3]])

intrinsic_expr <- as.matrix(dwd[row.names(dwd) %in% intrinsic_list,])
intrinsic_expr <- as.data.frame(intrinsic_expr - rowMedians(intrinsic_expr))
intrinsic_expr <- tibble::rownames_to_column(intrinsic_expr, var = "Gene")

write.table(intrinsic_expr,
          "output/intrinsic_cluster/n2486_Ensembl_dwd_symbols_intrinsic_genes.txt", 
          sep = "\t", row.names = F, col.names = T)

## Load expression data externally into Cluster 3.0 for hierarchical clustering:
# Cluster genes and arrays (using the "Correlation (centered)" Similarity metric),
# with the centroid linkage clustering method.

