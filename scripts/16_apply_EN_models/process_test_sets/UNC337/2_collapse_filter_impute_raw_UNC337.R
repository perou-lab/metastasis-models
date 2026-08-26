# Impute missing values, collapse each UNC337 dataset to remove duplicated ENSEMBL gene IDs

library(tidyverse)
library(impute)
library(annotate)
library("hgu133a.db")
source("source_files/functions.R")
set.seed(123)

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/"
dir.create(paste0(out_dir), recursive = T, showWarnings = F)

gpl_files <- c("GPL7504", "GPL6607", "GPL5325", "GPL1708", "GPL1390", "GPL887", "GPL885")

for (file_name in gpl_files) {
  # Collapse on duplicated ENSEMBL gene IDs (one row / gene ID)
  file <- paste0(out_dir, "/", file_name, "_mapped_ensembl_id.tsv")
  x <- readarray(file, hr=1, impute=F)
  write.table(signif(x$xd,4), paste0(out_dir, "/", file_name, "_mapped_ensembl_id_clpsd.txt"), 
            sep = "\t", col.names = NA)
  
  expr <- read.table(paste0(out_dir, "/", file_name, "_mapped_ensembl_id_clpsd.txt"), 
                      sep = "\t", header = T, row.names = 1)
  
  # Filter probes that are missing 30% of data: 
  expr <- expr[rowSums(!is.na(expr)) >= 0.7*ncol(expr),]
  
  # Impute missing data (k=10, default)
  expr <- impute.knn(as.matrix(expr))
  expr <- expr$data
  
  # Save final processed matrix with ENSEMBL gene IDs:
  write.table(expr, paste0(out_dir, "/", file_name, "_mapped_ensembl_id_clpsd_filt_imp.txt"), 
              sep = "\t", row.names = T, col.names = T)
}
