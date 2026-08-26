# Collapse Ensembl IDs to unique gene symbols

library(tidyverse)
library(devtools)
load_all("source_files/Breast-Cancer-Modules/")

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/final/"
dir.create(paste0(out_dir), recursive = T, showWarnings = F)

## Load in expression matrix with Ensembl IDs
ensembl_expr <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_primaries_combined_ensembl_id.tsv",
                      sep = "\t", header = T, row.names = 1, check.names = F)

## Filter for Ensembl IDs in gene map found in expr data
gene_mapping <- gene_map_human %>% 
  filter(ensembl_id %in% row.names(ensembl_expr))

## Identify gene IDs that are duplicated in mapping
duplicate_symbols <- unique(gene_mapping$gene_name[duplicated(gene_mapping$gene_name)])
duplicate_ensembl <- unique(gene_mapping$ensembl_id[gene_mapping$gene_name %in% duplicate_symbols])

## Collapse duplicated genes to unique Gene Symbols
symbol_expr <- matrix(NA, nrow = length(duplicate_symbols), 
                      ncol = ncol(ensembl_expr), 
                      dimnames = list(duplicate_symbols, 
                                      colnames(ensembl_expr)))

for(symbol in duplicate_symbols) {
  ids <- gene_mapping$ensembl_id[gene_mapping$gene_name == symbol] 
  symbol_expr[symbol,] <- colMeans(ensembl_expr[ids,])
}

gene_mapping <- gene_mapping %>% 
  filter(!(ensembl_id %in% duplicate_ensembl))

sub_expr <- ensembl_expr[gene_mapping$ensembl_id,]
row.names(sub_expr) <- gene_mapping$gene_name

## Combine final Gene Symbol expr data
symbol_expr <- rbind(sub_expr, symbol_expr)

write.table(cbind(Symbol = row.names(symbol_expr), symbol_expr), 
            paste0(out_dir, "UNC337_n239_primaries_combined_gene_symbols.tsv"), 
            sep = "\t", row.names = F, quote = F)
