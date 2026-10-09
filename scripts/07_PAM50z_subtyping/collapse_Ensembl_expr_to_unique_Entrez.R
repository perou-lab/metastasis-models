# Collapse Ensembl IDs to unique entrez IDs

library(tidyverse)
load("source_files/Breast-Cancer-Modules/data/gene_map_human.rda")
set.seed(123)

out_dir <- "output/PAM50z_subtyping/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Load in expression matrix with Ensembl IDs
ensembl_expr <- read.table("output/DWD_adjust/n2486_Ensembl_dwd.txt",
                      sep = "\t", header = T, row.names = 1, check.names = F)
ensembl_expr <- ensembl_expr[row.names(ensembl_expr) != "ENSG00000285188",] # doesn't map to Entrez ID

## Filter for Ensembl IDs in gene map found in expr data
gene_mapping <- gene_map_human %>% 
  filter(ensembl_id %in% row.names(ensembl_expr))

## Identify gene IDs that are duplicated in mapping
duplicate_entrez <- unique(gene_mapping$entrez_id[duplicated(gene_mapping$entrez_id)])
duplicate_ensembl <- unique(gene_mapping$ensembl_id[gene_mapping$entrez_id %in% duplicate_entrez])

## Collapse duplicated genes to unique Entrez IDs
entrez_expr <- matrix(NA, nrow = length(duplicate_entrez), 
                      ncol = ncol(ensembl_expr), 
                      dimnames = list(duplicate_entrez, 
                                      colnames(ensembl_expr)))

rownames(entrez_expr) <- as.character(rownames(entrez_expr))
duplicate_entrez <- as.character(duplicate_entrez)

for(entrez in duplicate_entrez) {
  ids <- gene_mapping$ensembl_id[!is.na(gene_mapping$entrez_id) & gene_mapping$entrez_id == as.integer(entrez)]
  ids <- ids[ids %in% rownames(ensembl_expr)]
  entrez_expr[entrez, ] <- colMeans(ensembl_expr[ids, , drop=FALSE])
}

gene_mapping <- gene_mapping %>% 
  filter(!(ensembl_id %in% duplicate_ensembl))

sub_expr <- ensembl_expr[gene_mapping$ensembl_id,]
row.names(sub_expr) <- gene_mapping$entrez_id

## Combine final Entrez ID expr data
entrez_expr <- rbind(sub_expr, entrez_expr)

write.table(cbind(Entrez = row.names(entrez_expr), entrez_expr), 
            paste0(out_dir, "n2486_Ensembl_dwd_entrez.tsv"), 
            sep = "\t", row.names = F, quote = F)
