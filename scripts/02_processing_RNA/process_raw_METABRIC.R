# Processing raw METABRIC data

library(tidyverse)
library(impute)
library(HGNChelper) # for modernizing gene symbols
source("source_files/functions.R")
load("source_files/Breast-Cancer-Modules/data/gene_map_human.rda")
set.seed(123)

out_dir <- "output/processing_RNA/METABRIC"
dir.create(out_dir, showWarnings = F, recursive = T)

## METABRIC expression data was obtained from cBioPortal on 9/2/2024 at the following link: 
# https://www.cbioportal.org/study/summary?id=brca_metabric

MB <- read.table("source_files/METABRIC (cBioPortal)/data_mrna_illumina_microarray.txt", 
                 sep = "\t", header = T)

## Map original gene symbols to ENSEMBL IDs

# Modernized gene symbols as needed
MB_genes <- MB %>% 
  filter(Hugo_Symbol != "NA") %>% 
  dplyr::select(-Entrez_Gene_Id) %>% 
  dplyr::rename(gene_name = Hugo_Symbol) %>% 
  dplyr::select(gene_name, everything())

MB_genes_input <- MB_genes$gene_name
updated_mapping <- checkGeneSymbols(MB_genes_input)
updated_mapping <- updated_mapping %>% 
  dplyr::rename(updated_name = Suggested.Symbol)

mapped_genes <- cbind(updated_mapping, MB_genes)
mapped_genes <- mapped_genes %>% 
  filter(Approved == "TRUE") %>% 
  dplyr::select(-x, -gene_name, -Approved) %>% 
  dplyr::rename(gene_name = updated_name)

# Map to ENSEMBL IDs (in Matthew's file)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id)

MB_ensembl <- right_join(gene_map_ensembl, mapped_genes,
                      by = "gene_name", relationship = "many-to-many")

MB_ensembl <- MB_ensembl %>% 
  dplyr::select(-gene_name) %>% 
  drop_na(ensembl_id)

write.table(MB_ensembl, paste0(out_dir, "/METABRIC_Ensembl.txt"), 
            sep = "\t", row.names = F, col.names = T)

# Collapse based on duplicated ENSEMBL IDs
file <- paste0(out_dir, "/METABRIC_Ensembl.txt")
x <- readarray(file, hr=1, impute=F)
write.table(signif(x$xd,4), paste0(out_dir, "/METABRIC_Ensembl_clpsd.txt"), 
            sep = "\t", col.names = NA)

MB <- read.table(paste0(out_dir, "/METABRIC_Ensembl_clpsd.txt"), 
                 sep = "\t", header = T, row.names = 1)

# Filter genes that are missing 30% of data:
MB <- MB[rowSums(!is.na(MB)) >= 0.7*ncol(MB),]

# Impute missing data (k=10, default)
MB <- impute.knn(as.matrix(MB))
MB <- MB$data

# Save final processed matrix with ENSEMBL:
write.table(MB, paste0(out_dir, "/METABRIC_Ensembl_clpsd_filt_imp.txt"),
            sep = "\t", row.names = T, col.names= T)
