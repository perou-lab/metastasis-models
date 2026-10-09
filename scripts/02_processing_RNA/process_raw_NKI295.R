# Processing raw NKI295 data

library(tidyverse)
library(impute)
library(seventyGeneData)
library(HGNChelper) # for modernizing gene symbols
source("source_files/functions.R")
library(devtools)
load_all("source_files/Breast-Cancer-Modules/")
set.seed(123)

out_dir <- "output/processing_RNA/NKI295/"
dir.create(out_dir, showWarnings = F, recursive = T)

# Load in NKI data stored in R packages
data(vanDeVijver)
vdv_expr <- vanDeVijver@assayData$exprs
vdv_expr <- tibble::rownames_to_column(vdv_expr, var = "Substance")

## Convert probes to ENTREZ ID or GENE SYMBOLs if no Entrez is listed ----------
# Starting with all genes that have ENTREZ ID mappings
probe_meta <- vanDeVijver@featureData@data
probe_meta <- probe_meta %>% 
  dplyr::rename(Entrez_ID = EntrezGene.ID) %>% 
  dplyr::rename(Gene_Symbol = Gene) %>% 
  mutate(Entrez_ID = gsub(" ", "", Entrez_ID)) %>% 
  dplyr::select(Substance, Entrez_ID, NCBI.gene.symbol, Gene_Symbol)

entrez_list <- probe_meta %>% 
  filter(!is.na(Entrez_ID)) %>% 
  dplyr::select(Substance, Entrez_ID)

# Modernize ENTREZ ID list 
entrez_list <- entrez_list %>% 
  mutate(entrez_id = as.integer(modernize_entrez_ids(Entrez_ID))) %>% 
  dplyr::select(-Entrez_ID)

# Map to ENSEMBL IDs (in Matthew's file)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(entrez_id, ensembl_id)

probe_entrez_ensembl <- right_join(gene_map_ensembl, entrez_list,
                            by = "entrez_id", relationship = "many-to-many")

probe_entrez_ensembl <- probe_entrez_ensembl %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-entrez_id)

## ONLY GENE SYMBOLS -----------------------------------------------------------
# 1750 additional genes that don't have Entrez ID/Gene Symbols
gene_symbol_extras <- probe_meta %>% 
  filter(is.na(Entrez_ID), 
         !is.na(Gene_Symbol)) %>% 
  dplyr::select(-Entrez_ID, -NCBI.gene.symbol)

# Update/modernize gene symbols
NKI_genes_input <- gene_symbol_extras$Gene_Symbol
updated_mapping <- checkGeneSymbols(NKI_genes_input)
## A bunch of genes have multiple mappings
updated_mapping <- updated_mapping %>% 
  rename(updated_name = Suggested.Symbol)

probe_updated_symbol <- cbind(gene_symbol_extras, updated_mapping) %>%
  filter(Approved == "TRUE") %>% 
  dplyr::select(Substance, updated_name) %>% 
  dplyr::rename(gene_name = updated_name)

# Map to ENSEMBL IDs (in Matthew's file)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id)

probe_symbol_ensembl <- right_join(gene_map_ensembl, probe_updated_symbol,
                                   by = "gene_name", relationship = "many-to-many")

probe_symbol_ensembl <- probe_symbol_ensembl %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-gene_name)

probe_mapping <- rbind(probe_symbol_ensembl, probe_entrez_ensembl)

## Map ENSEMBL ID back to probes (Substance)
vdv_expr <- merge(probe_mapping, vdv_expr, by = "Substance")
vdv_expr <- vdv_expr %>% 
  dplyr::select(-Substance) %>% 
  drop_na(ensembl_id)

write.table(vdv_expr, paste0(out_dir, "/NKI295_Ensembl.txt"), 
            sep = "\t", row.names = F, col.names = T)

## Collapse based on duplicated ENSEMBL IDs
file <- paste0(out_dir, "/NKI295_Ensembl.txt")
x <- readarray(file, hr=1, impute=F)
write.table(signif(x$xd,4), paste0(out_dir, "/NKI295_Ensembl_clpsd.txt"), 
            sep = "\t", col.names = NA)

NKI295 <- read.table(paste0(out_dir, "/NKI295_Ensembl_clpsd.txt"), 
                 sep = "\t", header = T, row.names = 1)

## Filter genes that are missing 30% of data:
NKI295 <- NKI295[rowSums(!is.na(NKI295)) >= 0.7*ncol(NKI295),]

## Impute missing data (k=10, default)
NKI295 <- impute.knn(as.matrix(NKI295))
NKI295 <- NKI295$data

## Save final processed matrix with ENSEMBL:
write.table(NKI295, paste0(out_dir, "/NKI295_Ensembl_clpsd_filt_imp.txt"),
            sep = "\t", row.names = T, col.names= T)