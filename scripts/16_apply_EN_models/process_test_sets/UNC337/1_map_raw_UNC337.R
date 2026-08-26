## Assemble UNC337 dataset from GSE18229
# (files downloaded 4.9.2026)

library(devtools)
library(tidyverse)
library(readxl)
library(HGNChelper)
library(rtracklayer)
load_all("source_files/Breast-Cancer-Modules/")

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/"
dir.create(paste0(out_dir), recursive = T, showWarnings = F)

## Read in 7 matrix files from GSE18229 ----------------------------------------
GPL1390 <- read.table("source_files/UNC337/GSE18229/GSE18229-GPL1390_series_matrix.txt",
                    sep = "\t", header = TRUE, comment.char = "!")

GPL1708 <- read.table("source_files/UNC337/GSE18229/GSE18229-GPL1708_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!")

GPL5325 <- read.table("source_files/UNC337/GSE18229/GSE18229-GPL5325_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!")

GPL6607 <- read.table("source_files/UNC337/GSE18229/GSE18229-GPL6607_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!")

GPL7504 <- read.table("source_files/UNC337/GSE18229/GSE18229-GPL7504_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!")

GPL885 <- read.table("source_files/UNC337/GSE18229/GSE18229-GPL885_series_matrix.txt", 
                             sep = "\t", header = T, comment.char = "!")

GPL887 <- read.table("source_files/UNC337/GSE18229/GSE18229-GPL887_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!")

## Read in 3 mapping files from GSE18229 ---------------------------------------
GPL1708_lines <- readLines("source_files/UNC337/GSE18229/GPL1708_old_annotations.txt")
GPL1708_mapping <- read.table(
  text = GPL1708_lines[!grepl("^(#|\\^)", GPL1708_lines)],
  sep = "\t",
  fill = TRUE,
  header = TRUE,
  quote = "")

GPL885_lines <- readLines("source_files/UNC337/GSE18229/GPL885_old_annotations.txt")
GPL885_mapping <- read.table(
  text = GPL885_lines[!grepl("^(#|\\^)", GPL885_lines)],
  sep = "\t",
  fill = TRUE,
  header = TRUE,
  quote = "")

GPL887_lines <- readLines("source_files/UNC337/GSE18229/GPL887_old_annotations.txt")
GPL887_mapping <- read.table(
  text = GPL887_lines[!grepl("^(#|\\^)", GPL887_lines)],
  sep = "\t",
  fill = TRUE,
  header = TRUE,
  quote = "")

rm(GPL1708_lines, GPL885_lines, GPL887_lines)


## Additional mapping files for other arrays -----------------------------------
# File downloaded on 4/9/2026 from GSE Platform GPL1390:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL1390&id=11300&db=GeoDb_blob92
GPL1390_mapping <- read.table("source_files/UNC337/array_mapping/GPL1390_probes.txt", 
                              sep = "\t", fill = TRUE, header = TRUE, quote = "", comment.char = "#")

# File downloaded on 4/9/2026 from GSE Platform GPL5325
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL5325
GPL5325_mapping <- read.table("source_files/UNC337/array_mapping/GPL5325_probes.txt", 
                              sep = "\t", fill = TRUE, header = TRUE, quote = "", comment.char = "#")

# File downloaded on 4/9/2026 from GSE Platform GPL6607
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL6607&id=8104&db=GeoDb_blob92
GPL6607_mapping <- read.table("source_files/UNC337/array_mapping/GPL6607_probes.txt", 
                              sep = "\t", fill = TRUE, header = TRUE, quote = "", comment.char = "#")

# File downloaded on 4/9/2026 from GSE Platform GPL7504
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL7504
GPL7504_mapping <- read.table("source_files/UNC337/array_mapping/GPL7504-10686.txt", 
                              sep = "\t", fill = TRUE, header = TRUE, quote = "", comment.char = "#")


## GPL1390 PROBE MAPPING -------------------------------------------------------

# Identify probes that are mapped to gene symbols (but not Ensembl IDs)
GPL1390_genes_wo_ensembl <- GPL1390_mapping %>% 
  filter(Ensembl.Gene.ID == "") %>% 
  pull(GENE_NAME)

# Update gene symbols (HGCNHelper)
updated_GPL1390_mapping_no_ensembl <- checkGeneSymbols(GPL1390_genes_wo_ensembl)

# Re-label columns
updated_GPL1390_mapping_no_ensembl <- updated_GPL1390_mapping_no_ensembl %>% 
  dplyr::rename(GENE_NAME = x) %>% 
  dplyr::rename(Updated_Gene_Name = Suggested.Symbol)

# Map new gene symbols back to the original map (using GENE_NAME)
GPL1390_genes_wo_ensembl <- left_join(GPL1390_mapping, updated_GPL1390_mapping_no_ensembl, 
                                      by = "GENE_NAME", relationship = "many-to-many") %>% 
  filter(Approved != "FALSE") %>% 
  dplyr::rename(provided_gene_name = GENE_NAME) %>% 
  dplyr::rename(gene_name = Updated_Gene_Name)

# Load in gene map (symbol to ensembl id)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id)

# Map new symbols to ensembl gene ids (with gene_map_ensembl)
GPL1390_map_new_ensembl <- left_join(GPL1390_genes_wo_ensembl, gene_map_ensembl, by = "gene_name") %>% 
  filter(!is.na(ensembl_id)) %>% 
  select(ID, ensembl_id) %>% 
  mutate(ensembl_id = as.character(ensembl_id))

# Select just probe IDs that map to ensembl ids
GPL1390_mapping_ensembl <- GPL1390_mapping %>% 
  select(ID, Ensembl.Gene.ID) %>% 
  left_join(GPL1390_map_new_ensembl, by = "ID") %>% 
  mutate(Ensembl.Gene.ID = na_if(Ensembl.Gene.ID, "")) %>% 
  mutate(final_ensembl_id = coalesce(Ensembl.Gene.ID, ensembl_id)) %>% 
  select(ID, final_ensembl_id) %>% 
  dplyr::rename(ensembl_id = final_ensembl_id) %>% 
  dplyr::rename(ID_REF = ID) %>% 
  distinct()

# Map gene IDs (form of numerical 1-N) to final ensembl ID list
GPL1390 <- right_join(GPL1390_mapping_ensembl, GPL1390, by = "ID_REF") %>% 
  filter(!is.na(ensembl_id)) %>% 
  select(-ID_REF)

write.table(GPL1390, paste0(out_dir, "GPL1390_mapped_ensembl_id.tsv"), 
            sep = "\t", row.names = F, col.names = T)

rm(GPL1390_genes_wo_ensembl, GPL1390_map_new_ensembl, GPL1390_mapping, updated_GPL1390_mapping_no_ensembl, GPL1390_mapping_ensembl, GPL1390)


## GPL1708 PROBE MAPPING -------------------------------------------------------

# Identify probes that are mapped to gene symbols (but not Ensembl transcripts IDs)
GPL1708_genes_wo_ensembl <- GPL1708_mapping %>% 
  filter(is.na(ENSEMBL_ID)) %>% 
  pull(GENE_SYMBOL)
## There are no genes uniquely mapped to GENE SYMBOLS

# Load gencode annotations
gtf.gr <- rtracklayer::import("/datastore/lbcfs/labs/bioinformatics/references/GRCh38.d1.vd1/star-2.7.6a-gencode.v36/gencode.v36.annotation.gtf")
gtf.df <- as.data.frame(gtf.gr) %>% 
  dplyr::select(transcript_id, gene_id) %>% 
  mutate(transcript_id = sub("\\..*", "", transcript_id), 
         gene_id = sub("\\..*", "", gene_id)) %>% 
  dplyr::rename(ENSEMBL_ID = transcript_id) %>% 
  dplyr::rename(ENSEMBL_gene_id = gene_id) %>% 
  filter(!is.na(ENSEMBL_ID)) %>% 
  distinct()

# Map probes (Ensembl transcripts) to Ensembl gene IDs
GPL1708_mapping <- left_join(GPL1708_mapping, gtf.df, by = "ENSEMBL_ID")

# For ENST that did not map to ENSG, use the gene symbol: 
GPL1708_genes <- GPL1708_mapping %>% 
  filter(is.na(ENSEMBL_gene_id) & !is.na(GENE_SYMBOL)) %>% 
  filter(GENE_SYMBOL != "") %>% 
  dplyr::select(ID, GENE_SYMBOL)

# Update gene symbols (HGCNHelper)
updated_GPL1708_genes <- checkGeneSymbols(GPL1708_genes$GENE_SYMBOL)

# Re-label columns
updated_GPL1708_genes <- updated_GPL1708_genes %>% 
  dplyr::rename(GENE_SYMBOL = x) %>% 
  dplyr::rename(Updated_Gene_Symbol = Suggested.Symbol)

# Map new gene symbols back to the original map (using GENE_NAME)
GPL1708_genes_wo_ensembl <- left_join(GPL1708_genes, updated_GPL1708_genes, 
                                      by = "GENE_SYMBOL", relationship = "many-to-many") %>% 
  filter(Approved != "FALSE") %>% 
  dplyr::rename(provided_gene_symbol = GENE_SYMBOL) %>% 
  dplyr::rename(gene_symbol = Updated_Gene_Symbol) %>% 
  distinct()

# Load in gene map (symbol to ensembl id)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id) %>% 
  dplyr::rename(gene_symbol = gene_name)

# Map new symbols to ensembl gene ids (with gene_map_ensembl)
GPL1708_map_symb_ensembl <- left_join(GPL1708_genes_wo_ensembl, gene_map_ensembl, by = "gene_symbol") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(ID, ensembl_id) %>% 
  mutate(ensembl_id = as.character(ensembl_id)) %>% 
  distinct()

# Combine symbol map with ensembl map
GPL1708_final_map <- left_join(GPL1708_mapping, GPL1708_map_symb_ensembl, by = "ID") %>% 
  mutate(final_ensembl_id = coalesce(ENSEMBL_gene_id, ensembl_id)) %>% 
  dplyr::select(ID, final_ensembl_id) %>% 
  dplyr::rename(ID_REF = ID, 
                ensembl_id = final_ensembl_id) %>% 
  distinct()

# Map gene IDs (form of numerical 1-N) to final ensembl ID list
GPL1708 <- right_join(GPL1708_final_map, GPL1708, by = "ID_REF") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-ID_REF)

write.table(GPL1708, paste0(out_dir, "GPL1708_mapped_ensembl_id.tsv"), 
            sep = "\t", row.names = F, col.names = T)

rm(GPL1708_gene_map, GPL1708_genes, GPL1708_mapping, GPL1708_map_symb_ensembl, GPL1708_genes_wo_ensembl, updated_GPL1708_genes, GPL1708_genes_wo_ensembl, GPL1708_lines)


## GPL5325 PROBE MAPPING -------------------------------------------------------

# Identify probes that are mapped to gene symbols
GPL5325_genes <- GPL5325_mapping %>% 
  filter(Blast.Gene.Symbol != "") %>% 
  pull(Blast.Gene.Symbol)

# Update gene symbols (HGCNHelper)
updated_GPL5325_genes <- checkGeneSymbols(GPL5325_genes)

# Re-label columns
updated_GPL5325_genes <- updated_GPL5325_genes %>% 
  dplyr::rename(Blast.Gene.Symbol = x) %>% 
  dplyr::rename(Updated_Gene_Symbol = Suggested.Symbol)

# Map new gene symbols back to the original map (using GENE_NAME)
GPL5325_mapping <- left_join(GPL5325_mapping, updated_GPL5325_genes,
                             by = "Blast.Gene.Symbol", relationship = "many-to-many") %>% 
  filter(Approved != "FALSE") %>% 
  dplyr::rename(provided_gene_name = Blast.Gene.Symbol) %>% 
  dplyr::rename(gene_name = Updated_Gene_Symbol)

# Load in gene map (symbol to ensembl id)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id)

# Map new symbols to ensembl gene ids (with gene_map_ensembl)
GPL5325_mapping_ensembl <- left_join(GPL5325_mapping, gene_map_ensembl, by = "gene_name") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(ID, ensembl_id) %>% 
  mutate(ensembl_id = as.character(ensembl_id)) %>% 
  dplyr::rename(ID_REF = ID) %>% 
  distinct()

# Map gene IDs (form of numerical 1-N) to final ensembl ID list
GPL5325 <- right_join(GPL5325_mapping_ensembl, GPL5325, by = "ID_REF") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-ID_REF)

write.table(GPL5325, paste0(out_dir, "GPL5325_mapped_ensembl_id.tsv"), 
            sep = "\t", row.names = F, col.names = T)

rm(GPL5325_mapping_ensembl, GPL5325_genes, updated_GPL5325_genes)


## GPL6607 PROBE MAPPING -------------------------------------------------------

# Identify probes that are mapped to gene symbols
GPL6607_genes <- GPL6607_mapping %>% 
  filter(Blast.Gene.Symbol != "") %>% 
  pull(Blast.Gene.Symbol)

# Update gene symbols (HGCNHelper)
updated_GPL6607_genes <- checkGeneSymbols(GPL6607_genes)

# Re-label columns
updated_GPL6607_genes <- updated_GPL6607_genes %>% 
  dplyr::rename(Blast.Gene.Symbol = x) %>% 
  dplyr::rename(Updated_Gene_Symbol = Suggested.Symbol)

# Map new gene symbols back to the original map (using GENE_NAME)
GPL6607_mapping <- left_join(GPL6607_mapping, updated_GPL6607_genes,
                             by = "Blast.Gene.Symbol", relationship = "many-to-many") %>% 
  filter(Approved != "FALSE") %>% 
  dplyr::rename(provided_gene_name = Blast.Gene.Symbol) %>% 
  dplyr::rename(gene_name = Updated_Gene_Symbol)

# Load in gene map (symbol to ensembl id)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id)

# Map new symbols to ensembl gene ids (with gene_map_ensembl)
GPL6607_mapping_ensembl <- left_join(GPL6607_mapping, gene_map_ensembl, by = "gene_name") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(ID, ensembl_id) %>% 
  mutate(ensembl_id = as.character(ensembl_id)) %>% 
  dplyr::rename(ID_REF = ID) %>% 
  distinct()

# Map gene IDs (form of numerical 1-N) to final ensembl ID list
GPL6607 <- right_join(GPL6607_mapping_ensembl, GPL6607, by = "ID_REF") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-ID_REF)

write.table(GPL6607, paste0(out_dir, "GPL6607_mapped_ensembl_id.tsv"), 
            sep = "\t", row.names = F, col.names = T)

rm(GPL6607_mapping_ensembl, GPL6607_genes, updated_GPL6607_genes)


## GPL7504 PROBE MAPPING -------------------------------------------------------

# Identify probes that are mapped to gene symbols
GPL7504_genes <- GPL7504_mapping %>% 
  filter(`Gene.Symbol` != "") %>% 
  pull(`Gene.Symbol`)

# Update gene symbols (HGCNHelper)
updated_GPL7504_genes <- checkGeneSymbols(GPL7504_genes)

# Re-label columns
updated_GPL7504_genes <- updated_GPL7504_genes %>% 
  dplyr::rename(Gene.Symbol = x) %>% 
  dplyr::rename(Updated_Gene_Symbol = Suggested.Symbol)

# Map new gene symbols back to the original map (using GENE_NAME)
GPL7504_mapping <- left_join(GPL7504_mapping, updated_GPL7504_genes,
                             by = "Gene.Symbol", relationship = "many-to-many") %>% 
  filter(Approved != "FALSE") %>% 
  dplyr::rename(provided_gene_name = Gene.Symbol) %>% 
  dplyr::rename(gene_name = Updated_Gene_Symbol)

# Load in gene map (symbol to ensembl id)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id)

# Map new symbols to ensembl gene ids (with gene_map_ensembl)
GPL7504_mapping_ensembl <- left_join(GPL7504_mapping, gene_map_ensembl, by = "gene_name") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(ID, ensembl_id) %>% 
  mutate(ensembl_id = as.character(ensembl_id)) %>% 
  dplyr::rename(ID_REF = ID) %>% 
  distinct()

# Map gene IDs (form of numerical 1-N) to final ensembl ID list
GPL7504 <- right_join(GPL7504_mapping_ensembl, GPL7504, by = "ID_REF") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-ID_REF)

write.table(GPL7504, paste0(out_dir, "GPL7504_mapped_ensembl_id.tsv"), 
            sep = "\t", row.names = F, col.names = T)

rm(GPL7504_mapping_ensembl, GPL7504_genes, updated_GPL7504_genes)


## GPL885 PROBE MAPPING --------------------------------------------------------

# Identify probes that are mapped to gene symbols (but not Ensembl transcripts IDs)
GPL885_genes_wo_ensembl <- GPL885_mapping %>% 
  filter(is.na(ENSEMBL_ID)) %>% 
  pull(GENE_SYMBOL)
## There are no genes uniquely mapped to GENE SYMBOLS

# Load gencode annotations
gtf.gr <- rtracklayer::import("/datastore/lbcfs/labs/bioinformatics/references/GRCh38.d1.vd1/star-2.7.6a-gencode.v36/gencode.v36.annotation.gtf")
gtf.df <- as.data.frame(gtf.gr) %>% 
  dplyr::select(transcript_id, gene_id) %>% 
  mutate(transcript_id = sub("\\..*", "", transcript_id), 
         gene_id = sub("\\..*", "", gene_id)) %>% 
  dplyr::rename(ENSEMBL_ID = transcript_id) %>% 
  dplyr::rename(ENSEMBL_gene_id = gene_id) %>% 
  filter(!is.na(ENSEMBL_ID)) %>% 
  distinct()

# Map probes (Ensembl transcripts) to Ensembl gene IDs
GPL885_mapping <- left_join(GPL885_mapping, gtf.df, by = "ENSEMBL_ID")

# For ENST that did not map to ENSG, use the gene symbol: 
GPL885_genes <- GPL885_mapping %>% 
  filter(is.na(ENSEMBL_gene_id) & !is.na(GENE_SYMBOL)) %>% 
  filter(GENE_SYMBOL != "") %>% 
  dplyr::select(ID, GENE_SYMBOL)

# Update gene symbols (HGCNHelper)
updated_GPL885_genes <- checkGeneSymbols(GPL885_genes$GENE_SYMBOL)

# Re-label columns
updated_GPL885_genes <- updated_GPL885_genes %>% 
  dplyr::rename(GENE_SYMBOL = x) %>% 
  dplyr::rename(Updated_Gene_Symbol = Suggested.Symbol)

# Map new gene symbols back to the original map (using GENE_NAME)
GPL885_genes_wo_ensembl <- left_join(GPL885_genes, updated_GPL885_genes, 
                                     by = "GENE_SYMBOL", relationship = "many-to-many") %>% 
  filter(Approved != "FALSE") %>% 
  dplyr::rename(provided_gene_symbol = GENE_SYMBOL) %>% 
  dplyr::rename(gene_symbol = Updated_Gene_Symbol) %>% 
  distinct()

# Load in gene map (symbol to ensembl id)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id) %>% 
  dplyr::rename(gene_symbol = gene_name)

# Map new symbols to ensembl gene ids (with gene_map_ensembl)
GPL885_map_symb_ensembl <- left_join(GPL885_genes_wo_ensembl, gene_map_ensembl, by = "gene_symbol") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(ID, ensembl_id) %>% 
  mutate(ensembl_id = as.character(ensembl_id)) %>% 
  distinct()

# Combine symbol map with ensembl map
GPL885_final_map <- left_join(GPL885_mapping, GPL885_map_symb_ensembl, by = "ID") %>% 
  mutate(final_ensembl_id = coalesce(ENSEMBL_gene_id, ensembl_id)) %>% 
  dplyr::select(ID, final_ensembl_id) %>% 
  dplyr::rename(ID_REF = ID, 
                ensembl_id = final_ensembl_id) %>% 
  distinct()

# Map gene IDs (form of numerical 1-N) to final ensembl ID list
GPL885 <- right_join(GPL885_final_map, GPL885, by = "ID_REF") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-ID_REF)

write.table(GPL885, paste0(out_dir, "GPL885_mapped_ensembl_id.tsv"), 
            sep = "\t", row.names = F, col.names = T)

rm(GPL885_genes, GPL885_mapping, GPL885_map_symb_ensembl, GPL885_genes_wo_ensembl, updated_GPL885_genes)


## GPL887 PROBE MAPPING --------------------------------------------------------

# Identify probes that are mapped to gene symbols (but not Ensembl transcripts IDs)
GPL887_genes_wo_ensembl <- GPL887_mapping %>% 
  filter(is.na(ENSEMBL_ID)) %>% 
  pull(GENE_SYMBOL)
## There are no genes uniquely mapped to GENE SYMBOLS

# Load gencode annotations
gtf.gr <- rtracklayer::import("/datastore/lbcfs/labs/bioinformatics/references/GRCh38.d1.vd1/star-2.7.6a-gencode.v36/gencode.v36.annotation.gtf")
gtf.df <- as.data.frame(gtf.gr) %>% 
  dplyr::select(transcript_id, gene_id) %>% 
  mutate(transcript_id = sub("\\..*", "", transcript_id), 
         gene_id = sub("\\..*", "", gene_id)) %>% 
  dplyr::rename(ENSEMBL_ID = transcript_id) %>% 
  dplyr::rename(ENSEMBL_gene_id = gene_id) %>% 
  filter(!is.na(ENSEMBL_ID)) %>% 
  distinct()

# Map probes (Ensembl transcripts) to Ensembl gene IDs
GPL887_mapping <- left_join(GPL887_mapping, gtf.df, by = "ENSEMBL_ID")

# For ENST that did not map to ENSG, use the gene symbol: 
GPL887_genes <- GPL887_mapping %>% 
  filter(is.na(ENSEMBL_gene_id) & !is.na(GENE_SYMBOL)) %>% 
  filter(GENE_SYMBOL != "") %>% 
  dplyr::select(ID, GENE_SYMBOL)

# Update gene symbols (HGCNHelper)
updated_GPL887_genes <- checkGeneSymbols(GPL887_genes$GENE_SYMBOL)

# Re-label columns
updated_GPL887_genes <- updated_GPL887_genes %>% 
  dplyr::rename(GENE_SYMBOL = x) %>% 
  dplyr::rename(Updated_Gene_Symbol = Suggested.Symbol)

# Map new gene symbols back to the original map (using GENE_NAME)
GPL887_genes_wo_ensembl <- left_join(GPL887_genes, updated_GPL887_genes, 
                                      by = "GENE_SYMBOL", relationship = "many-to-many") %>% 
  filter(Approved != "FALSE") %>% 
  dplyr::rename(provided_gene_symbol = GENE_SYMBOL) %>% 
  dplyr::rename(gene_symbol = Updated_Gene_Symbol) %>% 
  distinct()

# Load in gene map (symbol to ensembl id)
gene_map_ensembl <- gene_map_human %>% 
  dplyr::select(gene_name, ensembl_id) %>% 
  dplyr::rename(gene_symbol = gene_name)

# Map new symbols to ensembl gene ids (with gene_map_ensembl)
GPL887_map_symb_ensembl <- left_join(GPL887_genes_wo_ensembl, gene_map_ensembl, by = "gene_symbol") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(ID, ensembl_id) %>% 
  mutate(ensembl_id = as.character(ensembl_id)) %>% 
  distinct()

# Combine symbol map with ensembl map
GPL887_final_map <- left_join(GPL887_mapping, GPL887_map_symb_ensembl, by = "ID") %>% 
  mutate(final_ensembl_id = coalesce(ENSEMBL_gene_id, ensembl_id)) %>% 
  dplyr::select(ID, final_ensembl_id) %>% 
  dplyr::rename(ID_REF = ID, 
                ensembl_id = final_ensembl_id) %>% 
  distinct()

# Map gene IDs (form of numerical 1-N) to final ensembl ID list
GPL887 <- right_join(GPL887_final_map, GPL887, by = "ID_REF") %>% 
  filter(!is.na(ensembl_id)) %>% 
  dplyr::select(-ID_REF)

write.table(GPL887, paste0(out_dir, "GPL887_mapped_ensembl_id.tsv"), 
            sep = "\t", row.names = F, col.names = T)

