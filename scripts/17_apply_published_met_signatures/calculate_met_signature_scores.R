# Calculate signature scores on n2486 dataset using published metastasis signatures

library(tidyverse)
library(readxl)

out_dir <- "output/apply_published_met_signatures/"
dir.create(paste0(out_dir), recursive = T, showWarnings = F)

## Read in expression data
expr <- read.table("output/DWD_adjust/n2486_Ensembl_dwd.txt", 
                   sep = "\t", header = T, check.names = F)

# ------------------------------------------------------------------------------
## Calculate weighted 48 gene lung metastasis score from Minn et al.
symbol_expr <- read.table("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols.tsv", 
                   sep = "\t", header = T, check.names = F, row.names = 1)

## Read in gene list with weights
lung_classifier <- read_excel("source_files/published_met_signatures/STable3_PMC1283098_Minn_Nature_2005.xls",
                              sheet = 5, skip = 1) %>%
  rename(Gene_symbol = `Gene symbol`) %>% 
  rename(t_stat = `t-value`) %>% 
  select(Gene_symbol, t_stat)

weights <- setNames(lung_classifier$t_stat, lung_classifier$Gene_symbol)

## Identify common genes between classifier and expression data
common_genes <- intersect(lung_classifier$Gene_symbol, rownames(symbol_expr))

sub_weights <- weights[common_genes]
sub_expr <- symbol_expr[common_genes,]

## Calculate weighted lung met signature score
lung_met_score <- colSums(sub_expr * sub_weights)
lung_met_score <- lung_met_score * -1

lung_met_score_df <- data.frame(Sample_ID = names(lung_met_score), 
                                Lung_classifier_Minn = as.numeric(lung_met_score))


# ------------------------------------------------------------------------------
# genes from Bos et al. Supplementary Table 4 (PMC2698953); UP genes only
brain_classifier <- c("ANGPTL4", "PLOD2", "COL13A1", "PTGS2", "PELI1", "MMP1",
                      "B4GALT6", "HBEGF", "CSF3", "LTBP1", "FSCN1", "LAMA4")

## Subset expression based on 12 up genes from 17 gene brain classifier and calculate median signature score
sub_expr <- symbol_expr[rownames(symbol_expr) %in% brain_classifier,]
brain_12_gene_score <- as.data.frame(apply(sub_expr, 2, median))
colnames(brain_12_gene_score) <- "Brain_12_gene_Score_Bos_2009"
brain_12_gene_score <- tibble::rownames_to_column(brain_12_gene_score, var = "Sample_ID")


# ------------------------------------------------------------------------------
# 6 genes from Landemaine et al. 2008 (Cancer Research) PMID: 18676831 (Lung Mets)

genes <- c("DSC2", "TFCP2L1", "UGT8", "ITGB8", "ANP32E","FERMT1")

sub_expr <- symbol_expr[rownames(symbol_expr) %in% genes,]
landemaine_lung_score <- as.data.frame(apply(sub_expr, 2, median))
colnames(landemaine_lung_score) <- "Lung_score_Landemaine_2008"
landemaine_lung_score <- tibble::rownames_to_column(landemaine_lung_score, var = "Sample_ID")


# ------------------------------------------------------------------------------
# 5 gene expression signature-based nomogram from Zhao et al. PMC6346244 (Bone Mets)

model <- data.frame(
  genes = c("KRT23", "ALDH3B2", "REEP1", "SPIB", "GLDC"), 
  B_value = c(-1.058, 0.454, 0.513, -0.687, -0.703)
)

sub_expr <- symbol_expr[model$genes,]
sub_expr <- t(sub_expr)
sub_expr <- sub_expr[, model$genes]

Zhao_score <- as.numeric(as.matrix(sub_expr) %*% model$B_value)
names(Zhao_score) <- rownames(sub_expr)
Zhao_score <- as.data.frame(Zhao_score) %>% 
  rename(Zhao_Bone_score = Zhao_score) %>% 
  tibble::rownames_to_column(var = "Sample_ID")


# ------------------------------------------------------------------------------
## Combine scores from 5 signatures and save

met_scores <- lung_met_score_df %>% 
  left_join(brain_12_gene_score, by = "Sample_ID") %>% 
  left_join(landemaine_lung_score, by = "Sample_ID") %>% 
  left_join(Zhao_score, by = "Sample_ID")

write.table(met_scores, paste0(out_dir, "published_met_scores_applied_to_n2486.tsv"), 
            sep = "\t", col.names = T, row.names = F)
