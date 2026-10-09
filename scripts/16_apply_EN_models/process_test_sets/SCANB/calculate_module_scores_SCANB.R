# Calculate module scores using Breast-Cancer-Modules R package on SCANB test set

out_dir <- "output/apply_EN_models/process_test_sets/SCANB/calculate_modules/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

library(devtools)
library(tidyverse)
library(matrixStats)
library(SummarizedExperiment)
library(MatrixGenerics)
load_all("source_files/Breast-Cancer-Modules/")

## Load expression matrix ------------------------------------------------------
ensembl_expr <- readRDS("../metastatic_site_EN_manuscript/source_files/SCANB_MDS/human_scanb_ensembl_norm.RDS")
expr_mat <- assay(ensembl_expr)
scanb_expr <- as.data.frame(expr_mat)
rm(ensembl_expr, expr_mat)

set.seed(123)

## Calculate module scores -----------------------------------------------------
modules <- calculate_modules(scanb_expr, 
                             standardize_output = TRUE, 
                             filter_genes = TRUE, 
                             filter_gene_thresh = 0.7, 
                             row_median_center = TRUE, 
                             column_standardize = TRUE, 
                             log_transform = TRUE, 
                             include_median_modules = TRUE, 
                             include_single_genes = TRUE, 
                             include_special_models = TRUE, 
                             include_prognostic = TRUE,
                             input_species = "human")

write.table(modules, paste0(out_dir, "SCANB_Ensembl_938_published_modules.txt"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T)

## Add ER, PR, HER2 (IHC) status as features (0/1) -----------------------------
modules <- read.table(paste0(out_dir, "SCANB_Ensembl_938_published_modules.txt"),
                      sep = "\t", header = T, check.names = F)

receptor_status <- read.table("../datasets/SCANB_gitlab/SCAN-B_6329_Clinical_Data.txt",
                       sep = "\t", header = T, check.names = F) %>%
  mutate(
    clinical_ER_status = case_when(
      ER == "Positive" ~ 1, 
      ER == "Negative" ~ 0), 
    clinical_PR_status = case_when(
      PR == "Positive" ~ 1, 
      PR == "Negative" ~ 0), 
    clinical_HER2_status = case_when(
      HER2 == "Positive" ~ 1, 
      HER2 == "Negative" ~ 0)) %>%
  select(Merged, clinical_ER_status, clinical_PR_status, clinical_HER2_status)

modules <- as.data.frame(t(modules)) %>%
  tibble::rownames_to_column(var = "Merged")
modules_receptors <- left_join(modules, receptor_status, by="Merged")
modules_receptors <- tibble::column_to_rownames(modules_receptors, var = "Merged")
modules_receptors <- as.data.frame(t(modules_receptors))

write.table(modules_receptors, paste0(out_dir, "SCANB_Ensembl_938_published_modules_3_receptors.txt"),
            sep = "\t", quote = FALSE, row.names = T, col.names = T)
