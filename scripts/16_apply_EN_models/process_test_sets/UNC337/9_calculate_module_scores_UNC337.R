# Calculate module scores using Breast-Cancer-Modules R package on UNC337 test set

library(devtools)
library(tidyverse)
load_all("source_files/Breast-Cancer-Modules/")

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/calculate_modules/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Load expression matrix ------------------------------------------------------
unc337 <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_primaries_combined_ensembl_id.tsv", 
                     sep = "\t", header = T, row.names = 1, check.names = F)
set.seed(123)

## Calculate module scores -----------------------------------------------------
modules <- calculate_modules(unc337, 
                             standardize_output = TRUE, 
                             filter_genes = TRUE, 
                             filter_gene_thresh = 0.7, 
                             row_median_center = TRUE, 
                             column_standardize = TRUE, 
                             log_transform = FALSE, 
                             include_median_modules = TRUE, 
                             include_single_genes = TRUE, 
                             include_special_models = TRUE, 
                             include_prognostic = TRUE,
                             input_species = "human")

write.table(modules, paste0(out_dir, "UNC337_n239_primaries_Ensembl_938_published_modules.txt"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T)

## Add ER, PR, HER2 receptor status as features (0/1) --------------------------
receptor_status <- read.table("output/apply_EN_models/process_test_sets/UNC337/infer_ER_PR_HER2/UNC337_n239_Ensembl_dwd_inferred_clin_groups.txt", 
                              sep = "\t", header = T, check.names = F) %>% 
  mutate(ER_status_w_inferred = case_when(
    ER_status_w_inferred == "Positive" ~ 1,
    ER_status_w_inferred == "Negative" ~ 0)) %>% 
  mutate(PR_status_w_inferred = case_when(
    PR_status_w_inferred == "Positive" ~ 1,
    PR_status_w_inferred == "Negative" ~ 0)) %>% 
  mutate(HER2_status_w_inferred = case_when(
    HER2_status_w_inferred == "Positive" ~ 1,
    HER2_status_w_inferred == "Negative" ~ 0)) %>% 
  select(Sample_GEO_accession, ER_status_w_inferred, PR_status_w_inferred, HER2_status_w_inferred) %>% 
  rename(
    clinical_ER_status = ER_status_w_inferred, 
    clinical_PR_status = PR_status_w_inferred, 
    clinical_HER2_status = HER2_status_w_inferred)

modules <- as.data.frame(t(modules)) %>% 
  tibble::rownames_to_column(var = "Sample_GEO_accession")
modules_receptors <- left_join(modules, receptor_status, by="Sample_GEO_accession")
modules_receptors <- tibble::column_to_rownames(modules_receptors, var = "Sample_GEO_accession")
modules_receptors <- as.data.frame(t(modules_receptors))

write.table(modules_receptors, paste0(out_dir, "UNC337_n239_primaries_938_published_modules_3_receptors.txt"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T)
