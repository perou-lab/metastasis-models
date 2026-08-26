# Calculate module scores using Breast-Cancer-Modules R package

out_dir <- "output/calculate_modules/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

library(devtools)
library(tidyverse)
load_all("/datastore/nextgenout4/share/labs/peroulab/Breast-Cancer-Modules/")

## Load dwd expression matrix --------------------------------------------------
dwd <- read.table("output/DWD_adjust/n2486_Ensembl_dwd.txt",
                  sep = "\t", header = T, row.names = 1, check.names = F)

set.seed(123)

## Calculate module scores -----------------------------------------------------
modules <- calculate_modules(dwd, 
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

write.table(modules, paste0(out_dir, "n2486_Ensembl_dwd_938_published_modules.txt"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T)

## Add ER, PR, HER2 (IHC) status as features (0/1) -----------------------------
receptor_status <- read.table("output/infer_ER_PR_HER2/n2486_Ensembl_dwd_inferred_clin_groups.txt", 
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
  select(Sample_ID, ER_status_w_inferred, PR_status_w_inferred, HER2_status_w_inferred) %>% 
  rename(
    clinical_ER_status = ER_status_w_inferred, 
    clinical_PR_status = PR_status_w_inferred, 
    clinical_HER2_status = HER2_status_w_inferred)
  
modules <- as.data.frame(t(modules)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")
modules_receptors <- left_join(modules, receptor_status, by="Sample_ID")
modules_receptors <- tibble::column_to_rownames(modules_receptors, var = "Sample_ID")
modules_receptors <- as.data.frame(t(modules_receptors))

write.table(modules_receptors, paste0(out_dir, "n2486_Ensembl_dwd_941_published_modules_3_receptors.txt"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T)

## Remove 6 prognostic signatures ----------------------------------------------
modules_to_remove <- c("ROR_P_Model_JCO.2009_PMID.19204204", 
                       "GHI_RS_Model_NEJM.2004_PMID.15591335", 
                       "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860", 
                       "TNBC_Clinically_Relevant_Good.26_BCR.2011_PMID.21978456", 
                       "TNBC_Clinically_Relevant_Poor.230_BCR.2011_PMID.21978456", 
                       "TNBC_Clinically_Relevant_Poor.26_BCR.2011_PMID.21978456")

modules <- read.table(paste0(out_dir, "n2486_Ensembl_dwd_941_published_modules_3_receptors.txt"), 
                      sep = "\t", header = T) %>% 
  rownames_to_column(var = "rownames") %>% 
  filter(!rownames %in% modules_to_remove) %>%
  column_to_rownames(var = "rownames")

write.table(modules, paste0(out_dir, "n2486_Ensembl_dwd_935_published_modules_3_receptors.txt"),
            sep = "\t", quote = FALSE, row.names = T, col.names = T)


