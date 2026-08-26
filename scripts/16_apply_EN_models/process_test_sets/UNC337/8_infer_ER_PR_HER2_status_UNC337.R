# Infer ER, PR, HER2 status based on Youden's index for UNC337 (n239 primaries)

library(tidyverse)
library(cutpointr)
library(readxl)

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/infer_ER_PR_HER2/"
dir.create(out_dir, showWarnings = F, recursive = T)

# Read in UNC337 expression data
unc_expr <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_primaries_combined_ensembl_id.tsv", 
                       sep = "\t", header = T, check.names = F)

# Read in ER, PR, HER2 status for all available samples
metadata <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_metadata.tsv", 
                       sep = "\t", header = T, check.names = F) %>% 
  mutate(ER_status = case_when(
    ER_status == 1 ~ "Positive", 
    ER_status == 0 ~ "Negative")) %>% 
  mutate(PR_status = case_when(
    PR_status == 1 ~ "Positive",
    PR_status == 0 ~ "Negative")) %>% 
  mutate(HER2_status = case_when(
    HER2_status == 1 ~ "Positive", 
    HER2_status == 0 ~ "Negative"))

# AUC / Youden's Index ---------------------------------------------------------
receptor_data <- data.frame(Sample_GEO_accession = colnames(unc_expr), 
                            ER_expr = as.numeric(unc_expr["ENSG00000091831",]),
                            PR_expr = as.numeric(unc_expr["ENSG00000082175",]),
                            HER2_expr = as.numeric(unc_expr["ENSG00000141736",]))
receptor_data <- inner_join(metadata, receptor_data, by = "Sample_GEO_accession")

## ER Status -------------------------------------------------------------------
ER_opt_cut <- cutpointr(receptor_data, ER_expr, ER_status,
                        direction = ">=", pos_class = "Positive",
                        neg_class = "Negative", na.rm = TRUE,
                        method = maximize_metric, metric = youden)

png(paste0(out_dir, "UNC337_n239_opt_cut_ER.png"), width = 8, height = 4, units = 'in', res = 200)
plot(ER_opt_cut)
dev.off()

sink(paste0(out_dir, "UNC337_n239_opt_cut_ER_summary.txt"))
summary(ER_opt_cut)
sink()

## PR Status -------------------------------------------------------------------
PR_opt_cut <- cutpointr(receptor_data, PR_expr, PR_status,
                        direction = ">=", pos_class = "Positive",
                        neg_class = "Negative", na.rm = TRUE,
                        method = maximize_metric, metric = youden)

png(paste0(out_dir, "UNC337_n239_opt_cut_PR.png"), width = 8, height = 4, units = 'in', res = 200)
plot(PR_opt_cut)
dev.off()

sink(paste0(out_dir, "UNC337_n239_opt_cut_PR_summary.txt"))
summary(PR_opt_cut)
sink()

## HER2 Status -------------------------------------------------------------------
HER2_opt_cut <- cutpointr(receptor_data, HER2_expr, HER2_status,
                          direction = ">=", pos_class = "Positive",
                          neg_class = "Negative", na.rm = TRUE,
                          method = maximize_metric, metric = youden)

png(paste0(out_dir, "UNC337_n239_opt_cut_HER2.png"), width = 8, height = 4, units = 'in', res = 200)
plot(HER2_opt_cut)
dev.off()

sink(paste0(out_dir, "UNC337_n239_opt_cut_HER2_summary.txt"))
summary(HER2_opt_cut)
sink()

# Add receptor status to metadata
status <- receptor_data %>% 
  mutate(inferred_ER_status = ifelse(ER_expr >= ER_opt_cut$optimal_cutpoint, "Positive", "Negative"), 
         inferred_PR_status = ifelse(PR_expr >= PR_opt_cut$optimal_cutpoint, "Positive", "Negative"), 
         inferred_HER2_status = ifelse(HER2_expr >= HER2_opt_cut$optimal_cutpoint, "Positive", "Negative"), 
         ER_status_w_inferred = coalesce(ER_status, inferred_ER_status), 
         PR_status_w_inferred = coalesce(PR_status, inferred_PR_status), 
         HER2_status_w_inferred = coalesce(HER2_status, inferred_HER2_status), 
         IHC_subtype_w_inferred = case_when(
           ER_status_w_inferred == "Positive" & PR_status_w_inferred == "Positive" & HER2_status_w_inferred == "Negative" ~ "HR+/HER2-",
           ER_status_w_inferred == "Positive" & PR_status_w_inferred == "Negative" & HER2_status_w_inferred == "Negative" ~ "HR+/HER2-", 
           ER_status_w_inferred == "Negative" & PR_status_w_inferred == "Positive" & HER2_status_w_inferred == "Negative" ~ "HR+/HER2-",
           
           ER_status_w_inferred == "Negative" & PR_status_w_inferred == "Positive" & HER2_status_w_inferred == "Positive" ~ "HR+/HER2+",
           ER_status_w_inferred == "Positive" & PR_status_w_inferred == "Negative" & HER2_status_w_inferred == "Positive" ~ "HR+/HER2+",
           ER_status_w_inferred == "Positive" & PR_status_w_inferred == "Positive" & HER2_status_w_inferred == "Positive" ~ "HR+/HER2+",
           
           ER_status_w_inferred == "Negative" & PR_status_w_inferred == "Negative" & HER2_status_w_inferred == "Negative" ~ "TNBC", 
           ER_status_w_inferred == "Negative" & PR_status_w_inferred == "Negative" & HER2_status_w_inferred == "Positive" ~ "HR-/HER2+"
         )) %>% 
  dplyr::select(Sample_GEO_accession, ER_expr, PR_expr, HER2_expr,
                ER_status, PR_status, HER2_status,
                inferred_ER_status, inferred_PR_status, inferred_HER2_status,
                ER_status_w_inferred, PR_status_w_inferred, HER2_status_w_inferred, 
                IHC_subtype_w_inferred)

write.table(status, paste0(out_dir, "UNC337_n239_Ensembl_dwd_inferred_clin_groups.txt"), 
            sep = "\t", col.names = T, row.names = F)
