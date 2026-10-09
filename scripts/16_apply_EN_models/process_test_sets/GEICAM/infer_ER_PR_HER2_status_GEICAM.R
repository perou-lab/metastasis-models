# Infer ER, PR, HER2 status based on Youden's index for GEICAM (n204)

geicam_expr_path <- ""

library(tidyverse)
library(cutpointr)
library(readxl)

out_dir <- "output/apply_EN_models/process_test_sets/GEICAM/infer_ER_PR_HER2/"
dir.create(out_dir, showWarnings = F, recursive = T)

# Read in GEICAM expression data
geicam_expr <- read.table(geicam_expr_path,
                          sep = "\t", check.names = F, header = T, row.names = 1)
geicam_expr <- geicam_expr[,grepl("AP_", colnames(geicam_expr))]
colnames(geicam_expr) <- sub(".*[YCS]_", "", colnames(geicam_expr))

# Remove gene variants: 
ids <- gsub("\\..*$", "", rownames(geicam_expr))
expr <- rowsum(geicam_expr, group = ids)

# UQN and log2 transform
uqNormalize <- function(expr) {
  uqs <- apply(expr, 2, function(x){ quantile(x[x > 0 & !is.na(x)], 0.75) })
  y <- 1000
  expr <- apply(expr, 1, function(x, y){ x * y }, y / uqs)
  expr <- t(expr)
  return(expr)
}

uqn_expr <- log2(uqNormalize(expr)+1)

# Read in ER, PR, HER2 status for all available samples
## Metadata file was downloaded from Garcia Recio et al. Nature Cancer 2023 (PMC9886551)
metadata <- read_excel("source_files/AURORA_retrospective_Supplementary_Table.2.xlsx", sheet = 6) %>% 
  mutate(jUNCtion_ID = sub("_[A-Z0-9]+$", "", sub("-.*", "", `Sample Name`))) %>% 
  rename(Sample = `Short name`) %>% 
  rename(Set = Study) %>% 
  rename(Tissue_type = `Tissue type`) %>% 
  rename(Sample_type = PvsM) %>% 
  rename(PAM50_call = `PAM50 Call`) %>%
  filter(Tissue_type == "FFPE") %>% 
  rename(Site_original = `Anatomic Site`) %>% 
  mutate(Site_group = ifelse(`Anatomic Site Groups` %in% 
                               c("Brain", "Breast", "Liver",
                                 "Lung", "Lymph node"),
                             `Anatomic Site Groups`, "Other")) %>% 
  mutate(ER_status = case_when(
    ER == "Positive" ~ "Positive", 
    ER == "Negative" ~ "Negative", 
    ER == "NA" ~ NA_character_
  )) %>% 
  mutate(PR_status = case_when(
    PR == "Positive" ~ "Positive", 
    PR == "Negative" ~ "Negative", 
    PR == "NA" ~ NA_character_
  )) %>% 
  mutate(HER2_status = case_when(
    HER2 == "Positive" ~ "Positive", 
    HER2 == "Negative" ~ "Negative", 
    HER2 == "NA" ~ NA_character_
  ))

# AUC / Youden's Index ---------------------------------------------------------
receptor_data <- data.frame(jUNCtion_ID = colnames(uqn_expr), 
                            ER_expr = as.numeric(uqn_expr["ENSG00000091831",]),
                            PR_expr = as.numeric(uqn_expr["ENSG00000082175",]),
                            HER2_expr = as.numeric(uqn_expr["ENSG00000141736",]))
receptor_data <- inner_join(metadata, receptor_data, by = "jUNCtion_ID")

## ER Status -------------------------------------------------------------------
ER_opt_cut <- cutpointr(receptor_data, ER_expr, ER_status,
                        direction = ">=", pos_class = "Positive",
                        neg_class = "Negative", na.rm = TRUE,
                        method = maximize_metric, metric = youden)

png(paste0(out_dir, "GEICAM_n204_opt_cut_ER.png"), width = 8, height = 4, units = 'in', res = 200)
plot(ER_opt_cut)
dev.off()

sink(paste0(out_dir, "GEICAM_n204_opt_cut_ER_summary.txt"))
summary(ER_opt_cut)
sink()

## PR Status -------------------------------------------------------------------
PR_opt_cut <- cutpointr(receptor_data, PR_expr, PR_status,
                        direction = ">=", pos_class = "Positive",
                        neg_class = "Negative", na.rm = TRUE,
                        method = maximize_metric, metric = youden)

png(paste0(out_dir, "GEICAM_n204_opt_cut_PR.png"), width = 8, height = 4, units = 'in', res = 200)
plot(PR_opt_cut)
dev.off()

sink(paste0(out_dir, "GEICAM_n204_opt_cut_PR_summary.txt"))
summary(PR_opt_cut)
sink()

## HER2 Status -------------------------------------------------------------------
HER2_opt_cut <- cutpointr(receptor_data, HER2_expr, HER2_status,
                          direction = ">=", pos_class = "Positive",
                          neg_class = "Negative", na.rm = TRUE,
                          method = maximize_metric, metric = youden)

png(paste0(out_dir, "GEICAM_n204_opt_cut_HER2.png"), width = 8, height = 4, units = 'in', res = 200)
plot(HER2_opt_cut)
dev.off()

sink(paste0(out_dir, "GEICAM_n204_opt_cut_HER2_summary.txt"))
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
  dplyr::select(jUNCtion_ID, ER_expr, PR_expr, HER2_expr,
                ER_status, PR_status, HER2_status,
                inferred_ER_status, inferred_PR_status, inferred_HER2_status,
                ER_status_w_inferred, PR_status_w_inferred, HER2_status_w_inferred, 
                IHC_subtype_w_inferred)

write.table(status, paste0(out_dir, "GEICAM_n204_Ensembl_inferred_clin_groups.txt"), 
            sep = "\t", col.names = T, row.names = F)
