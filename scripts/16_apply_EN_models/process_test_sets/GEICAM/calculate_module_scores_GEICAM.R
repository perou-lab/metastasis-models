# Calculate module scores using Breast-Cancer-Modules R package on GEICAM test set

library(devtools)
library(tidyverse)
load_all("source_files/Breast-Cancer-Modules/")

out_dir <- "output/apply_EN_models/process_test_sets/GEICAM/calculate_modules/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Load expression matrix ------------------------------------------------------
geicam_meta <- read_excel("source_files/AURORA_retrospective_Supplementary_Table.2.xlsx", sheet = 6) %>% 
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
                             `Anatomic Site Groups`, "Other"))

geicam_expr <- read.table("/datastore/lbcfs/projects/nextflow/peroulab/aurora_123_rap96_geicam_204_harmony_186_sg/rna-seq-quant-2.5/hg38/gencode_v36/merge/summary/vqg3nccallvj8jbsz/salmon_ensemblgid.matrix.tsv",
                          sep = "\t", check.names = F, header = T, row.names = 1)
geicam_expr <- geicam_expr[,grepl("AP_", colnames(geicam_expr))]
colnames(geicam_expr) <- sub(".*[YCS]_", "", colnames(geicam_expr))

# Remove gene variants: 
ids <- gsub("\\..*$", "", rownames(geicam_expr))
expr <- rowsum(geicam_expr, group = ids)

# Upper-quartile normalize
uqNormalize <- function(expr) {
  uqs <- apply(expr, 2, function(x){ quantile(x[x > 0 & !is.na(x)], 0.75) })
  y <- 1000
  expr <- apply(expr, 1, function(x, y){ x * y }, y / uqs)
  expr <- t(expr)
  return(expr)
}

geicam_expr <- uqNormalize(expr)

set.seed(123)

## Calculate module scores -----------------------------------------------------
modules <- calculate_modules(geicam_expr, 
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

write.table(modules, paste0(out_dir, "GEICAM_n204_Ensembl_938_published_modules.txt"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T)

## Add ER, PR, HER2 (IHC) status as features (0/1) -----------------------------
receptor_status <- read.table("output/apply_EN_models/process_test_sets/GEICAM/infer_ER_PR_HER2/GEICAM_n204_Ensembl_inferred_clin_groups.txt", 
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
  select(jUNCtion_ID, ER_status_w_inferred, PR_status_w_inferred, HER2_status_w_inferred) %>% 
  rename(
    clinical_ER_status = ER_status_w_inferred, 
    clinical_PR_status = PR_status_w_inferred, 
    clinical_HER2_status = HER2_status_w_inferred)

modules <- as.data.frame(t(modules)) %>% 
  tibble::rownames_to_column(var = "jUNCtion_ID")
modules_receptors <- left_join(modules, receptor_status, by="jUNCtion_ID")
modules_receptors <- tibble::column_to_rownames(modules_receptors, var = "jUNCtion_ID")
modules_receptors <- as.data.frame(t(modules_receptors))

write.table(modules_receptors, paste0(out_dir, "GEICAM_n204_Ensembl_938_published_modules_3_receptors.txt"), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T)
