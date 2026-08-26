## Test performance of research version of commercial predictors on HR+/HER2- 
## tumors with known met sites
# 2026.08.12

library(tidyverse)
library(readxl)
library(ROCR)
library(ggplot2)

out_dir <- "output/compare_prognostic_assays/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Read in all metadata and modules --------------------------------------------
meta <- read.table("output/balanced_stratification/n2486/n2486_complete_meta_w_onlysite_stratification.txt",
                   sep = "\t", header = T, check.names = F)

meta <- meta %>% 
  mutate(No_v_Distant = factor(meta$No_v_Distant, levels = c("No_Relapse", "Distant_Relapse"))) %>% 
  filter(Clinical_Group_3 == "HR+/HER2-")

modules <- read.table("output/calculate_modules/n2486_Ensembl_dwd_938_published_modules.txt",
                      sep = "\t", header = T)

## Calculate AUC values for each group based on site of metastasis (brain, bone, liver, lung, multi-site)
auc_results <- data.frame(
  signature = character(),
  group = character(),
  auc = numeric(),
  stringsAsFactors = FALSE
)

sub_groups <- c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet")

prognostic_signatures <- c("ROR_P_Model_JCO.2009_PMID.19204204",
                           "GHI_RS_Model_NEJM.2004_PMID.15591335",
                           "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860")

for (signature in prognostic_signatures) {
  for (group in sub_groups) {
    
    # Subset metadata depending on group type
    if (group %in% c("Bone", "Brain", "Lung", "Liver", "Multimet")) {
      sub_meta <- meta %>% 
        filter(Relapse_Site %in% c("No_Relapse", group))
    } else {
      sub_meta <- meta
    }
    
    # Subset module columns based on sample IDs in sub_meta
    sample_ids <- intersect(sub_meta$Sample_ID, colnames(modules))
    
    sub_meta <- sub_meta %>%
      filter(Sample_ID %in% sample_ids) %>%
      slice(match(sample_ids, Sample_ID))
    
    sub_modules <- modules[, sample_ids, drop = FALSE]

    # Transpose for prediction
    x <- as.matrix(t(as.matrix(sub_modules)))
    y <- sub_meta$No_v_Distant
    
    # Only calculate AUC if column exists and y has 2 classes
    if (signature %in% colnames(x) && length(unique(y)) == 2) {
      if (signature == "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860") {
      pred <- prediction(
        predictions = x[, signature],
        labels = y,
        label.ordering = c("Distant_Relapse", "No_Relapse")
      ) } else {
        pred <- prediction(
          predictions = x[, signature],
          labels = y,
          label.ordering = c("No_Relapse", "Distant_Relapse")
        )   
      }
      
      auc <- signif(performance(pred, measure = "auc")@y.values[[1]][1], 4)
      
      # Store result
      auc_results <- rbind(
        auc_results,
        data.frame(signature = signature, group = group, auc = auc)
      )
    } else {
      warning(paste("Skipping ", signature, "in group", group, "- missing or invalid input"))
    }
  }
}


for (signature_name in unique(auc_results$signature)) {
  # Filter the results for this module
  signature_data <- auc_results %>%
    filter(signature == signature_name)
  
  short_name <- case_when(
    signature_name == "ROR_P_Model_JCO.2009_PMID.19204204" ~ "ROR_P_Model", 
    signature_name == "GHI_RS_Model_NEJM.2004_PMID.15591335" ~ "GHI_RS_Model", 
    signature_name == "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860" ~ "Pcorr_NKI70"
  )
  
  signature_data$group <- factor(signature_data$group, levels = unique(signature_data$group))
  
p <- ggplot(signature_data, aes(x = group, y = auc)) +
  geom_bar(stat = "identity") + 
  geom_col(color = "black", fill = "gray80", linewidth = 1) +
  geom_text(aes(label = round(auc,3)), vjust = -0.3, size = 4) + 
  theme_minimal() + 
  labs(x = "Group",
       y = "AUC", 
       title = paste0(signature_name)) +
#  theme(axis.text.x = element_text(hjust = 1)) + 
  scale_y_continuous(limits = c(0, 1))
ggsave(filename = paste0(out_dir, short_name, "_performance_comparison_HR+HER2-.png"),
       plot = p, width = 4, height = 4, units = "in", dpi = 350)

pdf(paste0(out_dir, short_name, "_performance_comparison_HR+HER2-.pdf"),
    width = 4, height = 4)
print(p)
dev.off()
}
