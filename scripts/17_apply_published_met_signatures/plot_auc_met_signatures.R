# Plot AUC curves for published metastasis signatures according to site specific 
# metastasis behavior

library(tidyverse)
library(ROCR)

out_dir <- "output/apply_published_met_signatures/AUC/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Read in metadata and modules ------------------------------------------------
metadata <- read.table("output/balanced_stratification/n2486/n2486_complete_meta_w_onlysite_stratification.txt",
                       sep = "\t", header = T) %>% 
  mutate(No_v_Distant = factor(No_v_Distant, levels = c("No_Relapse", "Distant_Relapse")))

signatures <- read.table("output/apply_published_met_signatures/published_met_scores_applied_to_n2486.tsv",
                      sep = "\t", header = T) %>% 
  tibble::column_to_rownames(var = "Sample_ID")

meta <- metadata[match(rownames(signatures), metadata$Sample_ID), ]

auc_results <- data.frame(
  signature = character(), 
  group = character(), 
  auc = numeric(), 
  stringsAsFactors = FALSE
)

sub_groups <- c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet",
                "Multimet w Bone", "Multimet wo Bone",
                "Multimet w Brain", "Multimet wo Brain", 
                "Multimet w Lung", "Multimet wo Lung")

met_signatures <- colnames(signatures)

meta <- meta %>% 
  mutate(multimet_with_bone = case_when(
    No_v_MultiMet == "Multimet_Relapse" & str_detect(Four_site_recurrence, "Bone") ~ "MM_w_bone", 
    No_v_MultiMet == "Multimet_Relapse" & !str_detect(Four_site_recurrence, "Bone") ~ "MM_wo_bone", 
    No_v_MultiMet == "No_Relapse" ~ "No_Relapse")) %>% 
  mutate(multimet_with_lung = case_when(
    No_v_MultiMet == "Multimet_Relapse" & str_detect(Four_site_recurrence, "Lung") ~ "MM_w_lung", 
    No_v_MultiMet == "Multimet_Relapse" & !str_detect(Four_site_recurrence, "Lung") ~ "MM_wo_lung", 
    No_v_MultiMet == "No_Relapse" ~ "No_Relapse")) %>% 
  mutate(multimet_with_brain = case_when(
    No_v_MultiMet == "Multimet_Relapse" & str_detect(Four_site_recurrence, "Brain") ~ "MM_w_brain", 
    No_v_MultiMet == "Multimet_Relapse" & !str_detect(Four_site_recurrence, "Brain") ~ "MM_wo_brain", 
    No_v_MultiMet == "No_Relapse" ~ "No_Relapse")) %>% 
  mutate(multimet_with_liver = case_when(
    No_v_MultiMet == "Multimet_Relapse" & str_detect(Four_site_recurrence, "Liver") ~ "MM_w_Liver", 
    No_v_MultiMet == "Multimet_Relapse" & !str_detect(Four_site_recurrence, "Liver") ~ "MM_wo_Liver", 
    No_v_MultiMet == "No_Relapse" ~ "No_Relapse"))

for (signature in met_signatures) {
  for(group in sub_groups) {
    
    # Subset metadata depending on group type
    if (group %in% c("Bone", "Brain", "Liver", "Lung", "Multimet")) {
      sub_meta <- meta %>% 
        filter(Relapse_Site %in% c("No_Relapse", group))
    } else if (group == "Multimet w Bone") { 
      sub_meta <- meta %>% 
        filter(multimet_with_bone %in% c("No_Relapse", "MM_w_bone")) 
    } else if (group == "Multimet wo Bone") {
      sub_meta <- meta %>% 
        filter(multimet_with_bone %in% c("No_Relapse", "MM_wo_bone")) 
    } else if (group == "Multimet w Lung") { 
      sub_meta <- meta %>% 
        filter(multimet_with_lung %in% c("No_Relapse", "MM_w_lung")) 
    } else if (group == "Multimet wo Lung") {
      sub_meta <- meta %>% 
        filter(multimet_with_lung %in% c("No_Relapse", "MM_wo_lung")) 
    } else if (group == "Multimet w Brain") { 
      sub_meta <- meta %>% 
        filter(multimet_with_brain %in% c("No_Relapse", "MM_w_brain")) 
    } else if (group == "Multimet wo Brain") {
      sub_meta <- meta %>% 
        filter(multimet_with_brain %in% c("No_Relapse", "MM_wo_brain")) 
    } else (sub_meta <- meta)
    
    # Subset signature rows based on sample IDs in sub_meta
    sub_signatures <- signatures[rownames(signatures) %in% sub_meta$Sample_ID,]
    
    # Set up for prediction
    x <- as.matrix(sub_signatures)
    y <- sub_meta$No_v_Distant
    
    # Only calculate AUC if column exists and y has 2 classes 
    if (signature %in% colnames(x) && length(unique(y)) == 2) {
      pred <- prediction(
        predictions = x[, signature],
        labels = y,
        label.ordering = c("No_Relapse", "Distant_Relapse")
      )
      
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
  # Filter results for a specific signature
  signature_data <- auc_results %>% 
    filter(signature == signature_name)
  
  short_name <- case_when(
    signature_name == "Lung_classifier_Minn" ~ "Lung_classifier_Minn", 
    signature_name == "Brain_12_gene_Score_Bos_2009" ~ "Brain_12_up_genes", 
    signature_name == "Lung_score_Landemaine_2008" ~ "Lung_Landemaine", 
    signature_name == "Zhao_Bone_score" ~ "Bone_Score_Zhao"
  )
  
  # Subset metadata depending on group type
  if (short_name %in% c("Bone_Score_Zhao")) {
    plot_data <- signature_data %>% 
      filter(group %in% c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet", "Multimet w Bone", "Multimet wo Bone")) %>% 
      mutate(group = factor(group, levels = c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet", "Multimet w Bone", "Multimet wo Bone")))
    
  } else if (short_name %in% c("Lung_classifier_Minn", "Lung_Landemaine")) { 
    plot_data <- signature_data %>% 
      filter(group %in% c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet", "Multimet w Lung", "Multimet wo Lung")) %>% 
      mutate(group = factor(group, levels = c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet", "Multimet w Lung", "Multimet wo Lung")))
    
  } else if (short_name == "Brain_12_up_genes") {
    plot_data <- signature_data %>% 
      filter(group %in% c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet", "Multimet w Brain", "Multimet wo Brain")) %>% 
      mutate(group = factor(group, levels = c("n2486", "Bone", "Brain", "Lung", "Liver", "Multimet", "Multimet w Brain", "Multimet wo Brain")))
    
  } else (plot_data <- signature_data)
  
  p <- ggplot(plot_data, aes(x = group, y = auc)) +
    geom_bar(stat = "identity") + 
    geom_col(color = "black", fill = "gray80", linewidth = 1) +
    geom_text(aes(label = round(auc,3)), vjust = -0.3, size = 4) + 
    theme_minimal() + 
    labs(x = "Group",
         y = "AUC", 
         title = paste0(signature_name)) +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(limits = c(0, 1))
  ggsave(filename = paste0(out_dir, short_name, "_performance_comparison.png"),
         plot = p, width = 6, height = 4.5, units = "in", dpi = 350)
  
  ggsave(filename = paste0(out_dir, short_name, "_performance_comparison.pdf"),
         plot = p, width = 6, height = 4.5, units = "in")
}

