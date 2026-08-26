## Compare bone AUC on prognostic models + luminal bone met EN model

library(tidyverse)
library(ggplot2)
library(matrixStats)
library(ggpubr)
library(survminer)
library(survival)

out_dir <- "output/compare_models_on_bone_mets/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Apply luminal bone metastasis model to all n2486 samples --------------------

metadata <- read.table("output/balanced_stratification/n2486/n2486_complete_meta_w_onlysite_stratification.txt", 
                       sep = "\t", header = T, check.names = F)
modules <- read.table("output/calculate_modules/n2486_Ensembl_dwd_935_published_modules_3_receptors.txt", 
                      sep = "\t", header = T, check.names = F)

lum_bone_model <- readRDS("output/train_EN_models/lumA_lumB_only/modules_935+3/NoRel-v-OnlyBone/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-OnlyBone_boot_accuracy.rds")
best_lambda <- lum_bone_model$bestTune$lambda

# extract model with weighted coefficients
coef_vec <- as.matrix(coef(lum_bone_model$finalModel, s = best_lambda))
nonzero_idx <- which(coef_vec != 0) # non zero model coefficients
coef_vec <- coef_vec[nonzero_idx, , drop = FALSE]

intercept <- coef_vec["(Intercept)", 1] # model intercept
coef_vec <- coef_vec[rownames(coef_vec) != "(Intercept)", , drop = FALSE]

## align test set gene list to model list
common_genes <- intersect(rownames(coef_vec), rownames(modules))
x <- t(modules[common_genes, ])
beta <- coef_vec[common_genes, 1]

## compute model scores
linear_pred <- x %*% beta + intercept
test_prediction <- 1 / (1 + exp(-linear_pred))
colnames(test_prediction)[1] <- "Luminal_Bone_EN_score"

test_prediction <- test_prediction %>% 
  as.data.frame() %>% 
  # 0 ~ No_Relapse, 1 ~ Relapse; threshold = 0.5
  tibble::rownames_to_column(var = "Sample_ID")

hist(test_prediction$Luminal_Bone_EN_score, breaks = 50)

metadata <- metadata %>% 
  left_join(test_prediction, by = "Sample_ID")

write.table(test_prediction, paste0(out_dir, "n2486_luminal_bone_module_EN_model_scores.tsv"), 
            sep = "\t", row.names = F, col.names = T)

## Identify test set samples: HR+/HER2 samples not included in luminal training set
luminal_meta <- read.table("output/balanced_stratification/lumA_lumB_only/lumA_lumB_only_complete_meta_w_onlysite_stratification.txt",
                           sep = "\t", header = T, check.names = F)

test_samples <- luminal_meta %>%
  filter(Train75_Test25 == "Test") %>%
  filter(Clinical_Group_3 == "HR+/HER2-") %>% 
  pull(Sample_ID)

## Read in modules and combine with metadata -----------------------------------
modules <- read.table("output/calculate_modules/n2486_Ensembl_dwd_938_published_modules.txt",
                      sep = "\t", header = T)

modules <- modules[c("ROR_P_Model_JCO.2009_PMID.19204204",
                     "GHI_RS_Model_NEJM.2004_PMID.15591335",
                     "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860"),
                   test_samples]

modules <- as.data.frame(t(modules)) %>% 
  tibble::rownames_to_column("Sample_ID")

lum_meta <- metadata %>% 
  filter(Sample_ID %in% test_samples) %>% 
  left_join(modules, by = "Sample_ID")

## Calculate AUC values --------------------------------------------------------

auc_results <- data.frame(
  signature = character(),
  group = character(),
  auc = numeric(),
  stringsAsFactors = FALSE
)

signatures <- c("ROR_P_Model_JCO.2009_PMID.19204204",
                "GHI_RS_Model_NEJM.2004_PMID.15591335",
                "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860", 
                "Luminal_Bone_EN_score")

for (signature in signatures) {
  group <- "Bone"
  sub_meta <- lum_meta %>% 
    filter(Relapse_Site %in% c("No_Relapse", group))
  
  # Check if the signature exists
  if (!signature %in% colnames(sub_meta)) {
    warning(paste("Skipping module", signature, "in group", group,
                  "signature column not found in test_meta"))
    next
  }
  
  # Remove samples with missing score or outcome
  sub_meta <- sub_meta %>% 
    filter(!is.na(.data[[signature]]), 
           !is.na(No_v_Distant))
  
  y <- sub_meta$No_v_Distant
  scores <- sub_meta[[signature]]
  
  # Only calculate AUC if column exists and y has 2 classes
  if (length(unique(y)) == 2) {
    if(signature == "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860") {
      pred <- prediction(
        predictions = scores,
        labels = y,
        label.ordering = c("Distant_Relapse", "No_Relapse")
    ) } else {
      pred <- prediction(
        predictions = scores,
        labels = y,
        label.ordering = c("No_Relapse", "Distant_Relapse")
    ) }
    
    auc <- signif(performance(pred, measure = "auc")@y.values[[1]][1], 4)
    
    # Store result
    auc_results <- rbind(
      auc_results,
      data.frame(signature = signature, group = group, auc = auc)
    )
  } else {
    warning(paste("Skipping module", signature, "- missing or invalid input"))
  }
}

auc_results <- auc_results %>% 
  mutate(short_name = case_when(
  signature == "Luminal_Bone_EN_score" ~ "Lum_Bone_EN_Model",
  signature == "ROR_P_Model_JCO.2009_PMID.19204204" ~ "ROR_P_Model", 
  signature == "GHI_RS_Model_NEJM.2004_PMID.15591335" ~ "GHI_RS_Model", 
  signature == "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860" ~ "Pcorr_NKI70")) %>% 
  mutate(short_name = factor(short_name,
                            levels = c("ROR_P_Model",
                                       "GHI_RS_Model",
                                       "Pcorr_NKI70",
                                       "Lum_Bone_EN_Model")))
  
p <- ggplot(auc_results, aes(x = short_name, y = auc)) +
  geom_bar(stat = "identity") + 
  geom_col(color = "black", fill = "gray80", linewidth = 1) +
  geom_text(aes(label = round(auc,3)), vjust = -0.3, size = 4) + 
  theme_minimal() + 
  labs(x = "Group",
       y = "AUC", 
       title = paste0("Comparison of bone met models")) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(limits = c(0, 1))
ggsave(filename = paste0(out_dir, "Bone_only_performance_comparison.png"),
       plot = p, width = 4, height = 4, units = "in", dpi = 350)

pdf(paste0(out_dir, "Bone_only_performance_comparison.pdf"),
    width = 4, height = 4)
print(p)
dev.off()

