# Calculate multi-site metastasis EN scores on SCANB dataset

library(tidyverse)
library(ggplot2)
library(matrixStats)
library(ggpubr)
library(survminer)
library(survival)

out_dir <- "output/apply_EN_models/SCANB/modules_935+3/multimet/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

pam50_colors <- c("Basal" = "#c91508", 
                  "Her2" = "#ed8ccb",
                  "LumA" = "#084dc4",
                  "LumB" = "#9cd0f0",
                  "Normal" = "#4fb30c", 
                  "Clow" = "goldenrod")

## Apply multi-site metastasis model -----
metadata <- read.table("../datasets/SCANB_gitlab/SCAN-B_6329_Clinical_Data.txt", 
                       sep = "\t", header = T, check.names = F)
scanb_expr <- read.table("output/apply_EN_models/process_test_sets/SCANB/calculate_modules/SCANB_Ensembl_938_published_modules_3_receptors.txt", 
                  sep = "\t", header = T, check.names = F)

mm_model <- readRDS("output/train_EN_models/n2486/modules_935+3/NoRel-v-MultiMet/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-MultiMet_boot_accuracy.rds")
best_lambda <- mm_model$bestTune$lambda

# extract model with weighted coefficients
coef_vec <- as.matrix(coef(mm_model$finalModel, s = best_lambda))
nonzero_idx <- which(coef_vec != 0) # non zero model coefficients
coef_vec <- coef_vec[nonzero_idx, , drop = FALSE]

intercept <- coef_vec["(Intercept)", 1] # model intercept
coef_vec <- coef_vec[rownames(coef_vec) != "(Intercept)", , drop = FALSE]

## align test set gene list to model list
common_genes <- intersect(rownames(coef_vec), rownames(scanb_expr))
x <- t(scanb_expr[common_genes, ])
beta <- coef_vec[common_genes, 1]

## compute model scores
linear_pred <- x %*% beta + intercept
test_prediction <- 1 / (1 + exp(-linear_pred))
colnames(test_prediction)[1] <- "MultiMet_EN_Score"

test_prediction <- test_prediction %>% 
  as.data.frame() %>% 
  # 0 ~ No_Relapse, 1 ~ Relapse; threshold = 0.5
  tibble::rownames_to_column(var = "Merged")

hist(test_prediction$MultiMet_EN_Score, breaks = 50)

metadata <- metadata %>% 
  left_join(test_prediction, by = "Merged")

write.table(test_prediction, paste0(out_dir, "SCANB_multimet_module_EN_model_scores.tsv"), 
            sep = "\t", row.names = F, col.names = T)
