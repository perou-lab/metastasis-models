# Calculate distant metastasis EN scores on UNC337 dataset

library(tidyverse)
library(ggplot2)
library(matrixStats)
library(ggpubr)
library(survminer)
library(survival)

out_dir <- "output/apply_EN_models/UNC337/genes_10814/distant/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

pam50_colors <- c("Basal" = "#c91508", 
                  "Her2" = "#ed8ccb",
                  "LumA" = "#084dc4",
                  "LumB" = "#9cd0f0",
                  "Normal" = "#4fb30c", 
                  "Claudin" = "goldenrod")

## Apply distant metastasis model -----
metadata <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_metadata.tsv", 
                       sep = "\t", header = T, check.names = F)
unc_expr <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_primaries_combined_ensembl_id.tsv", 
                       sep = "\t", header = T, check.names = F)

distant_model <- readRDS("output/train_EN_models/n2486/genes_10814/NoRel-v-Distant/EN_75Train_25Test_dwd-genes_glmnet_obj_s124022_NoRel-v-Distant_boot_accuracy.rds")
best_lambda <- distant_model$bestTune$lambda

# extract model with weighted coefficients
coef_vec <- as.matrix(coef(distant_model$finalModel, s = best_lambda))
nonzero_idx <- which(coef_vec != 0) # non zero model coefficients
coef_vec <- coef_vec[nonzero_idx, , drop = FALSE]

intercept <- coef_vec["(Intercept)", 1] # model intercept
coef_vec <- coef_vec[rownames(coef_vec) != "(Intercept)", , drop = FALSE]

## align test set gene list to model list
common_genes <- intersect(rownames(coef_vec), rownames(unc_expr))
x <- t(unc_expr[common_genes, ])
beta <- coef_vec[common_genes, 1]

## compute model scores
linear_pred <- x %*% beta + intercept
test_prediction <- 1 / (1 + exp(-linear_pred))
colnames(test_prediction)[1] <- "Distant_EN_score"

test_prediction <- test_prediction %>% 
  as.data.frame() %>% 
  # 0 ~ No_Relapse, 1 ~ Relapse; threshold = 0.5
  tibble::rownames_to_column(var = "Sample_GEO_accession")

hist(test_prediction$Distant_EN_score, breaks = 50)

metadata <- metadata %>% 
  left_join(test_prediction, by = "Sample_GEO_accession")

write.table(test_prediction, paste0(out_dir, "UNC337_distant_gene_EN_model_scores.tsv"), 
            sep = "\t", row.names = F, col.names = T)
