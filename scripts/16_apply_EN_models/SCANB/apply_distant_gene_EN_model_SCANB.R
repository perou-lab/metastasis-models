# Calculate distant metastasis EN scores on SCANB dataset

library(tidyverse)
library(matrixStats)

out_dir <- "output/apply_EN_models/SCANB/genes_10814/distant/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Read in distant metastasis model (genes) -----
distant_model <- readRDS("output/train_EN_models/n2486/genes_10814/NoRel-v-Distant/EN_75Train_25Test_dwd-genes_glmnet_obj_s124022_NoRel-v-Distant_boot_accuracy.rds")
best_lambda <- distant_model$bestTune$lambda

## Extract model features
model_features <- rownames(coef(distant_model$finalModel, s = best_lambda))
model_features <- setdiff(model_features, "(Intercept)")

## Apply distant metastasis model -----
metadata <- read.table("../datasets/SCANB_gitlab/SCAN-B_6329_Clinical_Data.txt", 
                       sep = "\t", header = T, check.names = F)
scanb_expr <- readRDS("../metastatic_site_EN_manuscript/source_files/SCANB_MDS/human_scanb_ensembl_norm.RDS")
scanb_mat <- assay(scanb_expr)
scanb_df <- as.data.frame(scanb_mat)

rm(scanb_expr, scanb_mat)

## Processing SCAN-B data (log2 transform, median center, column standardize)
scanb_df <- scanb_df[model_features, , drop = FALSE]
scanb <- as.matrix(log2(scanb_df + 1))
scanb <- scanb - rowMedians(scanb)
scanb <- scale(scanb)

newx <- t(scanb)

## Make predictions
scanb_pred <- predict(distant_model$finalModel,
                      newx = newx,
                      s = best_lambda,
                      type = "response")

scanb_predictions <- data.frame(Merged = rownames(newx),
                                Distant_EN_Score = as.numeric(scanb_pred))

metadata <- metadata %>% 
  left_join(scanb_predictions, by = "Merged")

write.table(metadata, paste0(out_dir, "SCANB_distant_gene_EN_model_scores.tsv"),
            sep = "\t", row.names = F, col.names = T)
