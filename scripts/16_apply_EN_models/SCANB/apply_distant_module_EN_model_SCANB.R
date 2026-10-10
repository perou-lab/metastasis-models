# Calculate distant metastasis EN scores on SCANB dataset

library(tidyverse)

out_dir <- "output/apply_EN_models/SCANB/modules_935+3/distant/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Apply distant metastasis model -----
metadata <- read.table("../datasets/SCANB_gitlab/SCAN-B_6329_Clinical_Data.txt", 
                       sep = "\t", header = T, check.names = F)
scanb_expr <- read.table("output/apply_EN_models/process_test_sets/SCANB/calculate_modules/SCANB_Ensembl_938_published_modules_3_receptors.txt", 
                         sep = "\t", header = T, check.names = F)

distant_model <- readRDS("output/train_EN_models/n2486/modules_935+3/NoRel-v-Distant/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-Distant_boot_accuracy.rds")
best_lambda <- distant_model$bestTune$lambda

## Extract model features
model_features <- rownames(coef(distant_model$finalModel, s = best_lambda))
model_features <- setdiff(model_features, "(Intercept)")

missing_features <- setdiff(model_features, rownames(scanb_expr))

newx <- t(as.matrix(scanb_expr[model_features, , drop = FALSE]))

## Make predictions
scanb_pred <- predict(distant_model$finalModel,
                      newx = newx,
                      s = best_lambda,
                      type = "response")

scanb_predictions <- data.frame(Merged = colnames(scanb_expr),
                                Distant_EN_Score = as.numeric(scanb_pred))

metadata <- metadata %>% 
  left_join(scanb_predictions, by = "Merged")

write.table(metadata, paste0(out_dir, "SCANB_distant_module_EN_model_scores.tsv"),
            sep = "\t", row.names = F, col.names = T)
