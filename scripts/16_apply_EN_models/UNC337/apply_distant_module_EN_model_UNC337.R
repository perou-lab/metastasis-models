# Calculate multi-site metastasis EN scores on UNC337 dataset

library(tidyverse)

out_dir <- "output/apply_EN_models/UNC337/modules_935+3/distant/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Apply multi-site metastasis model -----
metadata <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_metadata.tsv", 
                       sep = "\t", header = T, check.names = F)
unc_expr <- read.table("output/apply_EN_models/process_test_sets/UNC337/calculate_modules/UNC337_n239_primaries_938_published_modules_3_receptors.txt", 
                       sep = "\t", header = T, check.names = F)

distant_model <- readRDS("output/train_EN_models/n2486/modules_935+3/NoRel-v-Distant/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-Distant_boot_accuracy.rds")
best_lambda <- distant_model$bestTune$lambda

## Extract model features
model_features <- rownames(coef(distant_model$finalModel, s = best_lambda))
model_features <- setdiff(model_features, "(Intercept)")

missing_features <- setdiff(model_features, rownames(unc_expr))

newx <- t(as.matrix(unc_expr[model_features, , drop = FALSE]))

## Make predictions
unc_pred <- predict(distant_model$finalModel,
                    newx = newx,
                    s = best_lambda,
                    type = "response")

unc_predictions <- data.frame(Sample_GEO_accession = colnames(unc_expr),
                              Distant_EN_Score = as.numeric(unc_pred))

metadata <- metadata %>% 
  left_join(unc_predictions, by = "Sample_GEO_accession")

write.table(unc_predictions, paste0(out_dir, "UNC337_distant_module_EN_model_scores.tsv"), 
            sep = "\t", row.names = F, col.names = T)
