# Calculate luminal bone metastasis EN scores on UNC337 dataset

library(tidyverse)

out_dir <- "output/apply_EN_models/UNC337/modules_935+3/lum_bone/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Apply luminal bone metastasis model -----
metadata <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_metadata.tsv", 
                       sep = "\t", header = T, check.names = F)
unc_expr <- read.table("output/apply_EN_models/process_test_sets/UNC337/calculate_modules/UNC337_n239_primaries_938_published_modules_3_receptors.txt", 
                       sep = "\t", header = T, check.names = F)

lum_bone_model <- readRDS("output/train_EN_models/lumA_lumB_only/modules_935+3/NoRel-v-OnlyBone/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-OnlyBone_boot_accuracy.rds")
best_lambda <- lum_bone_model$bestTune$lambda

## Extract model features
model_features <- rownames(coef(lum_bone_model$finalModel, s = best_lambda))
model_features <- setdiff(model_features, "(Intercept)")

missing_features <- setdiff(model_features, rownames(unc_expr))

newx <- t(as.matrix(unc_expr[model_features, , drop = FALSE]))

## Make predictions
unc_pred <- predict(lum_bone_model$finalModel,
                    newx = newx,
                    s = best_lambda,
                    type = "response")

unc_predictions <- data.frame(Sample_GEO_accession = rownames(unc_pred),
                              Luminal_Bone_EN_score = as.numeric(unc_pred))

metadata <- metadata %>% 
  left_join(unc_predictions, by = "Sample_GEO_accession")

write.table(unc_predictions, paste0(out_dir, "UNC337_luminal_bone_module_EN_model_scores.tsv"), 
            sep = "\t", row.names = F, col.names = T)
