# Process raw SCAN-B dataset to apply EN models

library(tidyverse)

out_dir <- "output/apply_EN_models/process_test_sets/SCANB/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Read in SCAN-B dataset (ensembl ids) ----------------------------------------
metadata <- read.table("../datasets/SCANB_gitlab/SCAN-B_6329_Clinical_Data.txt", 
                       sep = "\t", header = T, check.names = F)

ensembl_expr <- readRDS("../metastatic_site_EN_manuscript/source_files/SCANB_MDS/human_scanb_ensembl_norm.RDS")
expr_mat <- assay(ensembl_expr)
expr_df <- as.data.frame(expr_mat)
rm(ensembl_expr, expr_mat)


## Read in EN models -----------------------------------------------------------
# Distant metastasis (genes) #
distant_model <- readRDS("output/train_EN_models/n2486/genes_10814/NoRel-v-Distant/EN_75Train_25Test_dwd-genes_glmnet_obj_s124022_NoRel-v-Distant_boot_accuracy.rds")
best_lambda <- distant_model$bestTune$lambda
coefs <- as.matrix(coef(distant_model$finalModel, s = best_lambda))
distant_model_genes <- rownames(coefs)[which(coefs!= 0)]
distant_model_genes <- setdiff(distant_model_genes, "(Intercept)")

# Luminal Bone metastasis (genes) #
lum_bone_model <- readRDS("output/train_EN_models/lumA_lumB_only/genes_10814/NoRel-v-OnlyBone/EN_75Train_25Test_dwd-genes_glmnet_obj_s124022_NoRel-v-OnlyBone_boot_accuracy.rds")
best_lambda <- lum_bone_model$bestTune$lambda
coefs <- as.matrix(coef(lum_bone_model$finalModel, s = best_lambda))
lum_bone_model_genes <- rownames(coefs)[which(coefs != 0)]
lum_bone_model_genes <- setdiff(lum_bone_model_genes, "(Intercept)")

# Multi-site metastasis (genes) #
mm_model <- readRDS("output/train_EN_models/n2486/genes_10814/NoRel-v-MultiMet/EN_75Train_25Test_dwd-genes_glmnet_obj_s124022_NoRel-v-MultiMet_boot_accuracy.rds")
best_lambda <- mm_model$bestTune$lambda
coefs <- as.matrix(coef(mm_model$finalModel, s = best_lambda))
mm_model_genes <- rownames(coefs)[which(coefs != 0)]
mm_model_genes <- setdiff(mm_model_genes, "(Intercept)")


## Subset SCAN-B expression by gene lists --------------------------------------
distant_expr <- expr_df[rownames(expr_df) %in% distant_model_genes,]
distant_expr <- log2(distant_expr + 1)
write.table(distant_expr, paste0(out_dir, "SCANB_No-v-Distant_genes_n144.tsv"), 
            sep = "\t", row.names = T, col.names = T)

lum_bone_expr <- expr_df[rownames(expr_df) %in% lum_bone_model_genes,]
lum_bone_expr <- log2(lum_bone_expr + 1)
write.table(lum_bone_expr, paste0(out_dir, "SCANB_LumAB_No-v-BoneOnly_genes_n511.tsv"), 
            sep = "\t", row.names = T, col.names = T)

mm_expr <- expr_df[rownames(expr_df) %in% mm_model_genes,]
mm_expr <- log2(mm_expr + 1)
write.table(mm_expr, paste0(out_dir, "SCANB_No-v-MultiMet_genes_n132.tsv"), 
            sep = "\t", row.names = T, col.names = T)
