# Apply luminal bone metastasis EN model to GEICAM samples

library(tidyverse)
library(ggplot2)
library(matrixStats)
library(ggpubr)
library(survminer)
library(survival)
library(readxl)

out_dir <- "output/apply_EN_models/GEICAM/modules_935+3/lum_bone/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

pam50_colors <- c("Basal" = "#c91508", 
                  "Her2" = "#ed8ccb",
                  "LumA" = "#084dc4",
                  "LumB" = "#9cd0f0",
                  "Normal" = "#4fb30c", 
                  "Claudin" = "goldenrod")

metastasis_site_colors <- c("Bone" = "#a18be8", "Brain" = "#a11045", "Liver" = "#a9c451", "Lung" = "#cc8014",
                            "Breast" = "deeppink3", "Endometrial" = "#CDAD00", "Lymph node" = "#1492cc", "Ovarian" = "slategray3", 
                            "Pericardium" = "#CD5555", "Peritoneum" = "lightgray",  "Pleura" = "#96CDCD", "Skin" = "darkblue")

## Apply distant metastasis model -----
## Read in GEICAM metadata
geicam_meta <- read_excel("source_files/AURORA_retrospective_Supplementary_Table.2.xlsx", sheet = 6) %>% 
  mutate(jUNCtion_ID = sub("_[A-Z0-9]+$", "", sub("-.*", "", `Sample Name`))) %>% 
  rename(Sample = `Short name`) %>% 
  rename(Set = Study) %>% 
  rename(Tissue_type = `Tissue type`) %>% 
  rename(Sample_type = PvsM) %>% 
  rename(PAM50_call = `PAM50 Call`) %>%
  filter(Tissue_type == "FFPE") %>% 
  rename(Site_original = `Anatomic Site`) %>% 
  mutate(Site_group = ifelse(`Anatomic Site Groups` %in% 
                               c("Brain", "Breast", "Liver",
                                 "Lung", "Lymph node"),
                             `Anatomic Site Groups`, "Other"))

geicam_expr <- read.table("output/apply_EN_models/process_test_sets/GEICAM/calculate_modules/GEICAM_n204_Ensembl_938_published_modules_3_receptors.txt",
                          sep = "\t", check.names = F, header = T, row.names = 1)

# Luminal bone EN model
lum_bone_model <- readRDS("output/train_EN_models/lumA_lumB_only/modules_935+3/NoRel-v-OnlyBone/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-OnlyBone_boot_accuracy.rds")
best_lambda <- lum_bone_model$bestTune$lambda

# extract model with weighted coefficients
coef_vec <- as.matrix(coef(lum_bone_model$finalModel, s = best_lambda))
nonzero_idx <- which(coef_vec != 0) # non zero model coefficients
coef_vec <- coef_vec[nonzero_idx, , drop = FALSE]

intercept <- coef_vec["(Intercept)", 1] # model intercept
coef_vec <- coef_vec[rownames(coef_vec) != "(Intercept)", , drop = FALSE]

## align test set gene list to model list
common_genes <- intersect(rownames(coef_vec), rownames(geicam_expr))
x <- t(geicam_expr[common_genes, ])
beta <- coef_vec[common_genes, 1]

## compute model scores
linear_pred <- x %*% beta + intercept
test_prediction <- 1 / (1 + exp(-linear_pred))
colnames(test_prediction)[1] <- "Luminal_Bone_EN_score"

test_prediction <- test_prediction %>% 
  as.data.frame() %>% 
  # 0 ~ No_Relapse, 1 ~ Relapse; threshold = 0.5
  tibble::rownames_to_column(var = "jUNCtion_ID")

hist(test_prediction$Luminal_Bone_EN_score, breaks = 50)

metadata <- geicam_meta %>% 
  left_join(test_prediction, by = "jUNCtion_ID")

write.table(test_prediction, paste0(out_dir, "GEICAM_luminal_bone_module_EN_model_scores.tsv"), 
            sep = "\t", row.names = F, col.names = T)

## Luminal patients only -------------------------------------------------------

## Plot signature expression between primaries with bone mets and those without bone mets 
metadata <- metadata %>% 
  mutate(Patient_ID = gsub(".*-(\\d+)$", "\\1", Sample))

bone_patients <- metadata %>% 
  filter(Site_original == "Bone") %>% 
  pull(Patient_ID)

primaries <- metadata %>% 
  filter(Sample_type == "Primary")

mets <- metadata %>% 
  filter(Sample_type == "Metastasis") %>% 
  filter(Site_group != "Breast") %>% 
  select(Patient_ID, Site_original) %>% 
  rename(Site_of_Metastasis = Site_original)

plot_meta <- primaries %>% 
  inner_join(mets, by = "Patient_ID") %>%
  mutate(Bone_mets = ifelse(Patient_ID %in% bone_patients, "Bone Met", "Other Met")) %>% 
  filter(Sample_type == "Primary") %>% 
  filter(PAM50_call %in% c("LumA", "LumB")) %>% 
  mutate(Bone_mets = factor(Bone_mets, c("Other Met", "Bone Met")))

plot <- ggplot(plot_meta, mapping = aes(x = Bone_mets, 
                                        y = Luminal_Bone_EN_score)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal(base_size = 14) + 
  stat_compare_means() + 
  geom_jitter(aes(color = Site_of_Metastasis),
              width = 0.2,
              size = 2) + 
  scale_color_manual(values = metastasis_site_colors,
                     na.value = "grey70") +
  xlab("Site of Metastasis") + 
  ylab("Signature score") + 
  ylim(c(0, 0.4)) +
  labs(title = "Luminal Bone EN score\nby site of metastasis")

ggsave(paste0(out_dir, "GEICAM_luminal_bone_module_EN_score_in_lum_primaries_based_on_met_site.png"),
       plot=plot, width = 6, height = 6, dpi = 300)

ggsave(paste0(out_dir, "GEICAM_luminal_bone_module_EN_score_in_lum_primaries_based_on_met_site.pdf"),
       plot=plot, width = 6, height = 6)
