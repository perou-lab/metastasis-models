## Plot top features for publication

library(tidyverse)
library(ggplot2)

## Distant metastasis (signatures) ---------------------------------------------

## Read in model with features and coefficients
plot_data <- read.table("output/train_EN_models/n2486/modules_935+3/NoRel-v-Distant/model_txt/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-Distant_boot_accuracy.txt", 
                       sep = "\t", header = T, check.names = F) %>% 
  mutate(feature = as.character(feature)) %>% 
  mutate(scaled_coef = as.numeric(scaled_coef)) %>% 
  dplyr::select(feature, scaled_coef)

## Identify top features (n <= 30 features up/down)
top_up <- plot_data %>%
  filter(scaled_coef > 0) %>%
  arrange(desc(scaled_coef)) %>%
  slice_head(n = 30)

top_down <- plot_data %>%
  filter(scaled_coef < 0) %>%
  arrange(scaled_coef) %>%
  slice_head(n = 30)

plot_data_plot <- bind_rows(top_up, top_down) %>% 
  mutate(feature = factor(feature, levels = unique(feature[order(scaled_coef, decreasing = FALSE)])))

## Plot top features
p <- ggplot(plot_data_plot, aes(x = feature, y = scaled_coef, fill = scaled_coef > 0)) + 
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = c("TRUE" = "darkblue", "FALSE" = "firebrick")) + 
  coord_flip() + 
  labs(title = paste("Distant Metastasis EN Model (Signatures)"),
       y = "Scaled B-coefficient") + 
  theme_minimal(base_size = 14)

png(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/n2486/modules_935+3/NoRel-v-Distant/top_features_plot/EN_NoRel-v-Distant_top_features_pub_modules_75Train_25Test_124022_FORMAT.png"),
    width = 10, height = 12, units = "in", res = 350)
print(p)
dev.off()

pdf(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/n2486/modules_935+3/NoRel-v-Distant/top_features_plot/EN_NoRel-v-Distant_top_features_pub_modules_75Train_25Test_124022_FORMAT.pdf"),
    width = 10, height = 12)
print(p)
dev.off()



## Luminal Bone metastasis (signatures) ----------------------------------------

## Read in model with features and coefficients
plot_data <- read.table("output/train_EN_models/lumA_lumB_only/modules_935+3/NoRel-v-OnlyBone/model_txt/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-OnlyBone_boot_accuracy.txt", 
                        sep = "\t", header = T, check.names = F) %>% 
  mutate(feature = as.character(feature)) %>% 
  mutate(scaled_coef = as.numeric(scaled_coef)) %>% 
  dplyr::select(feature, scaled_coef)

## Identify top features (n <= 30 features up/down)
top_up <- plot_data %>%
  filter(scaled_coef > 0) %>%
  arrange(desc(scaled_coef)) %>%
  slice_head(n = 30)

top_down <- plot_data %>%
  filter(scaled_coef < 0) %>%
  arrange(scaled_coef) %>%
  slice_head(n = 30)

plot_data_plot <- bind_rows(top_up, top_down) %>% 
  mutate(feature = factor(feature, levels = unique(feature[order(scaled_coef, decreasing = FALSE)])))

## Plot top features
p <- ggplot(plot_data_plot, aes(x = feature, y = scaled_coef, fill = scaled_coef > 0)) + 
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = c("TRUE" = "darkblue", "FALSE" = "firebrick")) + 
  coord_flip() + 
  labs(title = paste("Luminal Bone Metastasis EN Model (Signatures)"),
       y = "Scaled B-coefficient") + 
  theme_minimal(base_size = 16)

png(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/lumA_lumB_only/modules_935+3/NoRel-v-OnlyBone/top_features_plot/EN_NoRel-v-OnlyBone_top_features_pub_modules_75Train_25Test_124022_FORMAT.png"),
    width = 10, height = 12, units = "in", res = 350)
print(p)
dev.off()

pdf(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/lumA_lumB_only/modules_935+3/NoRel-v-OnlyBone/top_features_plot/EN_NoRel-v-OnlyBone_top_features_pub_modules_75Train_25Test_124022_FORMAT.pdf"),
    width = 10, height = 12)
print(p)
dev.off()


## Multi-site metastasis (signatures) ------------------------------------------

## Read in model with features and coefficients
plot_data <- read.table("output/train_EN_models/n2486/modules_935+3/NoRel-v-MultiMet/model_txt/EN_75Train_25Test_pub_modules-modules_glmnet_obj_s124022_NoRel-v-MultiMet_boot_accuracy.txt", 
                        sep = "\t", header = T, check.names = F) %>% 
  mutate(feature = as.character(feature)) %>% 
  mutate(scaled_coef = as.numeric(scaled_coef)) %>% 
  dplyr::select(feature, scaled_coef)

## Identify top features (n <= 30 features up/down)
top_up <- plot_data %>%
  filter(scaled_coef > 0) %>%
  arrange(desc(scaled_coef)) %>%
  slice_head(n = 30)

top_down <- plot_data %>%
  filter(scaled_coef < 0) %>%
  arrange(scaled_coef) %>%
  slice_head(n = 30)

plot_data_plot <- bind_rows(top_up, top_down) %>% 
  mutate(feature = factor(feature, levels = unique(feature[order(scaled_coef, decreasing = FALSE)])))

## Plot top features
p <- ggplot(plot_data_plot, aes(x = feature, y = scaled_coef, fill = scaled_coef > 0)) + 
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = c("TRUE" = "darkblue", "FALSE" = "firebrick")) + 
  coord_flip() + 
  labs(title = paste("Multi-Site Metastasis EN Model (Signatures)"),
       y = "Scaled B-coefficient") + 
  theme_minimal(base_size = 16)

png(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/n2486/modules_935+3/NoRel-v-MultiMet/top_features_plot/EN_NoRel-v-MultiMet_top_features_pub_modules_75Train_25Test_124022_FORMAT.png"),
    width = 10, height = 12, units = "in", res = 350)
print(p)
dev.off()

pdf(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/n2486/modules_935+3/NoRel-v-MultiMet/top_features_plot/EN_NoRel-v-MultiMet_top_features_pub_modules_75Train_25Test_124022_FORMAT.pdf"),
    width = 10, height = 12)
print(p)
dev.off()


## Distant metastasis (genes) --------------------------------------------------

## Read in model with features and coefficients
plot_data <- read.table("output/train_EN_models/n2486/genes_10814/NoRel-v-Distant/model_txt/EN_75Train_25Test_dwd-genes_glmnet_obj_s124022_NoRel-v-Distant_boot_accuracy.txt", 
                        sep = "\t", header = T, check.names = F) %>% 
  mutate(feature = as.character(GeneSymbol)) %>% 
  mutate(coef = as.numeric(unscaled_coef)) %>% 
  dplyr::select(feature, coef)

## Identify top features (n <= 30 features up/down)
top_up <- plot_data %>%
  filter(coef > 0) %>%
  arrange(desc(coef)) %>%
  slice_head(n = 30)

top_down <- plot_data %>%
  filter(coef < 0) %>%
  arrange(coef) %>%
  slice_head(n = 30)

plot_data_plot <- bind_rows(top_up, top_down) %>% 
  mutate(feature = factor(feature, levels = unique(feature[order(coef, decreasing = FALSE)])))

## Plot top features
p <- ggplot(plot_data_plot, aes(x = feature, y = coef, fill = coef > 0)) + 
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = c("TRUE" = "darkblue", "FALSE" = "firebrick")) + 
  coord_flip() + 
  labs(title = paste("Distant Metastasis EN Model (Genes)"),
       y = "B-coefficient") + 
  theme_minimal(base_size = 16)

png(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/n2486/genes_10814/NoRel-v-Distant/top_features_plot/EN_NoRel-v-Distant_top_features_dwd_75Train_25Test_124022_FORMAT.png"),
    width = 4, height = 12, units = "in", res = 350)
print(p)
dev.off()

pdf(paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/n2486/genes_10814/NoRel-v-Distant/top_features_plot/EN_NoRel-v-Distant_top_features_dwd_75Train_25Test_124022_FORMAT.pdf"),
    width = 4, height = 12)
print(p)
dev.off()