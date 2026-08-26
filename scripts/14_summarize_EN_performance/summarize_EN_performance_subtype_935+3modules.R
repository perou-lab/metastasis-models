# Evaluate performance of EN models (including subtype specific) using modules (935+3) as input

library(ggplot2)
library(tidyverse)

out_dir <- "output/summarize_EN_performance//"
dir.create(paste0(out_dir), recursive = T, showWarnings = F)

## Combine n2486 summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/n2486/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
    df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

n2486_data <- do.call(rbind, data_list)
n2486_data <- n2486_data %>% 
  mutate(subset = "n2486")


## Combine basal summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/basal_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

basal_data <- do.call(rbind, data_list)
basal_data <- basal_data %>% 
  mutate(subset = "Basal")


## Combine lumA summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/lumA_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

lumA_data <- do.call(rbind, data_list)
lumA_data <- lumA_data %>% 
  mutate(subset = "LumA")


## Combine lumB summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/lumB_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

lumB_data <- do.call(rbind, data_list)
lumB_data <- lumB_data %>% 
  mutate(subset = "LumB")


## Combine HER2E summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/HER2E_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

HER2E_data <- do.call(rbind, data_list)
HER2E_data <- HER2E_data %>% 
  mutate(subset = "HER2E")


## Combine LumA + LumB summary ## ----------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/lumA_lumB_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

LumAB_data <- do.call(rbind, data_list)
LumAB_data <- LumAB_data %>% 
  mutate(subset = "LumA_LumB")


## Combine HR+/HER2- summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/HR+HER2_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

HR_data <- do.call(rbind, data_list)
HR_data <- HR_data %>% 
  mutate(subset = "HR+/HER2-")


## Combine HER2+ summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/HER2+_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

HER2_data <- do.call(rbind, data_list)
HER2_data <- HER2_data %>% 
  mutate(subset = "HER2+")


## Combine TNBC summary ## ----------------------------------------------------
root_dir <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/TNBC_only/modules_935+3/"
sub_folders <- list.dirs(root_dir, full.names = T, recursive = F)

data_list <- list()

for (folder in sub_folders) {
  feature_summary_dir <- file.path(folder, "feature_summary")
  txt_files <- list.files(feature_summary_dir,
                          pattern = "\\.txt$",
                          full.names = T)
  
  if(length(txt_files) == 1) {
    df <- read.table(txt_files[1], header = T, sep = "\t", stringsAsFactors = F)
    data_list[[length(data_list) + 1]] <- df
  } else if (length(txt_files > 1)) {
    match_file <- txt_files[grepl("124022.txt", txt_files)]
    
    if (length(match_file) == 1) {
      df <- read.table(match_file, header = T, sep = "\t", stringsAsFactors = F)
      data_list[[length(data_list) + 1]] <- df
    } else {
      warning(paste0("File with seed = 124022 is not uniquely found in ", feature_summary_dir))
    }
  } else {
    warning(paste0("No .txt files found in ", feature_summary_dir))
  }
}

TNBC_data <- do.call(rbind, data_list)
TNBC_data <- TNBC_data %>% 
  mutate(subset = "TNBC")


## All combined ----------------------------------------------------------------        
all_data <- rbind(n2486_data, HR_data, HER2_data, TNBC_data, 
                  basal_data, HER2E_data, lumA_data, lumB_data, LumAB_data)

save_data <- all_data

write.table(all_data, paste0(out_dir, "Results_summary_s124022_935+3_modules.txt"),
            sep = "\t", row.names = F, col.names = T)


## Plot results ----------------------------------------------------------------
dataset_order <- rev(c("n2486", "HR+/HER2-", "HER2+", "TNBC", "Basal", "HER2E", "LumA", "LumB", "LumA_LumB"))
site_order <- c("Distant", "Distant_notBone", "OnlyBone", "OnlyBrain", "OnlyLiver", "OnlyLung", "MultiMet")

## Including all results (not filtering on sample size)
plot_data <- all_data %>%
  select(subset, Met_Site, Testing_AUC) %>%
  mutate(subset = factor(subset, levels = dataset_order)) %>%
  mutate(Met_Site = factor(Met_Site, levels = site_order)) %>%
  mutate(Fill_Color = ifelse(Testing_AUC < 0.5, 0.5, Testing_AUC))

p <- ggplot(plot_data, aes(x = Met_Site, y = subset, fill = Fill_Color)) +
  geom_tile(color = "black") +
  geom_text(data = filter(plot_data, !is.na(Fill_Color) & Fill_Color != 0.5),
            aes(label = round(Testing_AUC, 2)),
            size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "darkblue",
                      limits = c(0.5,1),
                      na.value = "gray",
                      name = "Test AUC") +
  labs(x = "Metastatic Site", y = "Dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Elastic Net Model Test AUC Performance (935+3 modules)") +
  coord_fixed()

png(paste0(out_dir, "EN_935+3modules_Test_AUC_subtype_grid_noNA.png"),
    width = 6.5, height = 6.5, units = "in", res = 300)
print(p)
dev.off()

pdf(paste0(out_dir, "EN_935+3modules_Test_AUC_subtype_grid_noNA.pdf"), 
    width = 6.5, height = 6.5)
print(p)
dev.off()


## Including gray squares for comparisons with low sample size
n2486_sample_sets <- read.table("output/balanced_stratification/n2486/n2486_complete_meta_w_onlysite_stratification.txt", 
                                sep = "\t", header = T, check.names = F)
basal_sample_sets <- read.table("output/balanced_stratification/basal_only/basal_only_complete_meta_w_onlysite_stratification.txt", 
                                sep = "\t", header = T, check.names = F)
her2e_sample_sets <- read.table("output/balanced_stratification/HER2E_only/HER2E_only_complete_meta_w_onlysite_stratification.txt", 
                                sep = "\t", header = T, check.names = F)
lumab_sample_sets <- read.table("output/balanced_stratification/lumA_lumB_only/lumA_lumB_only_complete_meta_w_onlysite_stratification.txt", 
                                sep = "\t", header = T, check.names = F)
luma_sample_sets <- read.table("output/balanced_stratification/lumA_only/lumA_only_complete_meta_w_onlysite_stratification.txt", 
                               sep = "\t", header = T, check.names = F)
lumb_sample_sets <- read.table("output/balanced_stratification/lumB_only/lumB_only_complete_meta_w_onlysite_stratification.txt", 
                               sep = "\t", header = T, check.names = F)
her2_sample_sets <- read.table("output/balanced_stratification/HER2+/HER2+_only_complete_meta_w_onlysite_stratification.txt", 
                               sep = "\t", header = T, check.names = F)
tnbc_sample_sets <- read.table("output/balanced_stratification/TNBC/TNBC_only_complete_meta_w_onlysite_stratification.txt", 
                               sep = "\t", header = T, check.names = F)
hr_sample_sets <- read.table("output/balanced_stratification/HR+HER2-/HR+HER2-_only_complete_meta_w_onlysite_stratification.txt", 
                             sep = "\t", header = T, check.names = F)

comparison_map <- tibble(
  Met_Site = c("Distant", "Distant_notBone", "OnlyLung", "OnlyBone", "OnlyBrain", "OnlyLiver", "Multimet"),
  comparison_col = c("No_v_Distant", "No_v_Distant", "No_v_OnlyLung", "No_v_OnlyBone", "No_v_OnlyBrain",
                     "No_v_OnlyLiver", "No_v_MultiMet"),
  positive_class = c("Distant_Relapse", "Distant_Relapse", "Lung_Relapse", "Bone_Relapse", "Brain_Relapse",
                     "Liver_Relapse", "Multimet_Relapse"))

sample_set_list <- list(
  n2486 = n2486_sample_sets,
  Basal = basal_sample_sets,
  HER2E = her2e_sample_sets,
  LumA_LumB = lumab_sample_sets,
  LumA = luma_sample_sets,
  LumB = lumb_sample_sets,
  `HER2+` = her2_sample_sets,
  TNBC = tnbc_sample_sets,
  `HR+/HER2-` = hr_sample_sets)


test_counts <- imap_dfr(sample_set_list, function(df, subset_name) {
  
  map_dfr(seq_len(nrow(comparison_map)), function(i) {
    
    dat <- df %>%
      filter(
        Train75_Test25 == "Test",
        .data[[comparison_map$comparison_col[i]]] ==
          comparison_map$positive_class[i]
      )
    
    # Special remove bone samples for Distant_notBone comparison
    if (comparison_map$Met_Site[i] == "Distant_notBone") {
      dat <- dat %>% filter(Relapse_Site != "Bone")
    }
    
    tibble(
      subset = subset_name,
      Met_Site = comparison_map$Met_Site[i],
      nTest = nrow(dat)
    )
  })
  
})

test_counts <- test_counts %>% 
  mutate(Met_Site = ifelse(Met_Site == "Multimet", "MultiMet", Met_Site)) %>% 
  mutate(merge_name = paste0(subset, "_", Met_Site))

plot_data <- all_data %>% 
  mutate(merge_name = paste0(subset, "_", Met_Site)) %>% 
  select(merge_name, Testing_AUC) %>% 
  right_join(test_counts, by = "merge_name") %>% 
  mutate(subset = factor(subset, levels = dataset_order)) %>% 
  mutate(Met_Site = factor(Met_Site, levels = site_order)) %>% 
  mutate(Fill_Color = ifelse(nTest < 10, NA_real_, Testing_AUC), 
         Fill_Color = ifelse(!is.na(Fill_Color) & Fill_Color < 0.5, 0.5, Fill_Color))

p <- ggplot(plot_data, aes(x = Met_Site, y = subset, fill = Fill_Color)) +
  geom_tile(color = "black") + 
  geom_text(data = filter(plot_data, !is.na(Fill_Color) & Fill_Color != 0.5), 
            aes(label = round(Testing_AUC, 2)), 
            size = 3, color = "black") + 
  scale_fill_gradient(low = "white", high = "darkblue",
                      limits = c(0.5,1), 
                      na.value = "gray", 
                      name = "Test AUC") + 
  labs(x = "Metastatic Site", y = "Dataset") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Elastic Net Model Test AUC Performance (935+3 modules)") + 
  coord_fixed()
  
png(paste0(out_dir, "EN_935+3modules_Test_AUC_subtype_grid.png"), 
    width = 6.5, height = 6.5, units = "in", res = 300)
print(p)
dev.off()

pdf(paste0(out_dir, "EN_935+3modules_Test_AUC_subtype_grid.pdf"), 
    width = 6.5, height = 6.5)
print(p)
dev.off()


## Plot number of test samples for each comparison

p <- ggplot(plot_data, aes(x = Met_Site, y = subset, fill = Fill_Color)) +
  geom_tile(color = "black") + 
  geom_text(data = filter(plot_data), 
            aes(label = nTest), 
            size = 3, color = "black") + 
  scale_fill_gradient(low = "white", high = "darkblue",
                      limits = c(0.5,1), 
                      na.value = "gray", 
                      name = "Test AUC") + 
  labs(x = "Metastatic Site", y = "Dataset") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Elastic Net Model Test AUC Performance\nwith number of test samples (935+3 modules)") + 
  coord_fixed()

png(paste0(out_dir, "EN_935+3modules_Test_n_samples_subtype_grid.png"), 
    width = 6.5, height = 6.5, units = "in", res = 300)
print(p)
dev.off()

pdf(paste0(out_dir, "EN_935+3modules_Test_n_samples_subtype_grid.pdf"), 
    width = 6.5, height = 6.5)
print(p)
dev.off()

