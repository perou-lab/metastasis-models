# table of clinical subtype + site of relapse

library(tidyverse)

out_dir <- "output/clinical_subtype_table/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Load in metadata with clinical subtype annotations
metadata <- read.table("output/processing_metadata/n2486_combined_clinical_data.tsv", 
                       sep = "\t", header = T, check.names = F) %>% 
  select(Sample_ID, Relapse_Site, No_v_Distant)

receptor_data <- read.table("output/infer_ER_PR_HER2/n2486_Ensembl_dwd_inferred_clin_groups.txt", 
                            sep = "\t", header = T) %>% 
  mutate(Clinical_Group_3 = case_when(
    (ER_status_w_inferred == "Positive" | PR_status_w_inferred == "Positive") & HER2_status_w_inferred == "Negative" ~ "HR+/HER2-", 
    ER_status_w_inferred == "Negative" & PR_status_w_inferred == "Negative" & HER2_status_w_inferred == "Negative" ~ "TNBC", 
    HER2_status_w_inferred == "Positive" ~ "HER2+"))

metadata <- left_join(metadata, receptor_data, by = "Sample_ID")

## All samples with clinical subtype associations
all_table <- as.data.frame(table(metadata$Clinical_Group_3)) %>% 
  rename(Clinical_Group = Var1) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2))

all_table_formatted <- all_table %>% 
  mutate(All_Samples = paste0(Freq, " (", Percent, "%)")) %>% 
  dplyr::select(Clinical_Group, All_Samples)

all_samples <- data.frame(Clinical_Group = "All_samples", 
                          All_Samples = "2486 (100%)")

all_table_wide <- rbind(all_table_formatted, all_samples) %>% 
  mutate(Clinical_Group = factor(Clinical_Group, levels = c("All_samples","HR+/HER2-", "HER2+", "TNBC"))) %>%
  arrange(Clinical_Group) %>%
  dplyr::select(Clinical_Group, All_Samples)

## Any distant relapse with clinical subtype associations
distant_table <- as.data.frame(table(metadata$No_v_Distant, metadata$Clinical_Group_3)) %>% 
  rename(Relapse_Site = Var1, 
         Clinical_Group = Var2) %>% 
  group_by(Relapse_Site) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>% 
  ungroup()

distant_table_all <- as.data.frame(table(metadata$No_v_Distant)) %>% 
  rename(Relapse_Site = Var1) %>% 
  mutate(Clinical_Group = "All_samples") %>% 
  group_by(Relapse_Site) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>% 
  ungroup() %>% 
  select(Relapse_Site, Clinical_Group, Freq, Percent)

distant_table_formatted <- rbind(distant_table, distant_table_all) %>% 
  mutate(Freq_Percent = paste0(Freq, " (", Percent, "%)")) %>%
  dplyr::select(Clinical_Group, Relapse_Site, Freq_Percent)

distant_table_wide <- distant_table_formatted %>%
  pivot_wider(names_from = Relapse_Site, values_from = Freq_Percent) %>% 
  mutate(Clinical_Group = factor(Clinical_Group, levels = c("All_samples", "HR+/HER2-", "HER2+", "TNBC"))) %>%
  arrange(Clinical_Group) %>%
  dplyr::select(Clinical_Group, Distant_Relapse)

## Site specific clinical subtype associations 
ihc_table <- as.data.frame(table(metadata$Relapse_Site, metadata$Clinical_Group_3)) %>%
  rename(Relapse_Site = Var1, 
         Clinical_Group = Var2) %>% 
  group_by(Relapse_Site) %>%
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>%
  ungroup()

ihc_table_all <- as.data.frame(table(metadata$Relapse_Site)) %>% 
  rename(Relapse_Site = Var1) %>% 
  mutate(Clinical_Group = "All_samples") %>% 
  group_by(Relapse_Site) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>% 
  ungroup() %>% 
  select(Relapse_Site, Clinical_Group, Freq, Percent)

ihc_table_formatted <- rbind(ihc_table, ihc_table_all) %>% 
  mutate(Freq_Percent = paste0(Freq, " (", Percent, "%)")) %>%
  dplyr::select(Clinical_Group, Relapse_Site, Freq_Percent)

ihc_table_wide <- ihc_table_formatted %>%
  pivot_wider(names_from = Relapse_Site, values_from = Freq_Percent) %>% 
  mutate(Clinical_Group = factor(Clinical_Group, levels = c("All_samples", "HR+/HER2-", "HER2+", "TNBC"))) %>%
  arrange(Clinical_Group) %>%
  dplyr::select(Clinical_Group, No_Relapse, Bone, Brain, Liver, Lung, Multimet)

table <- merge(all_table_wide, distant_table_wide, by = "Clinical_Group")
table <- merge(table, ihc_table_wide, by = "Clinical_Group")

table <- table %>% 
  select(Clinical_Group, All_Samples, No_Relapse, Distant_Relapse, Bone, Brain,
         Lung, Liver, Multimet)

write.table(table, paste0(out_dir, "/Clinical_group_Site_of_Relapse_table.txt"), 
            sep = "\t", col.names = T, row.names = F)
