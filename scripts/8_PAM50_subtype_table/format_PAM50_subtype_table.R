# table of PAM50 subtype + site of relapse

library(tidyverse)

out_dir <- "output/PAM50_subtype_table/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Load in metadata with PAM50 subtype annotations
metadata <- read.table("output/PAM50z_subtyping/n2486_combined_meta_w_PAM50_CLOW_subtypes.txt", 
                       sep = "\t", header = T, check.names = F)

## All samples with PAM50 subtype associations
all_table <- as.data.frame(table(metadata$PAM50z_CLOW_cluster)) %>% 
  rename(PAM50_Group = Var1) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2))

all_table_formatted <- all_table %>% 
  mutate(All_Samples = paste0(Freq, " (", Percent, "%)")) %>% 
  dplyr::select(PAM50_Group, All_Samples)

all_samples <- data.frame(PAM50_Group = "All_samples", 
                          All_Samples = "2486 (100%)")

all_table_wide <- rbind(all_table_formatted, all_samples) %>% 
  mutate(PAM50_Group = factor(PAM50_Group, levels = c("All_samples", "Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>%
  arrange(PAM50_Group) %>%
  dplyr::select(PAM50_Group, All_Samples)

## Any distant relapse with PAM50 subtype associations
distant_table <- as.data.frame(table(metadata$No_v_Distant, metadata$PAM50z_CLOW_cluster)) %>% 
  rename(Relapse_Site = Var1, 
         PAM50_Group = Var2) %>% 
  group_by(Relapse_Site) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>% 
  ungroup()

distant_table_all <- as.data.frame(table(metadata$No_v_Distant)) %>% 
  rename(Relapse_Site = Var1) %>% 
  mutate(PAM50_Group = "All_samples") %>% 
  group_by(Relapse_Site) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>% 
  ungroup() %>% 
  select(Relapse_Site, PAM50_Group, Freq, Percent)

distant_table_formatted <- rbind(distant_table, distant_table_all) %>% 
  mutate(Freq_Percent = paste0(Freq, " (", Percent, "%)")) %>%
  dplyr::select(PAM50_Group, Relapse_Site, Freq_Percent)

distant_table_wide <- distant_table_formatted %>%
  pivot_wider(names_from = Relapse_Site, values_from = Freq_Percent) %>% 
  mutate(PAM50_Group = factor(PAM50_Group, levels = c("All_samples", "Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>%
  arrange(PAM50_Group) %>%
  dplyr::select(PAM50_Group, Distant_Relapse)

## Site specific PAM50 subtype associations 
pam50_table <- as.data.frame(table(metadata$Relapse_Site, metadata$PAM50z_CLOW_cluster)) %>%
  rename(Relapse_Site = Var1, 
         PAM50_Group = Var2) %>% 
  group_by(Relapse_Site) %>%
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>%
  ungroup()

pam50_table_all <- as.data.frame(table(metadata$Relapse_Site)) %>% 
  rename(Relapse_Site = Var1) %>% 
  mutate(PAM50_Group = "All_samples") %>% 
  group_by(Relapse_Site) %>% 
  mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>% 
  ungroup() %>% 
  select(Relapse_Site, PAM50_Group, Freq, Percent)

pam50_table_formatted <- rbind(pam50_table, pam50_table_all) %>% 
  mutate(Freq_Percent = paste0(Freq, " (", Percent, "%)")) %>%
  dplyr::select(PAM50_Group, Relapse_Site, Freq_Percent)

pam50_table_wide <- pam50_table_formatted %>%
  pivot_wider(names_from = Relapse_Site, values_from = Freq_Percent) %>% 
  mutate(PAM50_Group = factor(PAM50_Group, levels = c("All_samples", "Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>%
  arrange(PAM50_Group) %>%
  dplyr::select(PAM50_Group, No_Relapse, Bone, Brain, Liver, Lung, Multimet)

table <- merge(all_table_wide, distant_table_wide, by = "PAM50_Group")
table <- merge(table, pam50_table_wide, by = "PAM50_Group")

table <- table %>% 
  select(PAM50_Group, All_Samples, No_Relapse, Distant_Relapse, Bone, Brain,
         Lung, Liver, Multimet)

write.table(table, paste0(out_dir, "/PAM50_Site_of_Relapse_table.txt"), 
            sep = "\t", col.names = T, row.names = F)
