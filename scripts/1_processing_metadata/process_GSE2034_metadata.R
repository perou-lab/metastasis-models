# Process and unionize metadata from GSE2034

library(tidyverse)
library(readxl)

out_dir <- "output/processing_metadata/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Published molecular data from Supplementary Material 1 of
# Harrell et al. 2012 (Breast Cancer Research and Treatment) PMC3303043
har_data <- read_excel("source_files/clinical_data/10549_2011_1619_MOESM1_ESM_Harrell_2012_SMaterial1.xlsx", sheet = 2, skip = 2) %>% 
  filter(`GEO#` == "GSE2034") %>% 
  rename(Patient_ID = `GEO array#`) %>% 
  mutate(Patient_ID_2 = Patient_ID,
         Sample_ID = Patient_ID, 
         NKI_ID = Patient_ID) %>% 
  mutate(Dataset = Cohort, 
         Dataset_5 = Cohort, 
         Cohort = "GSE2034") %>% 
  rename(MFS_Months = MFS, 
         Relapse_Event = `Any relapse`, 
         Stage = `T Stage`) %>% 
  mutate(Grade = NA_real_, 
         Size = NA_real_,
         ER_status = `ER clincial`, 
         PR_status = `PR clinical`, 
         HER2_status = `Her2 clinical`) %>% 
  rename(Bone_Relapse = `Bone relapse`, 
         Brain_Relapse = `Brain relapse`, 
         Lung_Relapse = `Lung relapse`, 
         Liver_Relapse = `Liver relapse`) %>% 
  mutate(Liver_Relapse = as.numeric(Liver_Relapse))

## Add additional site annotation columns
molecular_data <- har_data %>%
  mutate(No_v_Distant = case_when(
    Relapse_Event == 1 ~ "Distant_Relapse", 
    Relapse_Event == 0 ~ "No_Relapse")) %>% 
  mutate(Four_site_recurrence = case_when(
    Bone_Relapse == 1 & Brain_Relapse == 0 & Liver_Relapse == 0 & Lung_Relapse == 0 ~ "Bone", 
    Bone_Relapse == 0 & Brain_Relapse == 1 & Liver_Relapse == 0 & Lung_Relapse == 0 ~ "Brain",
    Bone_Relapse == 0 & Brain_Relapse == 0 & Liver_Relapse == 1 & Lung_Relapse == 0 ~ "Liver", 
    Bone_Relapse == 0 & Brain_Relapse == 0 & Liver_Relapse == 0 & Lung_Relapse == 1 ~ "Lung",
    Bone_Relapse == 1 & Brain_Relapse == 1 & Liver_Relapse == 0 & Lung_Relapse == 0 ~ "Bone+Brain",
    Bone_Relapse == 1 & Brain_Relapse == 0 & Liver_Relapse == 1 & Lung_Relapse == 0 ~ "Bone+Liver",
    Bone_Relapse == 1 & Brain_Relapse == 0 & Liver_Relapse == 0 & Lung_Relapse == 1 ~ "Bone+Lung",
    Bone_Relapse == 0 & Brain_Relapse == 1 & Liver_Relapse == 1 & Lung_Relapse == 0 ~ "Brain+Liver", 
    Bone_Relapse == 0 & Brain_Relapse == 1 & Liver_Relapse == 0 & Lung_Relapse == 1 ~ "Brain+Lung", 
    Bone_Relapse == 0 & Brain_Relapse == 0 & Liver_Relapse == 1 & Lung_Relapse == 1 ~ "Liver+Lung",
    Bone_Relapse == 1 & Brain_Relapse == 1 & Liver_Relapse == 1 & Lung_Relapse == 0 ~ "Bone+Brain+Liver", 
    Bone_Relapse == 1 & Brain_Relapse == 1 & Liver_Relapse == 0 & Lung_Relapse == 1 ~ "Bone+Brain+Lung", 
    Bone_Relapse == 1 & Brain_Relapse == 0 & Liver_Relapse == 1 & Lung_Relapse == 1 ~ "Bone+Liver+Lung", 
    Bone_Relapse == 0 & Brain_Relapse == 1 & Liver_Relapse == 1 & Lung_Relapse == 1 ~ "Brain+Liver+Lung", 
    Bone_Relapse == 1 & Brain_Relapse == 1 & Liver_Relapse == 1 & Lung_Relapse == 1 ~ "Bone+Brain+Liver+Lung")) %>% 
  mutate(
    n_sites = Bone_Relapse + Brain_Relapse + Liver_Relapse + Lung_Relapse,
    No_v_OnlyBrain = case_when(
      n_sites == 0 & Relapse_Event == 0 ~ "No_Relapse",
      Brain_Relapse == 1 & n_sites == 1 ~ "Brain_Relapse",
      TRUE ~ "Other_Relapse")) %>% 
  mutate(
    No_v_OnlyBone = case_when(
      n_sites == 0 & Relapse_Event == 0 ~ "No_Relapse",
      Bone_Relapse == 1 & n_sites == 1 ~ "Bone_Relapse",
      TRUE ~ "Other_Relapse")) %>% 
  mutate(
    No_v_OnlyLiver = case_when(
      n_sites == 0 & Relapse_Event == 0 ~ "No_Relapse",
      Liver_Relapse == 1 & n_sites == 1 ~ "Liver_Relapse",
      TRUE ~ "Other_Relapse")) %>% 
  mutate(
    No_v_OnlyLung = case_when(
      n_sites == 0 & Relapse_Event == 0 ~ "No_Relapse",
      Lung_Relapse == 1 & n_sites == 1 ~ "Lung_Relapse",
      TRUE ~ "Other_Relapse")) %>%
  mutate(
    Single_v_MultiMet = case_when(
      n_sites == 0 & Relapse_Event == 0 ~ "No_Relapse", 
      n_sites >= 2 ~ "Multimet_Relapse", 
      n_sites == 1 ~ "Single_Relapse",
      TRUE ~ "Other_Relapse")) %>% select(-n_sites) %>% 
  mutate(No_v_MultiMet = case_when(
    Single_v_MultiMet == "Multimet_Relapse" ~ "Multimet_Relapse", 
    Single_v_MultiMet == "No_Relapse" ~ "No_Relapse", 
    TRUE ~ "Other_Relapse")) %>% 
  mutate(Relapse_Site = case_when(
    Four_site_recurrence == "Bone" ~ "Bone", 
    Four_site_recurrence == "Brain" ~ "Brain", 
    Four_site_recurrence == "Lung" ~ "Lung", 
    Four_site_recurrence == "Liver" ~ "Liver", 
    !(Four_site_recurrence %in% c("Bone", "Brain", "Liver", "Lung")) & !is.na(Four_site_recurrence) ~ "Multimet",
    Relapse_Event == 0 ~ "No_Relapse", 
    !(Four_site_recurrence %in% c("Bone", "Brain", "Liver", "Lung")) & is.na(Four_site_recurrence) & Relapse_Event == 1 ~ "Others")) %>% 
  mutate(MultiMet_sites = case_when( 
    Relapse_Site == "Others" ~ "Others", 
    Relapse_Site == "No_Relapse" ~ "No_Relapse",
    TRUE ~ Four_site_recurrence)) %>% 
  select(Patient_ID, Patient_ID_2, Sample_ID, NKI_ID,
         Dataset, Dataset_5, Cohort,
         MFS_Months, Relapse_Event,
         Stage, Grade, Size, ER_status, PR_status, HER2_status,
         Bone_Relapse, Brain_Relapse, Liver_Relapse, Lung_Relapse,
         Four_site_recurrence, Relapse_Site, MultiMet_sites,
         No_v_Distant, No_v_OnlyBone, No_v_OnlyBrain, No_v_OnlyLiver, No_v_OnlyLung, No_v_MultiMet, Single_v_MultiMet)


metadata <- molecular_data %>% 
  filter(Relapse_Site %in% c("Bone", "Brain", "Liver", "Lung", "Multimet", "No_Relapse"))

write.table(metadata, paste0(out_dir, "GSE2034_clinical_data.tsv"), 
            sep = "\t", col.names = T, row.names = F)