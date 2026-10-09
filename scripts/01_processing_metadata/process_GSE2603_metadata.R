# Process and unionize metadata from GSE2603

library(tidyverse)
library(readxl)

out_dir <- "output/processing_metadata/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Published molecular data from Supplementary Data 3 of Minn et al. 2005 (Nature) PMC1283098
minn_data <- read_excel("source_files/clinical_data/Minn_PMC1283098_SData_STable3.xlsx", sheet = 1, skip = 1)  %>% 
  rename(Patient_ID_2 = `Sample Name`) %>% 
  mutate(GEO_set = "GSE2603",
         Dataset = "MSK82", 
         Dataset_5 = "MSK82",
         Cohort = "GSE2603",
         Patient_ID_2 = gsub("-", ".", Patient_ID_2), 
         Stage = NA_real_, 
         Grade = NA_real_, 
         Size = `Tumor Size (cm)`) %>% 
  mutate(ER_status = case_when(
    `Path ER status` == "P" ~ 1, 
    `Path ER status` == "N" ~ 0)) %>% 
  mutate(PR_status = case_when(
    `Path PR status` == "P" ~ 1, 
    `Path PR status` == "N" ~ 0)) %>% 
  mutate(HER2_status = case_when(
    `Her2 status` == "3+" ~ 1, 
    `Her2 status` != "3+" ~ 0))

## Published molecular data from Supplementary Table 8 of Bos et al. 2009 (Nature) PMC1283098
bos_data <- read_excel("source_files/clinical_data/Bos_PMC1283098_STable8.xlsx", sheet = 1) %>% 
  filter(`Patient ID` %in% minn_data$Patient_ID_2) %>% 
  rename(Patient_ID_2 = `Patient ID`, 
         Patient_ID = `GEO asscession number`) %>% 
  mutate(Sample_ID = Patient_ID, 
         NKI_ID = Patient_ID) %>% 
  rename(MFS_Months = `Metastasis-free-survival (months)`,
         Relapse_Event = `All relapses (1=yes, 0=no)`,
         Bone_Relapse = `Bone relapses (1=yes, 0=no)`,
         Brain_Relapse = `Brain relapses (1=yes, 0=no)`, 
         Lung_Relapse = `Lung relapses (1=yes, 0=no)`) %>% 
  mutate(Liver_Relapse = 0) # All samples are missing liver metastasis annotations

## Add additional site annotation columns
molecular_data <- minn_data %>% 
  inner_join(bos_data, by = "Patient_ID_2") %>%
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

write.table(metadata, paste0(out_dir, "GSE2603_clinical_data.tsv"), 
            sep = "\t", col.names = T, row.names = F)