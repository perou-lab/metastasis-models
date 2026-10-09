# Process and unionize metadata from METABRIC

library(tidyverse)

out_dir <- "output/processing_metadata/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Published molecular data from Supplementary Table 6 of Rueda et al. 2019 (Nature) PMC6647838
molecular_data <- read.table("source_files/clinical_data/41586_2019_1007_MOESM8_ESM _Rueda_2019_STable6.txt", 
                             sep = "\t", header = T, check.names = F) %>% 
  dplyr::select(METABRIC.ID, Grade, Size, Stage, TDR, DR, ER.Expr, Her2.Expr, PR.Expr, Complete.Rec.History) %>% 
  dplyr::rename(MFS_Days = TDR, # Time until last follow up or first distant relapse (days)
                Relapse_event = DR, # Distant relapse (1=YES, 0=NO)
                ER_status_by_expr = ER.Expr, 
                HER2_status_by_expr = Her2.Expr, 
                PR_status_by_expr = PR.Expr) %>% 
  mutate(ER_status_by_expr = case_when(ER_status_by_expr == "+" ~ 1, 
                                       ER_status_by_expr == "-" ~ 0, 
                                       is.na(ER_status_by_expr) ~ NA_integer_), 
         PR_status_by_expr = case_when(PR_status_by_expr == "+" ~ 1, 
                                       PR_status_by_expr == "-" ~ 0, 
                                       is.na(PR_status_by_expr) ~ NA_integer_), 
         HER2_status_by_expr = case_when(HER2_status_by_expr == "+" ~ 1, 
                                         HER2_status_by_expr == "-" ~ 0, 
                                         is.na(HER2_status_by_expr) ~ NA_integer_))

# Published recurrence data from Supplementary Table 7 of Rueda et al. 2019 (Nature) PMC6647838
metastasis_sites <- c("BONE", "BREAST", "LIVER", "LNS", "OVARY", "PLEURA", 
                      "ASCITES", "SKIN", "LUNG", "PERICARDIAL", "SOFT_TISSUES", 
                      "BRAIN", "MEDIASTINAL", "MENINGEAL", "ADRENALS", "ABDOMEN", 
                      "GI_TRACT", "BLADDER", "EYE", "GALLBLADDER",
                      "UTERUS_AND_FALLOPIAN_TUBE", "KIDNEY", "PERITONEUM",
                      "PANCREAS", "OTHER")

four_sites <- c("BONE", "BRAIN", "LIVER", "LUNG")

recurrence_data <- read.table("source_files/clinical_data/41586_2019_1007_MOESM9_ESM_Rueda_2019_STable7.txt",
                              sep = "\t", header = TRUE, check.names = FALSE) %>% 
  dplyr::select(METABRIC.ID, all_of(metastasis_sites)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(Relapse_site_all = { 
    sites <- metastasis_sites[c_across(all_of(metastasis_sites)) == 1]
    if (length(sites) == 0) NA_character_ else paste(sites, collapse = ", ")},
    
    Four_site_recurrence = {core_hits <- four_sites[c_across(all_of(four_sites)) == 1]
      if (length(core_hits) > 0) {
        paste(core_hits, collapse = "+")
      } else if (any(c_across(all_of(metastasis_sites)) == 1)) {
        "Other_site"
      } else {
        "Recurrent_cohort_no_site"
      }
    }) %>% 
  ungroup() %>% 
  dplyr::select(METABRIC.ID, Relapse_site_all, Four_site_recurrence, BONE, BRAIN, LIVER, LUNG) %>% 
  rename(Bone_Relapse = BONE, 
         Brain_Relapse = BRAIN, 
         Liver_Relapse = LIVER, 
         Lung_Relapse = LUNG)

## Identify patients with locoregional recurrence
lrr_patients <- read.table("source_files/clinical_data/41586_2019_1007_MOESM9_ESM_Rueda_2019_STable7.txt",
                           sep = "\t", header = TRUE, check.names = FALSE) %>% 
  filter(TYPE.RELAPSE == "LOCOREGIONAL") %>% 
  distinct() %>% 
  pull(METABRIC.ID)

## Join molecular and recurrence data
merged <- left_join(molecular_data, recurrence_data, by = "METABRIC.ID") %>% 
  rename(Complete_recurrence_history = Complete.Rec.History) %>% 
  mutate(Locoregional_recurrence = ifelse(METABRIC.ID %in% lrr_patients, "Local_Relapse", "No_Local_Relapse")) %>% 
  mutate(Relapse_site_4 = case_when(
    Four_site_recurrence == "BONE" ~ "Bone",
    Four_site_recurrence == "BRAIN" ~ "Brain",
    Four_site_recurrence == "LUNG" ~ "Lung", 
    Four_site_recurrence == "LIVER" ~ "Liver",
    grepl("\\+", Four_site_recurrence) ~ "Multimet",
    Four_site_recurrence == "Other_site" ~ "Other_Distant_Relapse",
    is.na(Four_site_recurrence) ~ "No_Relapse",
    Four_site_recurrence == "Recurrent_cohort_no_site" & Locoregional_recurrence == "Local_Relapse" ~ "Local_Relapse",
    TRUE ~ Relapse_site_all)) %>% 
  filter(Complete_recurrence_history == "YES")

## Subset for sites of interest: 
metadata <- merged %>% 
  filter(Relapse_site_4 %in% c("Bone", "Brain", "Liver", "Lung", "Multimet", "No_Relapse")) %>% 
  dplyr::rename(Patient_ID = METABRIC.ID) %>% 
  mutate(Dataset = "METABRIC", 
         GEO_set = "METABRIC",
         MFS_Months = MFS_Days/(365.25/12))

## Match formatting for HARRELL metadata
metadata <- metadata %>% 
  mutate(Patient_ID = gsub("-", ".", Patient_ID), 
         Patient_ID_2 = Patient_ID, 
         Sample_ID = Patient_ID, 
         NKI_ID = Patient_ID, 
         Dataset = "METABRIC", 
         Dataset_5 = "METABRIC", 
         Cohort = "METABRIC") %>% 
  rename(Relapse_Event = Relapse_event, 
         ER_status = ER_status_by_expr,
         PR_status = PR_status_by_expr, 
         HER2_status = HER2_status_by_expr, 
         Relapse_Site = Relapse_site_4) %>% 
  mutate(Bone_Relapse = ifelse(Relapse_Event == 0 & is.na(Bone_Relapse), 0, Bone_Relapse), 
         Brain_Relapse = ifelse(Relapse_Event == 0 & is.na(Brain_Relapse), 0, Brain_Relapse), 
         Lung_Relapse = ifelse(Relapse_Event == 0 & is.na(Lung_Relapse), 0, Lung_Relapse), 
         Liver_Relapse = ifelse(Relapse_Event == 0 & is.na(Liver_Relapse), 0, Liver_Relapse)) %>% 
  mutate(Four_site_recurrence = Four_site_recurrence %>%
           gsub("BONE", "Bone", ., fixed = TRUE) %>%
           gsub("BRAIN", "Brain", ., fixed = TRUE) %>%
           gsub("LIVER", "Liver", ., fixed = TRUE) %>%
           gsub("LUNG", "Lung", ., fixed = TRUE)) %>% 
  mutate(MultiMet_sites = case_when(
    Relapse_Site %in% c("Bone", "Brain", "Liver", "Lung") ~ Relapse_Site, 
    Relapse_Site == "Multimet" ~ Four_site_recurrence, 
    Relapse_Site == "No_Relapse" ~ "No_Relapse")) %>% 
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
  mutate(No_v_Distant = case_when(
    Relapse_Event == 1 ~ "Distant_Relapse", 
    Relapse_Event == 0 ~ "No_Relapse")) %>% 
  select(Patient_ID, Patient_ID_2, Sample_ID, NKI_ID,
         Dataset, Dataset_5, Cohort,
         MFS_Months, Relapse_Event,
         Stage, Grade, Size, ER_status, PR_status, HER2_status,
         Bone_Relapse, Brain_Relapse, Liver_Relapse, Lung_Relapse,
         Four_site_recurrence, Relapse_Site, MultiMet_sites,
         No_v_Distant, No_v_OnlyBone, No_v_OnlyBrain, No_v_OnlyLiver, No_v_OnlyLung, No_v_MultiMet, Single_v_MultiMet)

write.table(metadata, paste0(out_dir, "METABRIC_clinical_data.tsv"), 
            sep = "\t", col.names = T, row.names = F)
