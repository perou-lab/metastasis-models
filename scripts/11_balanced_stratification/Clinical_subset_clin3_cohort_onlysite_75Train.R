# Balanced Stratification of clinical cohorts based on the following parameters:
  ## Original dataset (n=5)
  ## Site of metastasis (only site)

library(tidyverse)
library(sampling)

out_dir <- "output/balanced_stratification/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

# Read in metadata -------------------------------------------------------------
metadata <- read.table("output/PAM50z_subtyping/n2486_combined_meta_w_PAM50_CLOW_subtypes.txt",
                       sep = "\t", header = T, check.names = F)

clinical_subtypes <- unique(metadata$Clinical_Group_3)

for (subtype in clinical_subtypes) {
  subtype_label <- case_when(
    subtype == "HR+/HER2-" ~ "HR+HER2-",
    subtype == "HER2+" ~ "HER2+",
    subtype == "TNBC" ~ "TNBC")
  
  sub_metadata <- metadata %>%
    mutate(cohort_G12276 = ifelse(Cohort == 'GSE12276', 1,0)) %>% 
    mutate(cohort_G2034 = ifelse(Cohort == 'GSE2034',1,0)) %>% 
    mutate(cohort_G2603 = ifelse(Cohort == 'GSE2603',1,0)) %>% 
    mutate(cohort_MB = ifelse(Cohort == 'METABRIC',1,0)) %>% 
    mutate(cohort_NKI = ifelse(Cohort == 'NKI295',1,0)) %>%
    mutate(site_strata = case_when(
      Relapse_Site == "No_Relapse" ~ 0,
      Relapse_Site == "Bone" ~ 1, 
      Relapse_Site == "Brain" ~ 2, 
      Relapse_Site == "Liver" ~ 3, 
      Relapse_Site == "Lung" ~ 4, 
      Relapse_Site == "Multimet" ~ 5)) %>% 
    mutate(site_norel = ifelse(No_v_Distant == "No_Relapse",1,0), 
           site_bone = ifelse(No_v_OnlyBone == "Bone_Relapse",1,0), 
           site_brain = ifelse(No_v_OnlyBrain == "Brain_Relapse",1,0), 
           site_liver = ifelse(No_v_OnlyLiver == "Liver_Relapse",1,0), 
           site_lung = ifelse(No_v_OnlyLung == "Lung_Relapse",1,0), 
           site_multimet = ifelse(No_v_MultiMet == "Multimet_Relapse",1,0)) %>% 
    filter(Clinical_Group_3 == subtype)
  
  
  ## Dataset stratification (Training on 75% of combined dataset) --------------

  # balance variables have to be categorical (0/1)
  balance_variables <- cbind(sub_metadata$cohort_G12276,
                             sub_metadata$cohort_G2034,
                             sub_metadata$cohort_G2603, 
                             sub_metadata$cohort_MB, 
                             sub_metadata$cohort_NKI,
                             sub_metadata$site_norel, 
                             sub_metadata$site_bone, 
                             sub_metadata$site_brain, 
                             sub_metadata$site_liver, 
                             sub_metadata$site_lung, 
                             sub_metadata$site_multimet)
  
  strata <- sub_metadata$site_strata
  # 75/100 is the split between training/testing; times = total number of samples in dataset
  pik <- rep(as.numeric(75/100), times=nrow(sub_metadata))
  
  set.seed(1) # for reproducibility
  
  balanced_split <- balancedstratification(balance_variables, strata, pik, TRUE)
  
  ## Append output to metadata
  meta75 <- sub_metadata %>% 
    mutate(balanced_split = balanced_split) %>% 
    mutate(Train75_Test25 = ifelse(balanced_split == 0, "Test", "Train")) %>% 
    rename(Split75 = balanced_split) %>% 
    dplyr::select(Sample_ID, Relapse_Site,
                  cohort_G12276, cohort_G2034, cohort_G2603, cohort_MB, cohort_NKI,
                  site_norel, site_bone, site_brain, site_liver, site_lung, site_multimet,
                  Clinical_Group_3, Split75, Train75_Test25)
  
  ## Write out metadata for future use -------------------------------------------
  dir.create(paste0(out_dir, subtype_label), showWarnings = F, recursive = T)
  write.table(meta75, paste0(out_dir, subtype_label, "/", subtype_label, "_only_train75_onlysite_stratification.txt"),
              sep = '\t', col.names = T, row.names = F)
  
  ## Combine metadata for input format required for train_EN_site_AM.R -----------
  meta75 <- meta75 %>% 
    select(Sample_ID, Train75_Test25)
  
  merge <- metadata %>% full_join(meta75, by = "Sample_ID")
  
  write.table(merge, paste0(out_dir, subtype_label, "/", subtype_label, "_only_complete_meta_w_onlysite_stratification.txt"),
              sep = "\t", row.names = F, col.names = T)
}

