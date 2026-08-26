# Balanced Stratification of EN cohort based on the following parameters:
  ## Original dataset (n=5)
  ## Clinical group (HR+/HER2-, HER2+, TNBC)
  ## Site of metastasis (no relapse, only site - brain, bone, liver, lung, multimet)

library(tidyverse)
library(sampling)

out_dir <- "output/balanced_stratification/n2486/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Read in metadata ------------------------------------------------------------
metadata <- read.table("output/PAM50z_subtyping/n2486_combined_meta_w_PAM50_CLOW_subtypes.txt",
                       sep = "\t", header = T, check.names = F)

metadata <- metadata %>%
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
  mutate(Clinical_HR = ifelse(Clinical_Group_3 == "HR+/HER2-",1,0)) %>% 
  mutate(Clinical_HER2 = ifelse(Clinical_Group_3 == "HER2+",1,0)) %>% 
  mutate(Clinical_TNBC = ifelse(Clinical_Group_3 == "TNBC",1,0)) %>% 
  mutate(site_norel = ifelse(No_v_Distant == "No_Relapse",1,0), 
         site_bone = ifelse(No_v_OnlyBone == "Bone_Relapse",1,0), 
         site_brain = ifelse(No_v_OnlyBrain == "Brain_Relapse",1,0), 
         site_liver = ifelse(No_v_OnlyLiver == "Liver_Relapse",1,0), 
         site_lung = ifelse(No_v_OnlyLung == "Lung_Relapse",1,0), 
         site_multimet = ifelse(No_v_MultiMet == "Multimet_Relapse",1,0))
                                                                                   

## Dataset stratification (Training on 75% of combined dataset) ----------------

# balance variables have to be categorical (0/1)
balance_variables <- cbind(metadata$cohort_G12276,
                           metadata$cohort_G2034,
                           metadata$cohort_G2603, 
                           metadata$cohort_MB, 
                           metadata$cohort_NKI,
                           metadata$Clinical_HR, 
                           metadata$Clinical_HER2, 
                           metadata$Clinical_TNBC, 
                           metadata$site_norel, 
                           metadata$site_bone, 
                           metadata$site_brain, 
                           metadata$site_liver, 
                           metadata$site_lung, 
                           metadata$site_multimet)

strata <- metadata$site_strata
# 75/100 is the split between training/testing; times = total number of samples in dataset
pik <- rep(as.numeric(75/100), times=nrow(metadata))

set.seed(1) # for reproducibility

balanced_split <- balancedstratification(balance_variables, strata, pik, TRUE)

## Append output to metadata
meta75 <- metadata %>% 
  mutate(balanced_split = balanced_split) %>% 
  mutate(Train75_Test25 = ifelse(balanced_split == 0, "Test", "Train")) %>% 
  rename(Split75 = balanced_split) %>% 
  dplyr::select(Sample_ID, Relapse_Site,
                cohort_G12276, cohort_G2034, cohort_G2603, cohort_MB, cohort_NKI,
                site_norel, site_bone, site_brain, site_liver, site_lung, site_multimet,
                Clinical_Group_3, Clinical_HR, Clinical_HER2, Clinical_TNBC, 
                Split75, Train75_Test25)

## Write out metadata for future use -------------------------------------------
write.table(meta75, paste0(out_dir, "/n2486_train75_onlysite_stratification.txt"),
                           sep = '\t', col.names = T, row.names = F)


## Combine metadata for input format required for train_EN_site_AM.R -----------
meta75 <- meta75 %>% 
  select(Sample_ID, Train75_Test25)

merge <- metadata %>% full_join(meta75, by = "Sample_ID")

write.table(merge, paste0(out_dir, "n2486_complete_meta_w_onlysite_stratification.txt"), 
            sep = "\t", row.names = F, col.names = T)
