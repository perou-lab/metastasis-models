# Adding in CLOW cluster sample labels to metadata

library(tidyverse)

out_dir <- "output/PAM50z_subtyping/"
dir.create(paste0(out_dir), recursive = T, showWarnings = F)

# Using the intrinsic cluster, look for expr of Cldn3, Cldn7, Elf3, choose a node
# correlation = 0.5000027: all samples are labeled as Claudin-low
CLOW_cluster_samples <- c("MB.0613", "MB.7078", "MB.0333", "MB.2845", "MB.5408", 
                          "MB.0644", "MB.3006", "MB.0285", "MB.0460", "MB.3707",
                          "MB.2963", "MB.0286", "MB.0155", "MB.7079", "MB.0442",
                          "MB.0367", "MB.0540", "MB.2724", "MB.0220", "MB.0476", 
                          "MB.5255", "MB.0489", "MB.0179", "GSM308364", "GSM36905",
                          "MB.7025", "MB.2821", "MB.0424", "MB.0617", "MB.0209", 
                          "MB.0893", "MB.0635", "MB.0432", "MB.0210", "MB.0502", 
                          "MB.2643", "MB.0289", "MB.0435", "MB.0499", "MB.0269", 
                          "MB.0876", "MB.5162", "MB.0403", "GSM50061", "GSM50066", 
                          "Sample.197", "GSM37050", "GSM50075", "GSM36889", 
                          "GSM308366", "GSM36794", "MB.0627", "GSM36816", "MB.2556",
                          "Sample.73", "Sample.350", "GSM36909", "MB.0885",
                          "MB.0153", "MB.0495", "MB.0168", "MB.7039", "MB.0211")

meta <- read.table("output/processing_metadata/n2486_combined_clinical_data.tsv", 
                   sep = "\t", header = T)

subtypes <- read.table("output/PAM50z_subtyping/n2486_Ensembl_dwd_PAM50z_output_pam50scores.txt",
                       sep = "\t", header = T, row.names = 1, check.names = F)

colnames(subtypes)[1:6] <- c("PAM50z_Basal", "PAM50z_Her2", "PAM50z_LumA", 
                             "PAM50z_LumB", "PAM50z_Normal", "PAM50z_Call")

subtypes <- subtypes %>% tibble::rownames_to_column(var = "Sample_ID")

meta <- meta %>%
  left_join(subtypes, by = "Sample_ID") %>% 
  mutate(CLOW_cluster = ifelse(Sample_ID %in% CLOW_cluster_samples, "Claudin-low", NA_character_)) %>% 
  mutate(PAM50z_CLOW_cluster = coalesce(CLOW_cluster, PAM50z_Call)) %>% 
  select(-ER_status, -PR_status, -HER2_status)

## Add clinical groups
receptor_data <- read.table("output/infer_ER_PR_HER2/n2486_Ensembl_dwd_inferred_clin_groups.txt", 
                            sep = "\t", header = T, check.names = F)

meta <- meta %>% 
  left_join(receptor_data, by = "Sample_ID") %>% 
  mutate(Clinical_Group_3 = case_when(
    (ER_status_w_inferred == "Positive" | PR_status_w_inferred == "Positive") & HER2_status_w_inferred == "Negative" ~ "HR+/HER2-", 
    ER_status_w_inferred == "Negative" & PR_status_w_inferred == "Negative" & HER2_status_w_inferred == "Negative" ~ "TNBC", 
    HER2_status_w_inferred == "Positive" ~ "HER2+")) %>% 
  select(Patient_ID, Patient_ID_2, Sample_ID, NKI_ID,
         Dataset, Dataset_5, Cohort, 
         MFS_Months, Relapse_Event, 
         Stage, Grade, Size, 
         ER_status, PR_status, HER2_status, ER_expr, PR_expr, HER2_expr, 
         inferred_ER_status, inferred_PR_status, inferred_HER2_status, 
         ER_status_w_inferred, PR_status_w_inferred, HER2_status_w_inferred, 
         IHC_subtype_w_inferred, Clinical_Group_3, 
         Bone_Relapse, Brain_Relapse, Liver_Relapse, Lung_Relapse,
         Four_site_recurrence, Relapse_Site, MultiMet_sites, 
         No_v_Distant, No_v_OnlyBone, No_v_OnlyBrain, No_v_OnlyLiver,
         No_v_OnlyLung, No_v_MultiMet, Single_v_MultiMet,
         PAM50z_Basal, PAM50z_Her2, PAM50z_LumA, PAM50z_LumB, PAM50z_Normal, 
         PAM50z_Call, Confidence, `ROR-S (Subtype Only)`, `ROR-S Group (Subtype Only)`, 
         `ROR-P (Subtype + Proliferation)`, `ROR-P Group (Subtype + Proliferation)`, 
         ER, Her2, CLOW_cluster, PAM50z_CLOW_cluster)

write.table(meta, paste0(out_dir, "/n2486_combined_meta_w_PAM50_CLOW_subtypes.txt"), 
            sep = "\t", row.names = F, col.names = T)
