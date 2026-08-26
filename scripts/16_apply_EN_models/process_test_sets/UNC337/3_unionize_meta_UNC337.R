# Process UNC-337 metadata
## (metadata was downloaded from GEO)

library(tidyverse)
library(readxl)

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/"
dir.create(paste0(out_dir, "final/"), recursive = T, showWarnings = F)


# GPL1708 ----------------------------------------------------------------------

GPL1708_meta <- read_excel("source_files/UNC337/GSE18229/extract_meta/GPL1708_transposed_meta.xlsx")

names(GPL1708_meta) <- make.unique(names(GPL1708_meta))

GPL1708_meta <- GPL1708_meta %>% 
  rename(Sample_GEO_accession = Sample_geo_accession, 
         Sample_name = Sample_title, 
         Age = Sample_characteristics_ch2...19, 
         ER_status = Sample_characteristics_ch2...20, 
         PR_status = Sample_characteristics_ch2...21,
         HER2_status = Sample_characteristics_ch2...22, 
         Node_status = Sample_characteristics_ch2...23, 
         Tumor_grade = Sample_characteristics_ch2...24, 
         Tumor_size = Sample_characteristics_ch2...25, 
         RFS_event = Sample_characteristics_ch2...26, 
         RFS_months = Sample_characteristics_ch2...27, 
         OS_event = Sample_characteristics_ch2...28, 
         OS_months = Sample_characteristics_ch2...29, 
         PAM50_call = Sample_characteristics_ch2...30) %>% 
  select(Sample_GEO_accession, Sample_name, Sample_description, Age, 
         ER_status, PR_status, HER2_status, Node_status, Tumor_grade, Tumor_size, 
         RFS_event, RFS_months, OS_event, OS_months, PAM50_call)

columns_to_clean <- c("Age", "ER_status", "PR_status", "HER2_status",
                      "Node_status", "Tumor_grade", "Tumor_size", "RFS_event",
                      "RFS_months", "OS_event", "OS_months", "PAM50_call")

GPL1708_meta <- GPL1708_meta %>% 
  mutate(across(all_of(columns_to_clean), 
                ~ str_trim(str_extract(.x, "(?<=: ).*$")))) %>% 
  mutate(across(everything(), ~na_if(.x, "--"))) %>% 
  mutate(Sample_type = case_when(
    Sample_description == "Breast Tumor" ~ "Primary", 
    str_detect(Sample_description, regex("Breast Tumor,", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("Normal breast", ignore_case = TRUE)) ~ "Normal Breast", 
    str_detect(Sample_description, regex("MMS", ignore_case = TRUE)) ~ "Mammospheres", 
    str_detect(Sample_name, regex("Met", ignore_case = TRUE)) ~ "Metastasis", 
    str_detect(Sample_description, regex("Primary", ignore_case = TRUE)) ~ "Primary"
  )) %>% 
  mutate(Dataset = "GPL1708")


# GPL1390 ----------------------------------------------------------------------

GPL1390_meta <- read_excel("source_files/UNC337/GSE18229/extract_meta/GPL1390_transposed_meta.xlsx")

names(GPL1390_meta) <- make.unique(names(GPL1390_meta))

GPL1390_meta <- GPL1390_meta %>% 
  rename(Sample_GEO_accession = Sample_geo_accession, 
         Sample_name = Sample_title, 
         Sample_description1 = Sample_source_name_ch2, 
         Age = Sample_characteristics_ch2...19, 
         ER_status = Sample_characteristics_ch2...20, 
         PR_status = Sample_characteristics_ch2...21,
         HER2_status = Sample_characteristics_ch2...22, 
         Node_status = Sample_characteristics_ch2...23, 
         Tumor_grade = Sample_characteristics_ch2...24, 
         Tumor_size = Sample_characteristics_ch2...25, 
         RFS_event = Sample_characteristics_ch2...26, 
         RFS_months = Sample_characteristics_ch2...27, 
         OS_event = Sample_characteristics_ch2...28, 
         OS_months = Sample_characteristics_ch2...29, 
         PAM50_call = Sample_characteristics_ch2...30) %>% 
  select(Sample_GEO_accession, Sample_name, Sample_description1, Age, 
         ER_status, PR_status, HER2_status, Node_status, Tumor_grade, Tumor_size, 
         RFS_event, RFS_months, OS_event, OS_months, PAM50_call) %>% 
  rename(Sample_description = Sample_description1)

columns_to_clean <- c("Age", "ER_status", "PR_status", "HER2_status",
                      "Node_status", "Tumor_grade", "Tumor_size", "RFS_event",
                      "RFS_months", "OS_event", "OS_months", "PAM50_call")

GPL1390_meta <- GPL1390_meta %>% 
  mutate(across(all_of(columns_to_clean), 
                ~ str_trim(str_extract(.x, "(?<=: ).*$")))) %>% 
  mutate(across(everything(), ~na_if(.x, "--"))) %>% 
  mutate(Sample_type = case_when(
    Sample_description == "Breast Tumor" ~ "Primary", 
    str_detect(Sample_description, regex("Breast Tumor,", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("Normal breast", ignore_case = TRUE)) ~ "Normal Breast", 
    str_detect(Sample_description, regex("MMS", ignore_case = TRUE)) ~ "Mammospheres", 
    str_detect(Sample_name, regex("Met", ignore_case = TRUE)) ~ "Metastasis", 
    str_detect(Sample_description, regex("Primary", ignore_case = TRUE)) ~ "Primary"
  )) %>% 
  mutate(Dataset = "GPL1390")

# GPL6607 ----------------------------------------------------------------------

GPL6607_meta <- read_excel("source_files/UNC337/GSE18229/extract_meta/GPL6607_transposed_meta.xlsx")

names(GPL6607_meta) <- make.unique(names(GPL6607_meta))

GPL6607_meta <- GPL6607_meta %>% 
  rename(Sample_GEO_accession = Sample_geo_accession, 
         Sample_name = Sample_title, 
         Age = Sample_characteristics_ch2...19, 
         ER_status = Sample_characteristics_ch2...20, 
         PR_status = Sample_characteristics_ch2...21,
         HER2_status = Sample_characteristics_ch2...22, 
         Node_status = Sample_characteristics_ch2...23, 
         Tumor_grade = Sample_characteristics_ch2...24, 
         Tumor_size = Sample_characteristics_ch2...25, 
         RFS_event = Sample_characteristics_ch2...26, 
         RFS_months = Sample_characteristics_ch2...27, 
         OS_event = Sample_characteristics_ch2...28, 
         OS_months = Sample_characteristics_ch2...29, 
         PAM50_call = Sample_characteristics_ch2...30) %>% 
  select(Sample_GEO_accession, Sample_name, Sample_description, Age, 
         ER_status, PR_status, HER2_status, Node_status, Tumor_grade, Tumor_size, 
         RFS_event, RFS_months, OS_event, OS_months, PAM50_call)

columns_to_clean <- c("Age", "ER_status", "PR_status", "HER2_status",
                      "Node_status", "Tumor_grade", "Tumor_size", "RFS_event",
                      "RFS_months", "OS_event", "OS_months", "PAM50_call")

GPL6607_meta <- GPL6607_meta %>% 
  mutate(across(all_of(columns_to_clean), 
                ~ str_trim(str_extract(.x, "(?<=: ).*$")))) %>% 
  mutate(across(everything(), ~na_if(.x, "--"))) %>% 
  mutate(Sample_type = case_when(
    Sample_description == "Breast Tumor" ~ "Primary", 
    str_detect(Sample_description, regex("Breast Tumor,", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("Normal breast", ignore_case = TRUE)) ~ "Normal Breast", 
    str_detect(Sample_description, regex("MMS", ignore_case = TRUE)) ~ "Mammospheres", 
    str_detect(Sample_name, regex("Met", ignore_case = TRUE)) ~ "Metastasis", 
    str_detect(Sample_description, regex("Primary", ignore_case = TRUE)) ~ "Primary"
  )) %>% 
  mutate(Dataset = "GPL6607")


# GPL7504 ----------------------------------------------------------------------

GPL7504_meta <- read_excel("source_files/UNC337/GSE18229/extract_meta/GPL7504_transposed_meta.xlsx")

names(GPL7504_meta) <- make.unique(names(GPL7504_meta))

GPL7504_meta <- GPL7504_meta %>% 
  rename(Sample_GEO_accession = Sample_geo_accession, 
         Sample_name = Sample_title) %>% 
  select(Sample_GEO_accession, Sample_name, Sample_description) %>% 
  mutate(Sample_type = case_when(
    str_detect(Sample_description, regex("MMS", ignore_case = TRUE)) ~ "Mammospheres",
    str_detect(Sample_description, regex("MCF7", ignore_case = TRUE)) ~ "Cell Line", 
    str_detect(Sample_description, regex("SUM149", ignore_case = TRUE)) ~ "Cell Line", 
    str_detect(Sample_description, regex("SUM159", ignore_case = TRUE)) ~ "Cell Line"
  )) %>% 
  mutate(Age = NA_character_, 
         ER_status = NA_character_, 
         PR_status = NA_character_, 
         HER2_status = NA_character_, 
         Node_status = NA_character_,
         Tumor_grade = NA_character_, 
         Tumor_size = NA_character_, 
         RFS_event = NA_character_, 
         RFS_months = NA_character_, 
         OS_event = NA_character_, 
         OS_months = NA_character_, 
         PAM50_call = NA_character_) %>% 
  mutate(Dataset = "GPL7504")


# GPL5325 ----------------------------------------------------------------------

GPL5325_meta <- read_excel("source_files/UNC337/GSE18229/extract_meta/GPL5325_transposed_meta.xlsx")

names(GPL5325_meta) <- make.unique(names(GPL5325_meta))

GPL5325_meta <- GPL5325_meta %>% 
  rename(Sample_GEO_accession = Sample_geo_accession, 
         Sample_name = Sample_title, 
         Age = Sample_characteristics_ch2...19, 
         ER_status = Sample_characteristics_ch2...20, 
         PR_status = Sample_characteristics_ch2...21,
         HER2_status = Sample_characteristics_ch2...22, 
         Node_status = Sample_characteristics_ch2...23, 
         Tumor_grade = Sample_characteristics_ch2...24, 
         Tumor_size = Sample_characteristics_ch2...25, 
         RFS_event = Sample_characteristics_ch2...26, 
         RFS_months = Sample_characteristics_ch2...27, 
         OS_event = Sample_characteristics_ch2...28, 
         OS_months = Sample_characteristics_ch2...29, 
         PAM50_call = Sample_characteristics_ch2...30) %>% 
  select(Sample_GEO_accession, Sample_name, Sample_description, Age, 
         ER_status, PR_status, HER2_status, Node_status, Tumor_grade, Tumor_size, 
         RFS_event, RFS_months, OS_event, OS_months, PAM50_call)

columns_to_clean <- c("Age", "ER_status", "PR_status", "HER2_status",
                      "Node_status", "Tumor_grade", "Tumor_size", "RFS_event",
                      "RFS_months", "OS_event", "OS_months", "PAM50_call")

GPL5325_meta <- GPL5325_meta %>% 
  mutate(across(all_of(columns_to_clean), 
                ~ str_trim(str_extract(.x, "(?<=: ).*$")))) %>% 
  mutate(across(everything(), ~na_if(.x, "--"))) %>% 
  mutate(Sample_type = case_when(
    Sample_description == "Breast Tumor" ~ "Primary", 
    str_detect(Sample_description, regex("Breast Tumor,", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("Normal breast", ignore_case = TRUE)) ~ "Normal Breast", 
    str_detect(Sample_description, regex("MMS", ignore_case = TRUE)) ~ "Mammospheres", 
    str_detect(Sample_description, regex("Primary", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("Metaplastic breast carcinoma", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("UNC0", ignore_case = TRUE)) ~ "Primary"
  )) %>% 
  # GSM455555 columns did not align correctly
  mutate(ER_status = ifelse(Sample_GEO_accession == "GSM455555", NA, ER_status)) %>% 
  mutate(Age = ifelse(Sample_GEO_accession == "GSM455555", NA, Age)) %>% 
  mutate(PAM50_call = ifelse(Sample_GEO_accession == "GSM455555", "Claudin", PAM50_call)) %>% 
  mutate(Dataset = "GPL5325")


# GPL887 -----------------------------------------------------------------------

GPL887_meta <- read_excel("source_files/UNC337/GSE18229/extract_meta/GPL887_transposed_meta.xlsx")

names(GPL887_meta) <- make.unique(names(GPL887_meta))

GPL887_meta <- GPL887_meta %>% 
  rename(Sample_GEO_accession = Sample_geo_accession, 
         Sample_name = Sample_title, 
         Sample_description = Sample_source_name_ch2,
         Age = Sample_characteristics_ch2...19, 
         ER_status = Sample_characteristics_ch2...20, 
         PR_status = Sample_characteristics_ch2...21,
         HER2_status = Sample_characteristics_ch2...22, 
         Node_status = Sample_characteristics_ch2...23, 
         Tumor_grade = Sample_characteristics_ch2...24, 
         Tumor_size = Sample_characteristics_ch2...25, 
         RFS_event = Sample_characteristics_ch2...26, 
         RFS_months = Sample_characteristics_ch2...27, 
         OS_event = Sample_characteristics_ch2...28, 
         OS_months = Sample_characteristics_ch2...29, 
         PAM50_call = Sample_characteristics_ch2...30) %>% 
  select(Sample_GEO_accession, Sample_name, Sample_description, Age, 
         ER_status, PR_status, HER2_status, Node_status, Tumor_grade, Tumor_size, 
         RFS_event, RFS_months, OS_event, OS_months, PAM50_call)

columns_to_clean <- c("Age", "ER_status", "PR_status", "HER2_status",
                      "Node_status", "Tumor_grade", "Tumor_size", "RFS_event",
                      "RFS_months", "OS_event", "OS_months", "PAM50_call")

GPL887_meta <- GPL887_meta %>% 
  mutate(across(all_of(columns_to_clean), 
                ~ str_trim(str_extract(.x, "(?<=: ).*$")))) %>% 
  mutate(across(everything(), ~na_if(.x, "--"))) %>% 
  mutate(Sample_type = case_when(
    Sample_description == "Breast Tumor" ~ "Primary", 
    str_detect(Sample_description, regex("Breast Tumor,", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("Normal breast", ignore_case = TRUE)) ~ "Normal Breast", 
    str_detect(Sample_description, regex("MMS", ignore_case = TRUE)) ~ "Mammospheres", 
    str_detect(Sample_name, regex("Met", ignore_case = TRUE)) ~ "Metastasis", 
    str_detect(Sample_description, regex("Lung Met", ignore_case = TRUE)) ~ "Metastasis", 
    str_detect(Sample_description, regex("Primary", ignore_case = TRUE)) ~ "Primary"
  )) %>% 
  mutate(Dataset = "GPL887")


# GPL885 -----------------------------------------------------------------------

GPL885_meta <- read_excel("source_files/UNC337/GSE18229/extract_meta/GPL885_transposed_meta.xlsx")

names(GPL885_meta) <- make.unique(names(GPL885_meta))

GPL885_meta <- GPL885_meta %>% 
  rename(Sample_GEO_accession = Sample_geo_accession, 
         Sample_name = Sample_title, 
         Sample_description1 = Sample_source_name_ch2,
         Age = Sample_characteristics_ch2...19, 
         ER_status = Sample_characteristics_ch2...20, 
         PR_status = Sample_characteristics_ch2...21,
         HER2_status = Sample_characteristics_ch2...22, 
         Node_status = Sample_characteristics_ch2...23, 
         Tumor_grade = Sample_characteristics_ch2...24, 
         Tumor_size = Sample_characteristics_ch2...25, 
         RFS_event = Sample_characteristics_ch2...26, 
         RFS_months = Sample_characteristics_ch2...27, 
         OS_event = Sample_characteristics_ch2...28, 
         OS_months = Sample_characteristics_ch2...29, 
         PAM50_call = Sample_characteristics_ch2...30) %>% 
  select(Sample_GEO_accession, Sample_name, Sample_description1, Age, 
         ER_status, PR_status, HER2_status, Node_status, Tumor_grade, Tumor_size, 
         RFS_event, RFS_months, OS_event, OS_months, PAM50_call) %>% 
  rename(Sample_description = Sample_description1)

columns_to_clean <- c("Age", "ER_status", "PR_status", "HER2_status",
                      "Node_status", "Tumor_grade", "Tumor_size", "RFS_event",
                      "RFS_months", "OS_event", "OS_months", "PAM50_call")

GPL885_meta <- GPL885_meta %>% 
  mutate(across(all_of(columns_to_clean), 
                ~ str_trim(str_extract(.x, "(?<=: ).*$")))) %>% 
  mutate(across(everything(), ~na_if(.x, "--"))) %>% 
  mutate(Sample_type = case_when(
    Sample_description == "Breast Tumor" ~ "Primary", 
    str_detect(Sample_description, regex("Breast Tumor,", ignore_case = TRUE)) ~ "Primary",
    str_detect(Sample_description, regex("Normal breast", ignore_case = TRUE)) ~ "Normal Breast", 
    str_detect(Sample_description, regex("MMS", ignore_case = TRUE)) ~ "Mammospheres", 
    str_detect(Sample_name, regex("Met", ignore_case = TRUE)) ~ "Metastasis", 
    str_detect(Sample_description, regex("Lung Met", ignore_case = TRUE)) ~ "Metastasis", 
    str_detect(Sample_description, regex("Primary", ignore_case = TRUE)) ~ "Primary", 
    str_detect(Sample_description, regex("recurrence", ignore_case = TRUE)) ~ "Primary"
  )) %>% 
  mutate(Dataset = "GPL885")


# Combine metadata -------------------------------------------------------------
combined_meta <- bind_rows(GPL1390_meta, GPL1708_meta, GPL5325_meta,
                           GPL6607_meta, GPL7504_meta, GPL885_meta, GPL887_meta)

write.table(combined_meta, paste0(out_dir, "final/UNC337_metadata_from_GEO_combined.tsv"), 
            sep = "\t", col.names = T, row.names = F)