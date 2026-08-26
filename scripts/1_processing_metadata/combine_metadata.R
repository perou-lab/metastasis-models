# Combine all metadata

library(tidyverse)

out_dir <- "output/processing_metadata/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Read in formatted metadata files
G2034_meta <- read.table("output/processing_metadata/GSE2034_clinical_data.tsv",
                       sep = "\t", header = T, check.names = F)
G2603_meta <- read.table("output/processing_metadata/GSE2603_clinical_data.tsv",
                         sep = "\t", header = T, check.names = F)
G12276_meta <- read.table("output/processing_metadata/GSE12276_clinical_data.tsv",
                         sep = "\t", header = T, check.names = F)
NKI_meta <- read.table("output/processing_metadata/NKI295_clinical_data.tsv",
                         sep = "\t", header = T, check.names = F)
MB_meta <- read.table("output/processing_metadata/METABRIC_clinical_data.tsv",
                      sep = "\t", header = T, check.names = F)

## Combine metadata
combined <- G2034_meta %>% 
  rbind(G2603_meta) %>% 
  rbind(G12276_meta) %>% 
  rbind(NKI_meta) %>% 
  rbind(MB_meta)

write.table(combined, paste0(out_dir, "n2486_combined_clinical_data.tsv"), 
            sep = "\t", col.names = T, row.names = F)
