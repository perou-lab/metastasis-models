# Common ENSEMBL IDs across all datasets

library(tidyverse)
library(matrixStats)

out_dir <- "output/processing_RNA/"

G2034 <- read.table("output/processing_RNA/GSE2034/G2034_Ensembl_clpsd_filt_imp.txt", 
                    sep = "\t", header = T, row.names = 1)
G2603 <- read.table("output/processing_RNA/GSE2603/G2603_Ensembl_clpsd_filt_imp.txt", 
                    sep = "\t", header = T, row.names = 1)
G12276 <- read.table("output/processing_RNA/GSE12276/G12276_Ensembl_clpsd_filt_imp.txt", 
                     sep = "\t", header = T, row.names = 1)
NKI <- read.table("output/processing_RNA/NKI295/NKI295_Ensembl_clpsd_filt_imp.txt", 
                  sep = "\t", header = T, row.names = 1)
MB <- read.table("output/processing_RNA/METABRIC/METABRIC_Ensembl_clpsd_filt_imp.txt", 
                 sep = "\t", header = T, row.names = 1)

# Identify common ENSEMBL IDs found in each dataset
list_ENSEMBL_IDs <- list(G2034, G2603, G12276, NKI, MB)
common_ENSEMBL_IDs <- Reduce(intersect, lapply(list_ENSEMBL_IDs, row.names))
paste0("There are ", length(common_ENSEMBL_IDs), " common ENSEMBL IDs across all 5 datasets")

# Read in metadata and filter for samples with defined metastatic sites: 
## bone only, brain only, lung only, liver only, multi-met, no relapse
metadata <- read.table(file = "output/processing_metadata/n2486_combined_clinical_data.tsv",
                       sep = "\t", header = T) %>% 
  filter(Relapse_Site %in% c("Bone", "Brain", "Liver", "Lung", "Multimet", "No_Relapse"))

# Subset each datasets to just the common ENSEMBL IDs
G2034 <- G2034[common_ENSEMBL_IDs,colnames(G2034) %in% metadata$Sample_ID]
G2603 <- G2603[common_ENSEMBL_IDs,colnames(G2603) %in% metadata$Sample_ID]
G12276 <- G12276[common_ENSEMBL_IDs,colnames(G12276) %in% metadata$Sample_ID]
NKI <- NKI[common_ENSEMBL_IDs,colnames(NKI) %in% metadata$Sample_ID]
MB <- MB[common_ENSEMBL_IDs,colnames(MB) %in% metadata$Sample_ID]

# Log2 transformation
G2034 <- as.matrix(log2(G2034))
G2603 <- as.matrix(log2(G2603))
G12276 <- as.matrix(log2(G12276))
NKI <- as.matrix(log2(10^NKI)) # Original NKI data is in log10 ratios
MB <- as.matrix(MB) # dataset is already log2 transformed

# Median center
G2034 <- G2034 - rowMedians(G2034)
G2603 <- G2603 - rowMedians(G2603)
G12276 <- G12276 - rowMedians(G12276)
NKI <- NKI - rowMedians(NKI)
MB <- MB - rowMedians(MB)

# Column standardization
G2034 <- scale(G2034)
G2603 <- scale(G2603)
G12276 <- scale(G12276)
NKI <- scale(NKI)
MB <- scale(MB)

# Save each dataset for DWD combining
write.table(G2034, paste0(out_dir, "/GSE2034/G2034_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt"), 
            sep = "\t")
write.table(G2603, paste0(out_dir, "/GSE2603/G2603_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt"), 
            sep = "\t")
write.table(G12276, paste0(out_dir, "/GSE12276/G12276_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt"), 
            sep = "\t")
write.table(NKI, paste0(out_dir, "/NKI295/NKI295_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt"), 
            sep = "\t")
write.table(MB, paste0(out_dir, "/METABRIC/METABRIC_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt"), 
            sep = "\t")

