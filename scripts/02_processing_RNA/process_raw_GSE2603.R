# Processing raw GSE2603 data

library(tidyverse)
library(impute)
library(annotate)
library("hgu133a.db")
source("source_files/functions.R")
set.seed(123)

out_dir <- "output/processing_RNA/GSE2603/"
dir.create(out_dir, showWarnings = F, recursive = T)

## GSE2603 data was obtained from GEO on 8/28/2025 at the following link:
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2603/matrix/
  
G2603 <- read.table("source_files/GSE2603_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!", row.names = 1)

# Convert probes to ENSEMBL IDs
affy_id_convert <- as.data.frame(row.names(G2603))
colnames(affy_id_convert) <- c("affy_id")
affy_id_convert <- AnnotationDbi::select(hgu133a.db, 
                                         keys = keys(hgu133a.db, keytype = "PROBEID"), 
                                         columns = "ENSEMBL", 
                                         keytype = "PROBEID")

G2603 <- tibble::rownames_to_column(as.data.frame(G2603), var = "PROBEID")
G2603 <- merge(affy_id_convert, G2603, by = "PROBEID")
G2603 <- G2603 %>% 
  drop_na(ENSEMBL) # remove all probes missing ENSEMBL mapping
G2603 <- G2603[,-which(colnames(G2603) == "PROBEID")]
write.table(G2603, paste0(out_dir, "/G2603_Ensembl.txt"), 
            sep = "\t", row.names = F, col.names = T)

# Collapse based on duplicated ENSEMBL IDs
file <- paste0(out_dir, "/G2603_Ensembl.txt")
x <- readarray(file, hr=1, impute=F)
write.table(signif(x$xd,4), paste0(out_dir, "/G2603_Ensembl_clpsd.txt"), 
            sep = "\t", col.names = NA)

G2603 <- read.table(paste0(out_dir, "/G2603_Ensembl_clpsd.txt"), 
                    sep = "\t", header = T, row.names = 1)

# Filter probes that are missing 30% of data: 
G2603 <- G2603[rowSums(!is.na(G2603)) >= 0.7*ncol(G2603),]

# Impute missing data (k=10, default)
G2603 <- impute.knn(as.matrix(G2603))
G2603 <- G2603$data

# Save final procesed matrix with ENSEMBL:
write.table(G2603, paste0(out_dir, "/G2603_Ensembl_clpsd_filt_imp.txt"), 
            sep = "\t", row.names = T, col.names = T)