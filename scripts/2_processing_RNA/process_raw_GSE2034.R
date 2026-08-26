# Processing raw GSE2034 data

library(tidyverse)
library(impute)
library(annotate)
library("hgu133a.db")
source("source_files/functions.R")
set.seed(123)

out_dir <- "output/processing_RNA/GSE2034/"
dir.create(out_dir, showWarnings = F, recursive = T)

## GSE2034 data was obtained from GEO on 8/28/2025 at the following link:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034
  
G2034 <- read.table("source_files/GSE2034_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!", row.names = 1)

# Convert probes to ENSEMBL IDs
affy_id_convert <- as.data.frame(row.names(G2034))
colnames(affy_id_convert) <- c("affy_id")
affy_id_convert <- AnnotationDbi::select(hgu133a.db, 
                                         keys = keys(hgu133a.db, keytype = "PROBEID"), 
                                         columns = "ENSEMBL", 
                                         keytype = "PROBEID")

G2034 <- tibble::rownames_to_column(as.data.frame(G2034), var = "PROBEID")
G2034 <- merge(affy_id_convert, G2034, by = "PROBEID")
G2034 <- G2034 %>% 
  drop_na(ENSEMBL) # remove all probes missing ENSEMBL mapping
G2034 <- G2034[,-which(colnames(G2034) == "PROBEID")]
write.table(G2034, paste0(out_dir, "/G2034_Ensembl.txt"), 
            sep = "\t", row.names = F, col.names = T)

# Collapse based on duplicated ENSEMBL IDs
file <- paste0(out_dir, "/G2034_Ensembl.txt")
x <- readarray(file, hr=1, impute=F)
write.table(signif(x$xd,4), paste0(out_dir, "/G2034_Ensembl_clpsd.txt"), 
            sep = "\t", col.names = NA)

G2034 <- read.table(paste0(out_dir, "/G2034_Ensembl_clpsd.txt"), 
                    sep = "\t", header = T, row.names = 1)

# Filter probes that are missing 30% of data: 
G2034 <- G2034[rowSums(!is.na(G2034)) >= 0.7*ncol(G2034),]

# Impute missing data (k=10, default)
G2034 <- impute.knn(as.matrix(G2034))
G2034 <- G2034$data

# Save final processed matrix with ENSEMBL:
write.table(G2034, paste0(out_dir, "/G2034_Ensembl_clpsd_filt_imp.txt"), 
            sep = "\t", row.names = T, col.names = T)