# Processing raw GSE12276 data

library(tidyverse)
library(impute)
library(annotate)
library(hgu133plus2.db)
source("source_files/functions.R")
set.seed(123)

out_dir <- "output/processing_RNA/GSE12276/"
dir.create(out_dir, showWarnings = F, recursive = T)

## GSE12276 data was obtained from GEO on 8/28/2025 at the following link:
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12276/matrix/
  
G12276 <- read.table("source_files/GSE12276_series_matrix.txt", 
                    sep = "\t", header = T, comment.char = "!", row.names = 1)

# Convert probes to ENSEMBL IDs
affy_id_convert <- as.data.frame(row.names(G12276))
colnames(affy_id_convert) <- c("affy_id")
affy_id_convert <- AnnotationDbi::select(hgu133plus2.db, 
                                         keys = keys(hgu133plus2.db, keytype = "PROBEID"), 
                                         columns = "ENSEMBL", 
                                         keytype = "PROBEID")

G12276 <- tibble::rownames_to_column(as.data.frame(G12276), var = "PROBEID")
G12276 <- merge(affy_id_convert, G12276, by = "PROBEID")
G12276 <- G12276 %>% 
  drop_na(ENSEMBL) # remove all probes missing ENSEMBL mapping
G12276 <- G12276[,-which(colnames(G12276) == "PROBEID")]
write.table(G12276, paste0(out_dir, "/G12276_Ensembl.txt"), 
            sep = "\t", row.names = F, col.names = T)

# Collapse based on duplicated ENSEMBL IDs
file <- paste0(out_dir, "/G12276_Ensembl.txt")
x <- readarray(file, hr=1, impute=F)
write.table(signif(x$xd,4), paste0(out_dir, "/G12276_Ensembl_clpsd.txt"), 
            sep = "\t", col.names = NA)

G12276 <- read.table(paste0(out_dir, "/G12276_Ensembl_clpsd.txt"), 
                     sep = "\t", header = T, row.names = 1)

# Filter probes that are missing 30% of data:
G12276 <- G12276[rowSums(!is.na(G12276)) >= 0.7*ncol(G12276),]

# Impute missing data (k=10, default)
G12276 <- impute.knn(as.matrix(G12276))
G12276 <- G12276$data

# Save final processed matrix with ENSEMBL:
write.table(G12276, paste0(out_dir, "/G12276_Ensembl_clpsd_filt_imp.txt"), 
            sep = "\t", row.names = T, col.names = T)
