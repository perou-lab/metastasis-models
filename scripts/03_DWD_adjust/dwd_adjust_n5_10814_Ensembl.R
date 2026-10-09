# Perform DWD adjustment in MatLab to generate dirvec files
## METABRIC is used as the reference
## All datasets are adjusted to METABRIC; METABRIC expr is not changed

# Using MatLab DWD output, adjust the datasets using the DWD direction vector. 
# Save adjusted files for combining and future analysis

out_dir <- "output/DWD_adjust/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

# Load datasets ----------------------------------------------------------------
G2034 <- read.table("output/processing_RNA/GSE2034/G2034_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt", 
                    sep = "\t", header = T, row.names = 1)
G2603 <- read.table("output/processing_RNA/GSE2603/G2603_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt", 
                    sep = "\t", header = T, row.names = 1)
G12276 <- read.table("output/processing_RNA/GSE12276/G12276_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt", 
                     sep = "\t", header = T, row.names = 1)
NKI <- read.table("output/processing_RNA/NKI295/NKI295_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt", 
                  sep = "\t", header = T, row.names = 1)
MB <- read.table("output/processing_RNA/METABRIC/METABRIC_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt", 
                 sep = "\t", header = T, row.names = 1)

# Load DWD direction vectors ---------------------------------------------------
MB_G2034_vec <- read.table("output/DWD_adjust/dirvec_METABRIC_G2034.txt", 
                           header = F)[,1]
MB_G2603_vec <- read.table("output/DWD_adjust/dirvec_METABRIC_G2603.txt", 
                           header = F)[,1]
MB_G12276_vec <- read.table("output/DWD_adjust/dirvec_METABRIC_G12276.txt", 
                            header = F)[,1]
MB_NKI_vec <- read.table("output/DWD_adjust/dirvec_METABRIC_NKI295.txt", 
                         header = F)[,1]

# METABRIC / G12276 ------------------------------------------------------------
# project data (two classes) of n1 vectors and n2 vectors 

METABRIC_p <- t(MB) %*% as.vector(MB_G12276_vec) # original set projection
G12276_p <- t(G12276) %*% as.vector(MB_G12276_vec) # new set projection

# calculate class means of projection (two classes) in R^d

METABRIC_pm <- mean(METABRIC_p) * MB_G12276_vec 
G12276_pm <- mean(G12276_p) * MB_G12276_vec 

# compute the difference (the starting point to the ending point)

METABRIC_G12276_delta <- METABRIC_pm - G12276_pm

# shift new dataset to the original set 

G12276_adj <- apply(G12276, 2, function(x) { x + METABRIC_G12276_delta })

# METABRIC / G2034 ------------------------------------------------------------
# project data (two classes) of n1 vectors and n2 vectors 

METABRIC_p <- t(MB) %*% as.vector(MB_G2034_vec) # original set projection
G2034_p <- t(G2034) %*% as.vector(MB_G2034_vec) # new set projection

# calculate class means of projection (two classes) in R^d

METABRIC_pm <- mean(METABRIC_p) * MB_G2034_vec 
G2034_pm <- mean(G2034_p) * MB_G2034_vec 

# compute the difference (the starting point to the ending point)

METABRIC_G2034_delta <- METABRIC_pm - G2034_pm

# shift new dataset to the original set 

G2034_adj <- apply(G2034, 2, function(x) { x + METABRIC_G2034_delta })

# METABRIC / G2603 ------------------------------------------------------------
# project data (two classes) of n1 vectors and n2 vectors 

METABRIC_p <- t(MB) %*% as.vector(MB_G2603_vec) # original set projection
G2603_p <- t(G2603) %*% as.vector(MB_G2603_vec) # new set projection

# calculate class means of projection (two classes) in R^d

METABRIC_pm <- mean(METABRIC_p) * MB_G2603_vec 
G2603_pm <- mean(G2603_p) * MB_G2603_vec 

# compute the difference (the starting point to the ending point)

METABRIC_G2603_delta <- METABRIC_pm - G2603_pm

# shift new dataset to the original set 

G2603_adj <- apply(G2603, 2, function(x) { x + METABRIC_G2603_delta })

# METABRIC / NKI295 ------------------------------------------------------------
# project data (two classes) of n1 vectors and n2 vectors 

METABRIC_p <- t(MB) %*% as.vector(MB_NKI_vec) # original set projection
NKI295_p <- t(NKI) %*% as.vector(MB_NKI_vec) # new set projection

# calculate class means of projection (two classes) in R^d

METABRIC_pm <- mean(METABRIC_p) * MB_NKI_vec 
NKI295_pm <- mean(NKI295_p) * MB_NKI_vec 

# compute the difference (the starting point to the ending point)

METABRIC_NKI295_delta <- METABRIC_pm - NKI295_pm

# shift new dataset to the original set 

NKI295_adj <- apply(NKI, 2, function(x) { x + METABRIC_NKI295_delta })

# Save adjusted datasets -------------------------------------------------------
write.table(G12276_adj, paste0(out_dir, 'G12276_DWDadj-to-MB.txt'), sep = '\t')
write.table(G2034_adj, paste0(out_dir, 'G2034_DWDadj-to-MB.txt'), sep = '\t')
write.table(G2603_adj, paste0(out_dir, 'G2603_DWDadj-to-MB.txt'), sep = '\t')
write.table(NKI295_adj, paste0(out_dir, 'NKI295_DWDadj-to-MB.txt'), sep = '\t')

# Combine datasets -------------------------------------------------------------
MB <- read.table("output/processing_RNA/METABRIC/METABRIC_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt", 
                 sep = "\t", header = T, row.names = 1)
G12276 <- read.table("output/DWD_adjust/G12276_DWDadj-to-MB.txt",
                     sep = "\t", header = T, row.names = 1)
G2603 <- read.table("output/DWD_adjust/G2603_DWDadj-to-MB.txt",
                    sep = "\t", header = T, row.names = 1)
G2034 <- read.table("output/DWD_adjust/G2034_DWDadj-to-MB.txt",
                    sep = "\t", header = T, row.names = 1)
NKI <- read.table("output/DWD_adjust/NKI295_DWDadj-to-MB.txt", 
                  sep = "\t", header = T, row.names = 1)
dwd <- cbind(MB, G12276, G2603, G2034, NKI)
write.table(dwd, paste0(out_dir, "n2486_Ensembl_dwd.txt"), 
            sep = "\t", row.names = T, col.names = T)
