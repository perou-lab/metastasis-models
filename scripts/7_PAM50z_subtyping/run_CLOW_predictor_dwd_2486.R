# Run Claudin-low predictor

out_dir <- "output/PAM50z_subtyping/"
dir.create(out_dir, showWarnings = F, recursive = T)

## Install genefu
# BiocManager::install(c("limma", "survcomp", "mclust"))
# BiocManager::install(c("AIMS", "iC10", "iC10TrainingData"))
# install.packages("/home/amich/genefu_2.40.0.tar.gz", repos = NULL, type = "source")

library(genefu)
library(tidyverse)
load("source_files/Breast-Cancer-Modules/data/gene_map_human.rda") # From Matthew

# Read in metadata
metadata <- read.table("output/infer_ER_PR_HER2/n2486_Ensembl_dwd_inferred_clin_groups.txt",
                       sep = "\t", header = T)

# Read in dwd-adjusted gene expression data (entrez IDs)
bc_expr <- read.table("output/PAM50z_subtyping/n2486_Ensembl_dwd_entrez.tsv",
                      sep = "\t", header = T, row.names = 1, check.names = F)
bc_expr <- bc_expr[,metadata$Sample_ID]
  
# Claudin-low predictor:

data(claudinLowData)

train <- claudinLowData
train$xd <- medianCtr(train$xd)

test <- medianCtr(bc_expr)

predout <- claudinLow(x=train$xd, classes=as.matrix(train$classes$Group,ncol=1), y=test)

saveRDS(predout, "output/PAM50z_subtyping/n2486_Ensembl_dwd_clow_predictor_out.rds")

predout <- readRDS("output/PAM50z_subtyping/n2486_Ensembl_dwd_clow_predictor_out.rds")

output <- predout$predictions

clow_samples <- output %>% filter(Call == "Claudin")

write.table(clow_samples, "output/PAM50z_subtyping/n2486_Ensembl_dwd_clow_predictor_out_samplelist.txt", 
            sep = '\t', row.names = F)
