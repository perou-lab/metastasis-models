## PAM50 subtyping for n5 combined set (dwd, ensembl, n = 2486)

library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(ctc)
library(impute)
library(matrixStats)
library(circlize)
library(heatmap.plus)
library(org.Hs.eg.db)
library(annotate)

set.seed(123)

out_dir <- "output/PAM50z_subtyping/"
dir.create(out_dir, showWarnings = F, recursive = T)

# Define colors ----------------------------------------------------------------

pam50_colors <- c("Basal" = "#c91508", 
									"Her2" = "#ed8ccb",
									"LumA" = "#084dc4",
									"LumB" = "#9cd0f0",
									"Normal" = "#4fb30c")

# Load data --------------------------------------------------------------------

# Read in metadata
metadata <- read.table("output/infer_ER_PR_HER2/n2486_Ensembl_dwd_inferred_clin_groups.txt", 
                       sep = "\t", header = T)

# Read in dwd-adjusted gene expression data
bc_expr <- read.table("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols.tsv",
											sep = "\t", header = T, row.names = 1, check.names = F)
bc_expr <- bc_expr[,metadata$Sample_ID]

# Read in intrinsic gene expr for clustering
intrinsic_expr <- read.table("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols_intrinsic_genes.txt",
														 header = T, sep = "\t", row.names = 1, check.names = F)

# Reorder expression and metadata to match order of intrinsic cluster
bc_expr <- bc_expr[, colnames(intrinsic_expr)]
metadata <- metadata[match(colnames(intrinsic_expr), metadata$Sample_ID), ]

# Convert the names of three genes in the RNA-seq data to match the  
# corresponding PAM50 gene name
row.names(bc_expr)[row.names(bc_expr) == "ORC6"] <- "ORC6L"
row.names(bc_expr)[row.names(bc_expr) == "NDC80"] <- "KNTC2"
row.names(bc_expr)[row.names(bc_expr) == "NUF2"] <- "CDCA1"

# Read in table of trained sigma values on PAM50 genes needed for Zhao normalization
sigmas <- read.table("../metastatic_site_EN/source_files/PAM50/SIGMAS_Zhao.txt",
										 sep = "\t", row.names = 1, header = T)

sigmas <- sigmas[row.names(sigmas) %in% row.names(bc_expr),]

# Subset expression matrix to only include PAM50 genes
pam50_expr <- bc_expr[row.names(sigmas),]

# Determine the Zhao class of each sample
metadata <- metadata %>% 
  mutate(Zhao_class = case_when(
    ER_status_w_inferred == "Positive" & HER2_status_w_inferred == "Positive" ~ "HER2pos_ERpos", 
    ER_status_w_inferred == "Positive" & HER2_status_w_inferred == "Negative" ~ "ERpos_HER2neg", 
    ER_status_w_inferred == "Negative" & HER2_status_w_inferred == "Positive" ~ "HER2pos_ERneg", 
    ER_status_w_inferred == "Negative" & HER2_status_w_inferred == "Negative" ~ "ERneg_HER2neg"
  ))

fn <- paste0(out_dir, "/n2486_Ensembl_dwd_inferred_Zhao_classes.txt")
write_tsv(metadata %>% dplyr::select(Sample_ID, Zhao_class), fn)

# Zhao subgroup-specific PAM50 gene centering ----------------------------------

# Function for quantile centering based on the quantile and group
quantile_centering <- function(expr_matrix, gene_quantile) {
	q_values <- diag(apply(expr_matrix, 1, function(x) quantile(x, percentile_group)))
	return(expr_matrix - q_values)
}

# Make sigma matrix PAM50 gene order match gene expression matrix
sigmas <- sigmas[row.names(pam50_expr),]

# Separately center the gene expression for each class
centered_pam50_expr <- pam50_expr
for (class in unique(metadata$Zhao_class)) {
	samples <- metadata$Sample_ID[metadata$Zhao_class == class]
	class_pam50_expr <- pam50_expr[,samples]
	percentile_group <- sigmas[,class]
	names(percentile_group) <- row.names(sigmas)
	centered_pam50_expr[,samples] <- quantile_centering(class_pam50_expr, percentile_group)
}

centered_pam50_df <- cbind(Sample_ID = row.names(centered_pam50_expr),
													 data.frame(centered_pam50_expr, check.names = F))

fn <- paste0(out_dir, "/n2486_Ensembl_dwd_PAM50genes_Zhao_normalized.txt")
write.table(centered_pam50_df, fn, sep = "\t", quote = F, row.names = F)

# PAM50 bioclassifier ----------------------------------------------------------

paramDir <- "../metastatic_site_EN/source_files/PAM50/bioclassifier_R"  # point to the bioclassifier_R folder
inputDir <- out_dir  # the folder where your gene expression array is located
inputFile <- "n2486_Ensembl_dwd_PAM50genes_Zhao_normalized.txt" # your input file with the 50 PAM50 genes
short <- "n2486_Ensembl_dwd_PAM50z_output"  # output name for files
calibrationParameters <- -1 # don't calibrate data (already Zhao centered)
hasClinical <- FALSE
collapseMethod <- "mean" # can be mean or iqr (probe with max iqr is selected)
source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))

# Add PAM50 results to metadata
fn <- paste0(out_dir, "/n2486_Ensembl_dwd_PAM50z_output_pam50scores.txt")
pam50_new <- read.table(fn, sep = "\t", header = T) %>% 
	dplyr::rename(Sample_ID = X)
colnames(pam50_new)[-1] <- paste0("PAM50z_", colnames(pam50_new)[-1])

compare <- inner_join(metadata, pam50_new)

# Intrinsic gene clustering, with labeled Zhao4 PAM50 calls:

row_clust <- xcluster2r("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols_intrinsic_genes.gtr",
												distance = "pearson")
col_clust <- xcluster2r("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols_intrinsic_genes.atr",
												distance = "pearson")

col_annot <- HeatmapAnnotation(PAM50z = compare$PAM50z_Call,
															 col = list(PAM50z= pam50_colors))

intrinsic_expr <- intrinsic_expr[,compare$Sample_ID]
p <- Heatmap(intrinsic_expr,
						 top_annotation = col_annot,
						 show_column_names = FALSE,
						 show_row_names = FALSE,
						 name = "expr",
						 column_title = "n2486 dwd combined PAM50z calls", 
						 row_title = "Intrinsic genes",
						 cluster_rows = row_clust,
						 cluster_columns = col_clust,
						 border = TRUE,
						 use_raster = FALSE, 
						 col = colorRamp2(c(-5, 0, 5), 
						 								 c("#1BB7E5", "black", "#FEFF00")))

fn <- paste0(out_dir, "/n2486_Ensembl_dwd_symbols_intrinsic_genes_heatmap_PAM50.pdf")
pdf(fn, width = 10, height = 7)
print(p)
dev.off()
