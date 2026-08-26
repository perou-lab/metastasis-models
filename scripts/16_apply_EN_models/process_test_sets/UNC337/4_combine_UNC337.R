# Combine UNC337 with common Ensembl IDs

library(tidyverse)
library(GGally)
library(matrixStats)

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/"
dir.create(paste0(out_dir), recursive = T, showWarnings = F)

# Read in files for each of the microarray sets (mapped to ensembl ids) - remove GPL7504 & GPL6607
gpl_files <- c("GPL5325", "GPL1708", "GPL1390", "GPL887", "GPL885")

expr_list <- list()

for (file_name in gpl_files) {
  expr_list[[file_name]] <- read.table(paste0(out_dir, "/", file_name, "_mapped_ensembl_id_clpsd_filt_imp.txt"),
                                       sep = "\t", header = TRUE, check.names = FALSE)
}

# Identify common ENSEMBL IDs found in each dataset
common_ENSEMBL_IDs <- Reduce(intersect, lapply(expr_list, row.names))
paste0("There are ", length(common_ENSEMBL_IDs), " common ENSEMBL IDs across all 5 datasets")

# Subset each datasets to just the common ENSEMBL IDs
expr_list_subset <- lapply(expr_list, function(df) {
  df[match(common_ENSEMBL_IDs, rownames(df)), , drop = FALSE]
})

# Median center each dataset individually
expr_list_medcent <- lapply(expr_list_subset, function(df) {
  matrix <- as.matrix(df)
  medians <- rowMedians(matrix)
  matrix - medians
})

# Column standardize each dataset individually
expr_list_colstd <- lapply(expr_list_medcent, function(df) {
  scale(df)
})

# Combine all columns using matched ordering of rows (genes)
names(expr_list_colstd) <- NULL
combined_expr <- do.call(cbind, expr_list_colstd)

# Read in metadata
meta <- read.table(paste0(out_dir, "final/UNC337_metadata_from_GEO_combined.tsv"),
                   sep = "\t", header = T, check.names = F)
# table(meta$Sample_GEO_accession %in% colnames(combined))

meta <- meta %>% 
  filter(Sample_type == "Primary") %>% 
  filter(Dataset != "GPL6607") %>% 
  filter(Dataset != "GPL7504") %>% 
  filter(!is.na(OS_months))

combined_expr <- combined_expr[,colnames(combined_expr) %in% meta$Sample_GEO_accession]

# Plot PCA
pc <- prcomp(t(combined_expr))
pca <- as.data.frame(pc$x, make.names = F)
pca <- pca[meta$Sample_GEO_accession,]
meta <- meta %>% 
  mutate(PC1 = pca$PC1) %>% 
  mutate(PC2 = pca$PC2) %>% 
  mutate(PC3 = pca$PC3) %>% 
  mutate(PC4 = pca$PC4) %>% 
  mutate(PC5 = pca$PC5)

ggpairs(meta, 
        columns = which(colnames(meta) %in% paste0("PC", 1:5)),
        upper = list(continuous = "points"),
        legend = 2,
        aes(color = Dataset,
            alpha = 0.7)) +
  theme_minimal() +
  theme(legend.position="right") +
  scale_alpha_continuous(guide="none") 

ggpairs(meta, 
        columns = which(colnames(meta) %in% paste0("PC", 1:5)),
        upper = list(continuous = "points"),
        legend = 2,
        aes(color = Sample_type,
            alpha = 0.7)) +
  theme_minimal() +
  theme(legend.position="right") +
  scale_alpha_continuous(guide="none") 

ER_colors <- c("1" = "#2a2d7c", "0" = "darkred")
ggpairs(meta, 
        columns = which(colnames(meta) %in% paste0("PC", 1:5)),
        upper = list(continuous = "points"),
        legend = 2,
        aes(color = ER_status,
            alpha = 0.7)) +
  scale_color_manual(values = ER_colors, na.value = "grey70", name = NULL) +
  scale_fill_manual(values = ER_colors, na.value = "grey70") + 
  theme_minimal() +
  theme(legend.position="right") +
  scale_alpha_continuous(guide="none") 

PAM50_colors <- c("Basal" = "#c82127", 
                  "Claudin" = "goldenrod", 
                  "LumA" = "#2f56a6", 
                  "LumB" = "#9ccfef", 
                  "Normal" = "#4fb348", 
                  "Her2" = "#e38ebc")
ggpairs(meta, 
        columns = which(colnames(meta) %in% paste0("PC", 1:5)),
        upper = list(continuous = "points"),
        legend = 2,
        aes(color = PAM50_call,
            alpha = 0.7)) +
  scale_color_manual(values = PAM50_colors, name = NULL) +
  scale_fill_manual(values = PAM50_colors) + 
  theme_minimal() +
  theme(legend.position="right") +
  scale_alpha_continuous(guide="none") 

write.table(combined_expr, paste0(out_dir, "final/UNC337_n239_primaries_combined_ensembl_id.tsv"), 
            sep = "\t", col.names = T, row.names = T)

write.table(meta, paste0(out_dir, "final/UNC337_n239_metadata.tsv"), 
            sep = "\t", row.names = F, col.names = T)