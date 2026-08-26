# plot intrinsic cluster using cluster files from Cluster 3.0

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

out_dir <- "output/apply_EN_models/process_test_sets/UNC337/intrinsic_cluster/"

cohort_colors <- c("GPL1390" = "#082a54",
                   "GPL1708" = "#59a89c", 
                   "GPL5325" = "#a559aa",
                   "GPL887" = "#cecece", 
                   "GPL885" = "#e02b35")

clin_colors_3 <- c("TNBC" = "darkred", 
                   "HER2+" = "#CC79A7", 
                   "HR+/HER2-" = "darkblue")

pam50_colors <- c("Basal" = "#c91508", 
                  "Her2" = "#ed8ccb",
                  "LumA" = "#084dc4",
                  "LumB" = "#9cd0f0",
                  "Normal" = "#4fb30c", 
                  "Claudin" = "goldenrod")

## Metadata 
metadata <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_metadata.tsv",
                       sep = "\t", header = T, check.names = F)

## Intrinsic Clustering
heatmap_data <- read.table("output/apply_EN_models/process_test_sets/UNC337/intrinsic_cluster/UNC337_n239_primaries_symbols_intrinsic_genes.txt",
                           sep = "\t", header = T, row.names = 1, check.names = F)
heatmap_data <- as.matrix(heatmap_data)

row_clust <- xcluster2r("output/apply_EN_models/process_test_sets/UNC337/intrinsic_cluster/UNC337_n239_primaries_symbols_intrinsic_genes.gtr", 
                        distance = "pearson")
col_clust <- xcluster2r("output/apply_EN_models/process_test_sets/UNC337/intrinsic_cluster/UNC337_n239_primaries_symbols_intrinsic_genes.atr", 
                          distance = "pearson")

heatmap_metadata <- metadata[match(colnames(heatmap_data), metadata$Sample_GEO_accession),]

## Heatmap Annotations
top_annot <- HeatmapAnnotation(Study = heatmap_metadata$Dataset, 
                               PAM50z_Subtype = heatmap_metadata$PAM50_call,
                               col = list(PAM50z_Subtype = pam50_colors,
                                          Study = cohort_colors),
                               show_legend = c(T, T),
                               annotation_legend_param = list( 
                                 PAM50z_Subtype = list(title = "PAM50 Subtype")))

heatmap_colors <- colorRamp2(c(-2, 0, 2), c("#0088ff", "#000000", "#fff047"))
p <- Heatmap(heatmap_data,
             cluster_rows = row_clust,
             cluster_columns = col_clust,
             top_annotation = top_annot,
             show_row_names = FALSE,
             show_column_names = FALSE,
             column_names_gp = gpar(fontsize = 6),
             name = "expr",
             row_title = "Intrinsic genes",
             #column_title = "n2486 DWD combined",
             col = heatmap_colors, 
             use_raster = F,
             raster_quality = 2, 
             raster_device = "png")

png(paste0(out_dir, "UNC337_n239_intrinsic_cluster.png"), 
    width = 15, height = 10, units = "in", res = 300)  # dimensions in inches
draw(p, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()


