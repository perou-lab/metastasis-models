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

out_dir <- "output/intrinsic_cluster/"

cohort_colors <- c("GSE12276" = "#082a54",
                   "GSE2034" = "#59a89c", 
                   "GSE2603" = "#a559aa",
                   "METABRIC" = "#cecece", 
                   "NKI295" = "#e02b35")

clin_colors_3 <- c("TNBC" = "darkred", 
                   "HER2+" = "#CC79A7", 
                   "HR+/HER2-" = "darkblue")

relapse_site_colors <- c("Bone" = "#a18be8", "Brain" = "#a11045",
                         "Liver" = "#a9c451", "Lung" = "#cc8014",
                         "Multi-Site" = "#1492cc", "No_Relapse" = "lightgray")

pam50_colors <- c("Basal" = "#c91508", 
                  "Her2" = "#ed8ccb",
                  "LumA" = "#084dc4",
                  "LumB" = "#9cd0f0",
                  "Normal" = "#4fb30c", 
                  "Claudin-low" = "goldenrod")

train_test_colors <- c("Train" = "#298c8c", "Test" = "#800074")

## Metadata 
metadata <- read.table("output/balanced_stratification/n2486/n2486_complete_meta_w_onlysite_stratification.txt",
                       sep = "\t", header = T, check.names = F) %>% 
  mutate(Relapse_Site = case_when(
    Relapse_Site == "Multimet" ~ "Multi-Site", 
    TRUE ~ as.character(Relapse_Site)
  ))

## Intrinsic Clustering
heatmap_data <- read.table("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols_intrinsic_genes.txt",
                           sep = "\t", header = T, row.names = 1, check.names = F)
heatmap_data <- as.matrix(heatmap_data)

row_clust <- xcluster2r("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols_intrinsic_genes.gtr", 
                        distance = "pearson")
col_clust <- xcluster2r("output/intrinsic_cluster/n2486_Ensembl_dwd_symbols_intrinsic_genes.atr", 
                          distance = "pearson")

heatmap_metadata <- metadata[match(colnames(heatmap_data), metadata$Sample_ID),]

## Heatmap Annotations
top_annot <- HeatmapAnnotation(Study = heatmap_metadata$Cohort, 
                               Clinical_Group = heatmap_metadata$Clinical_Group_3,
                               PAM50z_Subtype = heatmap_metadata$PAM50z_CLOW_cluster,
                               Site_of_Metastasis = heatmap_metadata$Relapse_Site,
                               n2486_Train_Test = heatmap_metadata$Train75_Test25,
                               col = list(Site_of_Metastasis = relapse_site_colors, 
                                          Clinical_Group = clin_colors_3,
                                          PAM50z_Subtype = pam50_colors,
                                          Study = cohort_colors, 
                                          n2486_Train_Test = train_test_colors),
                               show_legend = c(T, T, T, T),
                               annotation_legend_param = list(
                                 Clinical_Group = list(title = "Clinical Subtype"), 
                                 PAM50z_Subtype = list(title = "PAM50 Subtype"), 
                                 Site_of_Metastasis = list(title = "Site of Metastasis"), 
                                 n2486_Train_Test = list(title = "Train/Test")))

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

png(paste0(out_dir, "n2486_dwd_intrinsic_cluster_n2486.png"), 
    width = 15, height = 10, units = "in", res = 300)  # dimensions in inches
draw(p, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

pdf(paste0(out_dir, "n2486_dwd_intrinsic_cluster_n2486.pdf"),
  width = 15, height = 10)

draw(p,
     padding = unit(c(10, 10, 10, 10), "mm"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

dev.off()
