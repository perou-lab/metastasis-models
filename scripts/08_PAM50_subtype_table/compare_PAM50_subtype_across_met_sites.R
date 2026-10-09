# Compare metastatic sites across PAM50 subtypes

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

out_dir <- "output/PAM50_subtype_table/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Load in metadata with PAM50 subtype annotations
metadata <- read.table("output/PAM50z_subtyping/n2486_combined_meta_w_PAM50_CLOW_subtypes.txt", 
                       sep = "\t", header = T, check.names = F) %>% 
  filter(Relapse_Site != "No_Relapse") %>% 
  filter(!(PAM50z_CLOW_cluster %in% c("Claudin-low", "Normal"))) %>% 
  mutate(Relapse_Site = factor(Relapse_Site, levels = c("Bone", "Brain", "Liver", "Lung", "Multimet")))

## Generate contingency table
cont_table <- table(metadata$PAM50z_CLOW_cluster, metadata$Relapse_Site)
cont_table

## Test for independence (chi-square)
csq_res <- chisq.test(cont_table)

## Report chi-square p-value
sprintf("%.3g", csq_res$p.value)

## Standardized residuals
stdres <- csq_res$stdres
round(stdres,2)

## Visualize standard residuals
hm <- Heatmap(stdres,
              name = "Std\nResidual",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 18),
              column_names_gp = gpar(fontsize = 18),
              column_names_rot = 0,
              col = colorRamp2(c(-8, 0, 8), c("#2166AC", "white", "#B2182B")),
              rect_gp = gpar(col = "black"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(round(stdres[i, j], 1), x = x, y = y, gp = gpar(fontsize = 20))},
              heatmap_legend_param = list(
                title = "Residual",
                at = c(-8, -4, 0, 4, 8)))

pdf(paste0(out_dir, "PAM50_subtype_chisq_heatmap.pdf"), 
    width = 7, height = 4.5)
draw(hm)
dev.off()

png(paste0(out_dir, "PAM50_subtype_chisq_heatmap.png"), 
    width = 7, height = 4.5, units = "in", res = 300)
draw(hm)
dev.off()
