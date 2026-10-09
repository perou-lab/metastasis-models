# Compare metastatic sites across clinical subtypes

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

out_dir <- "output/clinical_subtype_table/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Load in metadata with clinical subtype annotations
metadata <- read.table("output/processing_metadata/n2486_combined_clinical_data.tsv", 
                       sep = "\t", header = T, check.names = F) %>% 
  select(Sample_ID, Relapse_Site)

receptor_data <- read.table("output/infer_ER_PR_HER2/n2486_Ensembl_dwd_inferred_clin_groups.txt", 
                            sep = "\t", header = T) %>% 
  mutate(Clinical_Group_3 = case_when(
    (ER_status_w_inferred == "Positive" | PR_status_w_inferred == "Positive") & HER2_status_w_inferred == "Negative" ~ "HR+/HER2-", 
    ER_status_w_inferred == "Negative" & PR_status_w_inferred == "Negative" & HER2_status_w_inferred == "Negative" ~ "TNBC", 
    HER2_status_w_inferred == "Positive" ~ "HER2+"))

metadata <- left_join(metadata, receptor_data, by = "Sample_ID") %>% 
  filter(Relapse_Site != "No_Relapse") %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = c("HR+/HER2-", "HER2+", "TNBC")))

## Generate contingency table
cont_table <- table(metadata$Clinical_Group_3, metadata$Relapse_Site)
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
                grid.text(round(stdres[i, j], 1), x = x, y = y, gp = gpar(fontsize = 14))},
              heatmap_legend_param = list(title = "Residual", at = c(-8, -4, 0, 4, 8)))

pdf(paste0(out_dir, "Clinical_subtype_chisq_heatmap.pdf"), 
    width = 7, height = 3.2)
draw(hm)
dev.off()

png(paste0(out_dir, "Clinical_subtype_chisq_heatmap.png"), 
    width = 7, height = 3.2, units = "in", res = 300)
draw(hm)
dev.off()
