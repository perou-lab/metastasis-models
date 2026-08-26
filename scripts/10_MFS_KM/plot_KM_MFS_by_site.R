## Site specific MFS by subtype KM plots
# Plot by site of metastasis

library(tidyverse)
library(matrixStats)
library(survival)
library(survminer)

out_dir <- "output/MFS_KM/Site_of_Metastasis/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Color annotations -----------------------------------------------------------
cohort_colors <- c("GSE12276" = "#bbb0f8",
                  "GSE2034" = "#F8766D", 
                  "GSE2603" = "#aad6aa",
                  "METABRIC" = "#fdd246", 
                  "NKI295" = "#426080")

clin_colors <- c("TNBC" = "darkred",
                 "HR-/HER2+" = "#CC79A7", 
                 "HR+/HER2-" = "darkblue", 
                 "HR+/HER2+" = "darkgreen")

clin_colors_3 <- c("TNBC" = "darkred", 
                   "HER2+" = "#CC79A7", 
                   "HR+/HER2-" = "darkblue")

relapse_site_colors <- c("Bone" = "#a18be8", "Brain" = "#a11045",
                         "Liver" = "#a9c451", "Lung" = "#cc8014",
                         "Multimet" = "#1492cc", "No_Relapse" = "lightgray")

pam50_colors <- c("Basal" = "#c91508", 
									"Her2" = "#ed8ccb",
									"LumA" = "#084dc4",
									"LumB" = "#9cd0f0",
									"Normal" = "#4fb30c", 
									"Claudin-low" = "goldenrod")

## Load in metadata ------------------------------------------------------------
metadata <- read.table("output/PAM50z_subtyping/n2486_combined_meta_w_PAM50_CLOW_subtypes.txt", 
                       sep = "\t", header = T, check.names = F) %>% 
  rename(PAM50_CLOW = PAM50z_CLOW_cluster) %>% 
  mutate(MFS_Years = MFS_Months/12)


## Plotting Metastasis Free Survival by Site -----------------------------------
site <- "Any Metastasis"
file_name <- "By_Site"
title <- paste0("Any Distant Metastasis")

plot_meta <- metadata %>% 
  filter(Relapse_Site != "No_Relapse") %>% 
  mutate(Relapse_Site = factor(Relapse_Site, levels = c("Bone", "Brain", "Liver", "Lung", "Multimet"))) %>% 
  mutate(MFS_Years = MFS_Months/12)

surv_object <- Surv(time = plot_meta$MFS_Years, event = plot_meta$Relapse_Event)
site <- plot_meta$Relapse_Site
fit1 <- survfit(surv_object ~ site, data = plot_meta)

site_levels <- levels(factor(plot_meta$Relapse_Site))
site_palette <- relapse_site_colors[site_levels]

names(site_palette) <- paste0("site=", names(site_palette))
#color_vector <- unname(relapse_site_colors[strata_names])

## KM Plot
p <- ggsurvplot(
  fit1,
  data = plot_meta,
 # pval = TRUE,
  title = title,
  palette = site_palette, 
  font.tickslab = c(15),
  break.x.by = 5,
  xlab = "Years",
  font.x = c(20),
  xlim = c(0,15),
  break.y.by = 0.2,
  ylab = "Percent MFS",
  font.y = c(15),
  linetype = 1,
  ggtheme = theme_classic(),
  axes.offset = FALSE, 
)

# Log Rank Test
surv_diff <- survdiff(surv_object ~ site, data = plot_meta)
pval <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
pval_text <- paste0("Log-rank p = ", signif(pval, 4))
p$plot <- p$plot +
  annotate("text", x = 2.2, y = 0.05, label = pval_text, size = 4) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 40, 10, 10)) + 
  coord_cartesian(xlim = c(0, 15))

# Save plot as .pdf and .png files
ggsave(filename = paste0(out_dir, "MFS_Site_of_Relapse.pdf"),
       plot = p$plot, 
       device = "pdf", 
       width = 8, 
       height = 6, 
       units = "in")

ggsave(filename = paste0(out_dir, "MFS_Site_of_Relapse.png"),
       plot = p$plot,
       device = "png",
       dpi = 300,
       width = 8,
       height = 6,
       units = "in")


