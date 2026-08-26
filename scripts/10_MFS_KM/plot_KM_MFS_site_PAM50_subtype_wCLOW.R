## Site specific MFS by subtype KM plots
# Plot by PAM50 subtypes (with claudin-low samples)

library(tidyverse)
library(matrixStats)
library(survival)
library(survminer)

out_dir <- "output/MFS_KM/PAM50_w_CLOW/"
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


## Plotting Distant Metastasis Free Survival -----------------------------------
site <- "Any Metastasis"
file_name <- "Any_Site"
title <- paste0("Any Distant Metastasis")

plot_meta <- metadata %>% 
  mutate(PAM50_CLOW = factor(PAM50_CLOW, levels = c("Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal")))

surv_object <- Surv(time = plot_meta$MFS_Years, event = plot_meta$Relapse_Event)
subtype <- plot_meta$PAM50_CLOW
fit1 <- survfit(surv_object ~ subtype, data = plot_meta)

subtype_levels <- levels(factor(plot_meta$PAM50_CLOW))
subtype_palette <- pam50_colors[subtype_levels]

strata_name <- gsub("subtype=","",names(fit1$strata))
color_vector <- unname(pam50_colors[strata_name])

# KM Plot
p <- ggsurvplot(
  fit1,
  data = plot_meta,
  pval = F,
  title = title,
  palette = color_vector, 
  font.tickslab = c(15),
  break.x.by = 5,
  xlab = "Years",
  font.x = c(20),
  break.y.by = 0.1,
  ylab = "% Distant Relapse",
  font.y = c(20),
  linetype = 1,
  ylim = c(0.4, 1),
  ggtheme = theme_classic(),
  axes.offset = FALSE, 
  font.main = c(25), 
  legend = "none"
)

# Log Rank Test
surv_diff <- survdiff(surv_object ~ subtype, data = plot_meta)
pval <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
pval_text <- paste0("Log-rank p = ", signif(pval, 4))
p$plot <- p$plot +
  annotate("text", x = 7, y = 0.45, label = pval_text, size = 6) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 40, 10, 10))

# Save plot as .pdf and .png files
ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".pdf"),
       plot = p$plot, 
       device = "pdf", 
       width = 8, 
       height = 6, 
       units = "in")

ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".png"),
       plot = p$plot,
       device = "png",
       dpi = 300,
       width = 8,
       height = 6,
       units = "in")


## Plotting Bone Only Metastasis Free Survival -----------------------------------
site_label <- "Bone"
title <- paste0("Bone Only Metastasis")
y_lab_label <- paste0("% ", site_label, " Only Relapse")
file_name <- "Bone_Only"

plot_meta <- metadata %>% 
  mutate(PAM50_CLOW = factor(PAM50_CLOW, levels = c("Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>% 
  mutate(Site = ifelse(Relapse_Site == site_label, 1, 0))

surv_object <- Surv(time = plot_meta$MFS_Years, event = plot_meta$Site)
subtype <- plot_meta$PAM50_CLOW
fit1 <- survfit(surv_object ~ subtype, data = plot_meta)

subtype_levels <- levels(factor(plot_meta$PAM50_CLOW))
subtype_palette <- pam50_colors[subtype_levels]

strata_name <- gsub("subtype=","",names(fit1$strata))
color_vector <- unname(pam50_colors[strata_name])

# KM Plot
p <- ggsurvplot(
  fit1,
  data = plot_meta,
  pval = F,
  title = title,
  palette = color_vector, 
  font.tickslab = c(15),
  break.x.by = 5,
  xlab = "Years",
  font.x = c(20),
  break.y.by = 0.1,
  ylab = "% Bone Only Relapse",
  font.y = c(20),
  linetype = 1,
  ylim = c(0.4, 1),
  ggtheme = theme_classic(),
  axes.offset = FALSE, 
  font.main = c(25), 
  legend = "none"
)

# Log Rank Test
surv_diff <- survdiff(surv_object ~ subtype, data = plot_meta)
pval <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
pval_text <- paste0("Log-rank p = ", signif(pval, 4))
p$plot <- p$plot +
  annotate("text", x = 7, y = 0.45, label = pval_text, size = 6) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 40, 10, 10))

# Save plot as .pdf and .png files
ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".pdf"),
       plot = p$plot, 
       device = "pdf", 
       width = 8, 
       height = 6, 
       units = "in")

ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".png"),
       plot = p$plot,
       device = "png",
       dpi = 300,
       width = 8,
       height = 6,
       units = "in")


## Plotting Brain Only Metastasis Free Survival -----------------------------------
site_label <- "Brain"
title <- paste0("Brain Only Metastasis")
y_lab_label <- paste0("% ", site_label, " Only Relapse")
file_name <- "Brain_Only"

plot_meta <- metadata %>% 
  mutate(PAM50_CLOW = factor(PAM50_CLOW, levels = c("Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>% 
  mutate(Site = ifelse(Relapse_Site == site_label, 1, 0))

surv_object <- Surv(time = plot_meta$MFS_Years, event = plot_meta$Site)
subtype <- plot_meta$PAM50_CLOW
fit1 <- survfit(surv_object ~ subtype, data = plot_meta)

subtype_levels <- levels(factor(plot_meta$PAM50_CLOW))
subtype_palette <- pam50_colors[subtype_levels]

strata_name <- gsub("subtype=","",names(fit1$strata))
color_vector <- unname(pam50_colors[strata_name])

# KM Plot
p <- ggsurvplot(
  fit1,
  data = plot_meta,
  pval = F,
  title = title,
  palette = color_vector, 
  font.tickslab = c(15),
  break.x.by = 5,
  xlab = "Years",
  font.x = c(20),
  break.y.by = 0.1,
  ylab = "% Brain Only Relapse",
  font.y = c(20),
  linetype = 1,
  ylim = c(0.4, 1),
  ggtheme = theme_classic(),
  axes.offset = FALSE, 
  font.main = c(25), 
  legend = "none"
)

# Log Rank Test
surv_diff <- survdiff(surv_object ~ subtype, data = plot_meta)
pval <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
pval_text <- paste0("Log-rank p = ", signif(pval, 4))
p$plot <- p$plot +
  annotate("text", x = 7, y = 0.45, label = pval_text, size = 6) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 40, 10, 10))

# Save plot as .pdf and .png files
ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".pdf"),
       plot = p$plot, 
       device = "pdf", 
       width = 8, 
       height = 6, 
       units = "in")

ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".png"),
       plot = p$plot,
       device = "png",
       dpi = 300,
       width = 8,
       height = 6,
       units = "in")


## Plotting Liver Only Metastasis Free Survival -----------------------------------
site_label <- "Liver"
title <- paste0("Liver Only Metastasis")
y_lab_label <- paste0("% ", site_label, " Only Relapse")
file_name <- "Liver_Only"

plot_meta <- metadata %>% 
  mutate(PAM50_CLOW = factor(PAM50_CLOW, levels = c("Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>% 
  mutate(Site = ifelse(Relapse_Site == site_label, 1, 0))

surv_object <- Surv(time = plot_meta$MFS_Years, event = plot_meta$Site)
subtype <- plot_meta$PAM50_CLOW
fit1 <- survfit(surv_object ~ subtype, data = plot_meta)

subtype_levels <- levels(factor(plot_meta$PAM50_CLOW))
subtype_palette <- pam50_colors[subtype_levels]

strata_name <- gsub("subtype=","",names(fit1$strata))
color_vector <- unname(pam50_colors[strata_name])

# KM Plot
p <- ggsurvplot(
  fit1,
  data = plot_meta,
  pval = F,
  title = title,
  palette = color_vector, 
  font.tickslab = c(15),
  break.x.by = 5,
  xlab = "Years",
  font.x = c(20),
  break.y.by = 0.1,
  ylab = "% Liver Only Relapse",
  font.y = c(20),
  linetype = 1,
  ylim = c(0.4, 1),
  ggtheme = theme_classic(),
  axes.offset = FALSE, 
  font.main = c(25), 
  legend = "none"
)

# Log Rank Test
surv_diff <- survdiff(surv_object ~ subtype, data = plot_meta)
pval <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
pval_text <- paste0("Log-rank p = ", signif(pval, 4))
p$plot <- p$plot +
  annotate("text", x = 7, y = 0.45, label = pval_text, size = 6) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 40, 10, 10))

# Save plot as .pdf and .png files
ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".pdf"),
       plot = p$plot, 
       device = "pdf", 
       width = 8, 
       height = 6, 
       units = "in")

ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".png"),
       plot = p$plot,
       device = "png",
       dpi = 300,
       width = 8,
       height = 6,
       units = "in")


## Plotting Lung Only Metastasis Free Survival -----------------------------------
site_label <- "Lung"
title <- paste0("Lung Only Metastasis")
y_lab_label <- paste0("% ", site_label, " Only Relapse")
file_name <- "Lung_Only"

plot_meta <- metadata %>% 
  mutate(PAM50_CLOW = factor(PAM50_CLOW, levels = c("Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>% 
  mutate(Site = ifelse(Relapse_Site == site_label, 1, 0))

surv_object <- Surv(time = plot_meta$MFS_Years, event = plot_meta$Site)
subtype <- plot_meta$PAM50_CLOW
fit1 <- survfit(surv_object ~ subtype, data = plot_meta)

subtype_levels <- levels(factor(plot_meta$PAM50_CLOW))
subtype_palette <- pam50_colors[subtype_levels]

strata_name <- gsub("subtype=","",names(fit1$strata))
color_vector <- unname(pam50_colors[strata_name])

# KM Plot
p <- ggsurvplot(
  fit1,
  data = plot_meta,
  pval = F,
  title = title,
  palette = color_vector, 
  font.tickslab = c(15),
  break.x.by = 5,
  xlab = "Years",
  font.x = c(20),
  break.y.by = 0.1,
  ylab = "% Lung Only Relapse",
  font.y = c(20),
  linetype = 1,
  ylim = c(0.4, 1),
  ggtheme = theme_classic(),
  axes.offset = FALSE, 
  font.main = c(25), 
  legend = "none"
)

# Log Rank Test
surv_diff <- survdiff(surv_object ~ subtype, data = plot_meta)
pval <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
pval_text <- paste0("Log-rank p = ", signif(pval, 4))
p$plot <- p$plot +
  annotate("text", x = 7, y = 0.45, label = pval_text, size = 6) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 40, 10, 10))

# Save plot as .pdf and .png files
ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".pdf"),
       plot = p$plot, 
       device = "pdf", 
       width = 8, 
       height = 6, 
       units = "in")

ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".png"),
       plot = p$plot,
       device = "png",
       dpi = 300,
       width = 8,
       height = 6,
       units = "in")


## Plotting MultiMet Free Survival ---------------------------------------------
site_label <- "Multimet"
title <- paste0("Multi-Site Metastasis")
y_lab_label <- paste0("% ", site_label, " Only Relapse")
file_name <- "Multimet"

plot_meta <- metadata %>% 
  mutate(PAM50_CLOW = factor(PAM50_CLOW, levels = c("Basal", "Claudin-low", "Her2", "LumA", "LumB", "Normal"))) %>% 
  mutate(Site = ifelse(Relapse_Site == site_label, 1, 0))

surv_object <- Surv(time = plot_meta$MFS_Years, event = plot_meta$Site)
subtype <- plot_meta$PAM50_CLOW
fit1 <- survfit(surv_object ~ subtype, data = plot_meta)

subtype_levels <- levels(factor(plot_meta$PAM50_CLOW))
subtype_palette <- pam50_colors[subtype_levels]

strata_name <- gsub("subtype=","",names(fit1$strata))
color_vector <- unname(pam50_colors[strata_name])

# KM Plot
p <- ggsurvplot(
  fit1,
  data = plot_meta,
  pval = F,
  title = title,
  palette = color_vector, 
  font.tickslab = c(15),
  break.x.by = 5,
  xlab = "Years",
  font.x = c(20),
  break.y.by = 0.1,
  ylab = "% Multi-Site Relapse",
  font.y = c(20),
  linetype = 1,
  ylim = c(0.4, 1),
  ggtheme = theme_classic(),
  axes.offset = FALSE, 
  font.main = c(25), 
  legend = "none"
)

# Log Rank Test
surv_diff <- survdiff(surv_object ~ subtype, data = plot_meta)
pval <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
pval_text <- paste0("Log-rank p = ", signif(pval, 4))
p$plot <- p$plot +
  annotate("text", x = 7, y = 0.45, label = pval_text, size = 6) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 40, 10, 10))

# Save plot as .pdf and .png files
ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".pdf"),
       plot = p$plot, 
       device = "pdf", 
       width = 8, 
       height = 6, 
       units = "in")

ggsave(filename = paste0(out_dir, "MFS_PAM50_w_CLOW_", file_name, ".png"),
       plot = p$plot,
       device = "png",
       dpi = 300,
       width = 8,
       height = 6,
       units = "in")


