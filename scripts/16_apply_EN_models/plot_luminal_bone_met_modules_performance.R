# Forest plot
  # Datasets = UNC-337, SCAN-B
  # Model: Luminal Bone Metastasis
  # Input: Modules

library(survival)
library(survminer)
library(ggplot2)
library(tidyverse)
library(patchwork)

out_dir <- "output/apply_EN_models/"

## SCAN-B data / scores:
scanb_meta <- read.table("../datasets/SCANB_gitlab/SCAN-B_6329_Clinical_Data.txt", 
                       sep = "\t", header = T, check.names = F)

scanb_bone_modules <- read.table("output/apply_EN_models/SCANB/modules_935+3/lum_bone/SCANB_luminal_bone_module_EN_model_scores.tsv", 
                             sep = "\t", header = T, check.names = F) %>% 
  rename(Module_EN_Score = Luminal_Bone_EN_score)

scanb_score <- scanb_meta %>% 
  left_join(scanb_bone_modules, by = "Merged") %>% 
  filter(PAM50z_Call_CLOW %in% c("LumA", "LumB"))

## UNC-337 data / scores: 
unc337_meta <- read.table("output/apply_EN_models/process_test_sets/UNC337/final/UNC337_n239_metadata.tsv", 
                          sep = "\t", header = T, check.names = F)

unc337_bone_modules <- read.table("output/apply_EN_models/UNC337/modules_935+3/lum_bone/UNC337_luminal_bone_module_EN_model_scores.tsv",
                                    sep = "\t", header = T, check.names = F) %>% 
  rename(Module_EN_Score = Luminal_Bone_EN_score)

unc337_score <- unc337_meta %>% 
  left_join(unc337_bone_modules, by = "Sample_GEO_accession") %>% 
  mutate(RFS_years = RFS_months/12) %>% 
  mutate(OS_years = OS_months/12) %>% 
  filter(PAM50_call %in% c("LumA", "LumB"))

## Standardize model scores within datasets: 
scanb_score$score_z <- as.numeric(scale(scanb_score$Module_EN_Score))
unc337_score$score_z <- as.numeric(scale(unc337_score$Module_EN_Score))

## Calculate Coxph values for each score
scanb_os_cox <- coxph(Surv(OS_years, OS_event) ~ score_z,
             data = scanb_score)

unc_rfs_cox <- coxph(Surv(RFS_years, RFS_event) ~ score_z, 
                 data = unc337_score)

unc_os_cox <- coxph(Surv(OS_years, OS_event) ~ score_z, 
                    data = unc337_score)

## Assemble data for plotting
scanb_hr <- exp(coef(scanb_os_cox))
scanb_ci <- exp(confint(scanb_os_cox))[1,]
scanb_p <- summary(scanb_os_cox)$coefficients[,"Pr(>|z|)"]

unc_rfs_hr <- exp(coef(unc_rfs_cox))
unc_rfs_ci <- exp(confint(unc_rfs_cox))[1,]
unc_rfs_p <- summary(unc_rfs_cox)$coefficients[,"Pr(>|z|)"]

unc_os_hr <- exp(coef(unc_os_cox))
unc_os_ci <- exp(confint(unc_os_cox))[1,]
unc_os_p <- summary(unc_os_cox)$coefficients[,"Pr(>|z|)"]

df_forest <- data.frame(Study = c("SCAN-B OS", "UNC-337 OS", "UNC-337 RFS"),
                        HR = c(scanb_hr, unc_os_hr, unc_rfs_hr),
                        lower = c(scanb_ci[1], unc_os_ci[1], unc_rfs_ci[1]),
                        upper = c(scanb_ci[2], unc_os_ci[2], unc_rfs_ci[2]), 
                        p_value = c(scanb_p, unc_os_p, unc_rfs_p))

df_forest <- df_forest %>% 
  mutate(label = paste0(round(HR,2), " (", round(lower,2), "-", round(upper,2), "), p = ", signif(p_value,3))) %>% 
  mutate(Study = factor(Study, levels = c("UNC-337 RFS", "UNC-337 OS", "SCAN-B OS")))

## Plot forest plot + table
p_forest <- ggplot(df_forest, aes(x = HR, y = Study)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.15) +
  geom_point(size = 2.8) +
  geom_vline(xintercept = 1,
             linetype = "dashed") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 5, 10)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 15)
        ) +
  xlab("Hazard Ratio (HR per SD)") +
  ylab("")

forest_table <- ggplot(df_forest, aes(y = Study)) + 
  geom_text(aes(x = 1, label = label),
            hjust = 0, 
            size = 5) + 
  xlim(0.8, 2) +
  theme_void() +
  ggtitle("")

plot <- p_forest + forest_table + 
  plot_layout(width = c(0.5,1.5))

ggsave(paste0(out_dir, "/Luminal_bone_met_modules_forest_plot_FORMAT.png"), 
       plot, 
       width = 8, height = 2, units = "in", dpi = 300)

pdf(paste0(out_dir, "/Luminal_bone_met_modules_forest_plot_FORMAT.pdf"),
    width = 8, height = 2)
print(plot)
dev.off()
