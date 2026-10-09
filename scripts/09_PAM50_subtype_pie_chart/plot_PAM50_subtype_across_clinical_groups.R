# Plot pie charts of intrinsic and molecular subtypes

library(tidyverse)
library(ggplot2)

out_dir <- "output/PAM50_subtype_pie_chart/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

subtype_colors <- c("Basal" = "#c91508","Her2" = "#ed8ccb",
                    "LumA" = "#084dc4", "LumB" = "#9cd0f0",
                    "Normal" = "#4fb30c", "Claudin-low" = "goldenrod")

clin_colors <- c("TNBC" = "darkred", "HR-/HER2+" = "#CC79A7", 
                 "HR+/HER2-" = "darkblue", "HR+/HER2+" = "darkgreen")

clin_colors_3 <- c("TNBC" = "darkred", 
                   "HER2+" = "#CC79A7", 
                   "HR+/HER2-" = "darkblue")

metadata <- read.table("output/PAM50z_subtyping/n2486_combined_meta_w_PAM50_CLOW_subtypes.txt", 
                   sep = "\t", header = T, check.names = F)

## All samples and clinical subtypes (n=2486) ----------------------------------

pie_data <- metadata %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n2486_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes (n=2486)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## All samples and clinical subtypes (n=2486) ----------------------------------

pie_data <- metadata %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n2486_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes (n=2486)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## TNBC only (n = 436) ---------------------------------------------------------

tnbc_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "TNBC") %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n436_TNBC_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(tnbc_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in TNBC samples (n=436)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## HR+/HER2- only (n = 1699) ---------------------------------------------------

hr_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "HR+/HER2-") %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n1699_HR+HER2-_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(hr_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in HR+/HER2- samples (n=1699)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## HR+/HER2- only (n = 1699) ---------------------------------------------------

hr_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "HR+/HER2-") %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n1699_HR+HER2-_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(hr_pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes in HR+/HER2- samples (n=1699)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()



## HER2+ only (n = 351) --------------------------------------------------------

her2_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "HER2+") %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n351_HER2+_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(her2_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in HER2+ samples (n=351)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## LumA + LumB only (n = 1515) ---------------------------------------------------

lum_pie_data <- metadata %>% 
  filter(PAM50z_CLOW_cluster %in% c("LumA", "LumB")) %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n1515_LumA+B_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(lum_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in LumA+LumB samples (n=1515)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## LumA + LumB only (n = 1515) ---------------------------------------------------

lum_pie_data <- metadata %>% 
  filter(PAM50z_CLOW_cluster %in% c("LumA", "LumB")) %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "n1515_LumA+B_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(lum_pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in LumA+LumB samples (n=1515)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## n2486 (No/Bone) ------

#### PAM50 subtype 
pie_data <- metadata %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Bone")) %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "All_no-bone_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in all samples\n(no relapse or bone only)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

#### Clinical subtype 
pie_data <- metadata %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Bone")) %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "All_no-bone_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes in all samples\n(no relapse or bone only)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## HR+/HER2- (No/Bone) ------

#### PAM50 subtype 
hr_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "HR+/HER2-") %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Bone")) %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "HR+HER2-_no-bone_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(hr_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in HR+/HER2-\nsamples (no relapse or bone only)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

#### Clinical subtype 
hr_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "HR+/HER2-") %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Bone")) %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "HR+HER2-_no-bone_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(hr_pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes in HR+/HER2-\nsamples (no relapse or bone only)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

## LumA+LumB (No/Bone) ------

#### PAM50 subtype 
lum_pie_data <- metadata %>% 
  filter(PAM50z_CLOW_cluster %in% c("LumA", "LumB")) %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Bone")) %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "LumA+B_no-bone_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(lum_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in LumA+B\nsamples (no relapse or bone only)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

#### Clinical subtype 
lum_pie_data <- metadata %>% 
  filter(PAM50z_CLOW_cluster %in% c("LumA", "LumB")) %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Bone")) %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "LumA+B_no-bone_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(lum_pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes in LumA+B\nsamples (no relapse or bone only)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off

## n2486 (No/Multimet) ------

#### PAM50 subtype 
pie_data <- metadata %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Multimet")) %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "All_no-multimet_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in all samples\n(no relapse or multimet only)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

#### Clinical subtype 
pie_data <- metadata %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Multimet")) %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "All_no-multimet_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes in all samples\n(no relapse or multimet only)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()


## HR+/HER2- (No/Multimet) ------

#### PAM50 subtype 
hr_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "HR+/HER2-") %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Multimet")) %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "HR+HER2-_no-multimet_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(hr_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in HR+/HER2-\nsamples (no relapse or multimet only)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

#### Clinical subtype 
hr_pie_data <- metadata %>% 
  filter(Clinical_Group_3 == "HR+/HER2-") %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Multimet")) %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "HR+HER2-_no-multimet_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(hr_pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes in HR+/HER2-\nsamples (no relapse or multimet only)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

## LumA+LumB (No/Multimet) ------

#### PAM50 subtype 
lum_pie_data <- metadata %>% 
  filter(PAM50z_CLOW_cluster %in% c("LumA", "LumB")) %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Multimet")) %>% 
  count(PAM50z_CLOW_cluster) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(PAM50z_CLOW_cluster, "\n", "(", perc, "%)")) %>% 
  mutate(PAM50z_CLOW_cluster = factor(PAM50z_CLOW_cluster, levels = names(subtype_colors))) %>% 
  arrange(desc(PAM50z_CLOW_cluster)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "LumA+B_no-multimet_only_PAM50_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(lum_pie_data, aes(x = "", y = n, fill = PAM50z_CLOW_cluster)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of PAM50z subtypes in LumA+B\nsamples (no relapse or multimet only)", fill = "Subtype") + 
  scale_fill_manual(values = subtype_colors) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()

#### Clinical subtype 
lum_pie_data <- metadata %>% 
  filter(PAM50z_CLOW_cluster %in% c("LumA", "LumB")) %>% 
  filter(Relapse_Site %in% c("No_Relapse", "Multimet")) %>% 
  count(Clinical_Group_3) %>% 
  mutate(perc = round(n / sum(n) * 100, 1), 
         label = paste0(Clinical_Group_3, "\n", "(", perc, "%)")) %>% 
  mutate(Clinical_Group_3 = factor(Clinical_Group_3, levels = names(clin_colors_3))) %>% 
  arrange(desc(Clinical_Group_3)) %>% 
  mutate(cumulative = cumsum(n),
         midpoint = cumulative - n / 2, 
         label_pos = midpoint / sum(n) * 2 * pi )

pdf(paste0(out_dir, "LumA+B_no-multimet_only_clinical_subtype_pie_chart.pdf"), 
    width = 6, height = 4)
ggplot(lum_pie_data, aes(x = "", y = n, fill = Clinical_Group_3)) + 
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar("y") + 
  labs(title = "Distribution of clinical subtypes in LumA+B\nsamples (no relapse or multimet only)", fill = "Subtype") + 
  scale_fill_manual(values = clin_colors_3) + 
  theme_void() + 
  geom_text(aes(x = 1.4, # moves label outside of pie
                y = midpoint, 
                label = label),
            size = 4, 
            inherit.aes = FALSE)
dev.off()