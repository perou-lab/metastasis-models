# Evaluate the performance of breast cancer prognostic signatures

library(tidyverse)
library(ROCR)

out_dir <- "output/compare_prognostic_assays/"
dir.create(paste0(out_dir), showWarnings = F, recursive = T)

## Functions -------------------------------------------------------------------
plot_ROC2 <- function(perf1, a1, set1) {
  par(mai = c(0.24, 0.24, 0.05, 0.05), cex.axis = 0.3)
  plot(perf1, col = "red", lwd = 1.2, xlab = "", ylab = "", box.lwd=0.8,
       xaxis.xaxt = "n", yaxis.yaxt = "n")
  abline(a=0, b=1, lwd = 0.8)
  axis(1, tck=(-0.02), lwd = 0.8)
  axis(2, tck=(-0.02), lwd = 0.8)
  mtext(side = 1, text = seq(0, 1, by=0.2), at = seq(0, 1, by=0.2), cex = 0.35, line = (-0.3))
  mtext(side = 2, text = seq(0, 1, by=0.2), at = seq(0, 1, by=0.2), cex = 0.35, line = 0.1)
  mtext(side = 2, text = "True positive rate", at = 0.5, cex = 0.4, line = 0.6)
  mtext(side = 1, text = "False positive rate", at = 0.5, cex = 0.4, line = 0.2)
  legend("bottomright", legend = c(paste0(set1, " (AUC = ", round(a1, 2), ")")),
         lty = c(1,1), lwd = c(1,1),
         col = c('red'), cex = 0.35, bty = "n")
}

## Read in metadata and modules ------------------------------------------------
metadata <- read.table("output/balanced_stratification/n2486/n2486_complete_meta_w_onlysite_stratification.txt",
                   sep = "\t", header = T) %>% 
  mutate(No_v_Distant = factor(No_v_Distant, levels = c("No_Relapse", "Distant_Relapse")))

modules <- read.table("output/calculate_modules/n2486_Ensembl_dwd_938_published_modules.txt",
                      sep = "\t", header = T)

metadata <- metadata[match(colnames(modules), metadata$Sample_ID), ]

x <- as.matrix(t(as.matrix(modules)))
y <- metadata$No_v_Distant

## Evaluate ROR_S_Model performance (Prosigna) ---------------------------------
pred_prosigna <- x[, "ROR_P_Model_JCO.2009_PMID.19204204"]
pred_prosigna <- prediction(pred_prosigna,
                            labels = y, 
                            label.ordering = c("No_Relapse", "Distant_Relapse"))
perf_prosigna <- performance(pred_prosigna,
                             measure = "tpr", 
                             x.measure = "fpr")
auc_prosigna <- signif(performance(pred_prosigna, measure = 'auc')@y.values[[1]][1], 4)

png(paste0(out_dir, "AUC_n2486_ROR_P_Model.png"), width = 1.5, height = 1.5, units = "in", res = 350)
plot_ROC2(perf_prosigna, auc_prosigna, "ROR_P_Model")
dev.off()

pdf(paste0(out_dir, "AUC_n2486_ROR_P_Model.pdf"), width = 1.5, height = 1.5)
plot_ROC2(perf_prosigna, auc_prosigna, "ROR_P_Model")
dev.off()

## GHI_RS_Model_NEJM.2004_PMID.15591335 (OncotypeDX) ---------------------------
pred_oncDX <- x[, "GHI_RS_Model_NEJM.2004_PMID.15591335"]
pred_oncDX <- prediction(pred_oncDX,
                         labels = y, 
                         label.ordering = c("No_Relapse", "Distant_Relapse"))
perf_oncDX <- performance(pred_oncDX,
                          measure = "tpr", 
                          x.measure = "fpr")
auc_oncDX <- signif(performance(pred_oncDX, measure = 'auc')@y.values[[1]][1], 4)

png(paste0(out_dir, "AUC_n2486_GHI_RS_Model.png"), width = 1.5, height = 1.5, units = "in", res = 350)
plot_ROC2(perf_oncDX, auc_oncDX, "GHI_RS_Model")
dev.off()

pdf(paste0(out_dir, "AUC_n2486_GHI_RS_Model.pdf"), width = 1.5, height = 1.5)
plot_ROC2(perf_oncDX, auc_oncDX, "GHI_RS_Model")
dev.off()

## Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860 (Mammaprint)
pred_mammaprint <- x[, "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860"]
pred_mammaprint <- prediction(pred_mammaprint,
                              labels = y, 
                              label.ordering = c("Distant_Relapse", "No_Relapse"))
perf_mammaprint <- performance(pred_mammaprint,
                               measure = "tpr",
                               x.measure = "fpr")
auc_mammaprint <- signif(performance(pred_mammaprint, measure = 'auc')@y.values[[1]][1], 4)

png(paste0(out_dir, "AUC_n2486_PCorr_NKI70.png"), width = 1.5, height = 1.5, units = "in", res = 350)
plot_ROC2(perf_mammaprint, auc_mammaprint, "PCorr_NKI70")
dev.off()

pdf(paste0(out_dir, "AUC_n2486_PCorr_NKI70.pdf"), width = 1.5, height = 1.5)
plot_ROC2(perf_mammaprint, auc_mammaprint, "PCorr_NKI70")
dev.off()