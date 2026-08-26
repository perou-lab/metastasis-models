#' train_EN_site_AM.R
#'  
#' Bootstrapping workflow to train elastic net model on RNA features, and IHC
#' receptor status (ER/PR/HER2) using combined n5 (HARRELL + METABRIC) dataset
#' where training on % of data, and testing on % of data
#' 
#' Predicting site specific metastatic behavior in breast cancer

.libPaths("/home/amich/conda/envs/metastatic_site_EN/lib/R/library")

library(caret)
library(glmnet)
library(ROCR)
library(tidyverse)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(optparse)
library(ggplot2)
library(cutpointr)

## Input parameters ------------------------------------------------------------
option_list <- list(
  make_option("--seed", type = "numeric", default=124022, 
              help = "Seed to be used for EN model; default = 124022",
              metavar = "character"), 
  make_option("--percent_train", type = "character", default="75", 
              help = "Percentage of the data that will be used for training; 
              default = 75",
              metavar = "character"), 
  make_option("--percent_test", type = "character", default="25", 
              help = "Percentage of the data that will be used for testing;
              default = 25",
              metavar = "character"), 
  make_option("--relapse_site_1", type = "character", default=NULL, 
              help = "First relapse group for comparison (usually No Relapse).
              Options are: NoRel, NoOtherRel, OtherRel",
              metavar = "character"), 
  make_option("--relapse_site_2", type = "character", default=NULL, 
              help = "Second relapse group for comparison (usually distant or a 
              specific site). Options are: Distant, OnlyLung, OnlyBone, OnlyLiver,
              OnlyBrain, MultiMet, AnyLung, AnyBone, AnyLiver, AnyBrain,
              Distant_notBone", 
              metavar = "character"),
  make_option("--input_file", type = "character", default=NULL, 
              help = "Version of the RNA data or features to be included.
              Options are limma, dwd, pub_modules, all_modules",
              metavar = "character"),
  make_option("--resampling_method", type = "character", default="boot", 
              help = "Resampling method to be used for training the model. 
              Options are CV or boot; default is boot",
              metavar = "character"), 
  make_option("--metric", type = "character", default= "Accuracy", 
              help = "Metric used to evaluate performance of the model. 
              Options are Accuracy or Kappa",
              metavar = "character"), 
  make_option("--n_bootstraps", type = "numeric", default = 1000, 
              help = "Number of boostraps for training the model;
              default = 1000",
              metavar = "character"),
  make_option("--alphas", type = "character", default="0.1,0.9,0.1",
              help = "Alpha values to use, will generate based on 
							seq(from = from, to = to, by = by). Specify 'from,to,by' with 
							comma separation. Default is 0.1,0.9,0.1",
							metavar="character"),
  make_option("--n_lambdas", type = "numeric", default = 100,
              help = "Number of lambda values to consider; will generate based
              on the values of alpha; default is 100.",
              metavar="character"), 
  make_option("--n_features", type = "numeric", default = 30, 
              help = "Number of features to include in the bar plot.
              Final plot will have 2*n_features or all features if < 2*n_features
              are chosen in the final model; default is 30 (60 total)", 
              metavar="character"),
  make_option("--output_directory", type = "character", default = NULL, 
              help = "Path to output directory:
              e.g. ///Anna/metastatic_site_EN/output/Modules/n2486/", 
              metavar="character"),
  make_option("--metadata_directory", type = "character", default = NULL,
              help = "Path to directory holding metadata file: e.g.
              //datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_2/output/HER2+_only_analysis/balanced_stratification/"),
  make_option("--metadata_file", type = "character", default = NULL, 
              help = "Metadata file name: e.g.
              n364_HER2+_only_combined_metadata_balanced_stratification.txt")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# seed <- 124022
# percent_train <- "75"
# percent_test <- "25"
# relapse_site1 <- "NoRel"
# relapse_site2 <- "Distant"
# file <- "pub_modules"
# file <- "dwd"
# method <- "boot"
# train_metric <- "Accuracy"
# n_bootstraps <- 100
# alphas <- seq("0.1,0.9,0.1")
# n_lambdas <- 100
# output_directory <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_manuscript/output/train_EN_models/n2486/modules_933+3/"
# output_directory <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_manuscript/output/train_EN_models/n2486/genes_10814/OLD_DWD/"
# top_features <- 30
# meta_source_directory <- "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_manuscript/output/balanced_stratification/n2486/"
# metadata_file <- "n2486_complete_meta_w_onlysite_stratification.txt"

seed <- opt$seed
percent_train <- opt$percent_train
percent_test <- opt$percent_test
relapse_site1 <- opt$relapse_site_1
relapse_site2 <- opt$relapse_site_2
file <- opt$input_file
method <- opt$resampling_method
train_metric <- opt$metric
n_bootstraps <- opt$n_bootstraps
alpha_parameters <- as.numeric(unlist(str_split(opt$alphas, ",")))
alphas <- seq(alpha_parameters[1], alpha_parameters[2], by = alpha_parameters[3])
n_lambdas <- opt$n_lambdas
output_directory <- opt$output_directory
top_features <- opt$n_features

## Define variables
seed_name <- paste0(seed)
RNA <- case_when(file == "dwd" ~ "genes", 
                 file == "pub_modules" ~ "modules")
metric <- ifelse(train_metric == "Accuracy", "accuracy", 
                 ifelse(train_metric == "Kappa", "kappa", NA))
level1 <- ifelse(relapse_site1 == "NoRel", "No_Relapse", NA)
level2 <- case_when(relapse_site2 == "Distant" ~ "Distant_Relapse", 
                    relapse_site2 == "OnlyBrain" ~ "Brain_Relapse",
                    relapse_site2 == "OnlyBone" ~ "Bone_Relapse",
                    relapse_site2 == "OnlyLiver" ~ "Liver_Relapse",
                    relapse_site2 == "OnlyLung" ~ "Lung_Relapse",
                    relapse_site2 == "MultiMet" ~ "Multimet_Relapse", 
                    relapse_site2 == "Distant_notBone" ~ "Distant_Relapse")
comparison_name <- case_when(relapse_site1 == "NoRel" & relapse_site2 == "Distant" ~ "No_v_Distant",
                             relapse_site1 == "NoRel" & relapse_site2 == "Distant_notBone" ~ "No_v_Distant",
                             relapse_site1 == "NoRel" & relapse_site2 == "OnlyBone" ~ "No_v_OnlyBone", 
                             relapse_site1 == "NoRel" & relapse_site2 == "OnlyBrain" ~ "No_v_OnlyBrain", 
                             relapse_site1 == "NoRel" & relapse_site2 == "OnlyLiver" ~ "No_v_OnlyLiver",
                             relapse_site1 == "NoRel" & relapse_site2 == "OnlyLung" ~ "No_v_OnlyLung", 
                             relapse_site1 == "NoRel" & relapse_site2 == "MultiMet" ~ "No_v_MultiMet")
comparison_name_2 <- case_when(relapse_site1 == "NoRel" & relapse_site2 == "Distant_notBone" ~ "No_v_Distant_notBone")

test_train_name <- paste0(percent_train, 'Train_', percent_test, 'Test')
train_name <- paste0(percent_train, '% Training')
test_name <- paste0(percent_test, '% Testing')
method_metric <- paste(method, metric, sep='_')
data <- paste0(file, '-', RNA)
relapse_name <- paste0(relapse_site1, '-v-', relapse_site2)
gene_file <- paste0(file, '_file')
test_train_column <- if (grepl("any", relapse_name, ignore.case = T)) {
  paste0("Train", percent_train, "_Test", percent_test, '_', relapse_site2)
} else if (grepl("only", relapse_name, ignore.case = T)) {
  paste0('Train', percent_train, '_Test', percent_test)
} else if (grepl("MultiMet", relapse_name, ignore.case = T)) {
  paste0("Train", percent_train, "_Test", percent_test)
} else if (grepl("distant", relapse_site2, ignore.case = T)) {
  paste0("Train", percent_train, "_Test", percent_test)}

## Source files and metadata
source_directory <- case_when(
  file == "dwd" ~ "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/DWD_adjust/", 
  file == "pub_modules" ~ "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/calculate_modules/")
meta_source_directory <- opt$metadata_directory
metadata_file <- opt$metadata_file
dwd_file <- "n2486_Ensembl_dwd.txt"
pub_modules_file <- "n2486_Ensembl_dwd_935_published_modules_3_receptors.txt"

dir.create(paste0(output_directory, relapse_name), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_directory, relapse_name, "/ROC_plot"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_directory, relapse_name, "/feature_summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_directory, relapse_name, "/top_features_plot"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_directory, relapse_name, "/comparisons"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_directory, relapse_name, "/model_txt"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_directory, relapse_name, "/Youden"), recursive = TRUE, showWarnings = FALSE)

set.seed(seed)


## Prepare Data ----------------------------------------------------------------
## Load in metadata
setwd(meta_source_directory)
if (relapse_site1 %in% c("NoRel", "OtherRel")) { 
  
  # Based on comparison for model, exclude the appropriate samples:
  
  if (level1 == "No_Relapse" &&
    relapse_site2 %in% c("OnlyBone", "OnlyBrain", "OnlyLiver", "OnlyLung", "MultiMet")) {
  metadata_RNA <- read.table(paste0(meta_source_directory, metadata_file), sep = "\t", header = T) %>% 
    filter(.data[[comparison_name]] != "Other_Relapse")
  
  } else if (relapse_site2 == "Distant_notBone") {
    metadata_RNA <- read.table(paste0(meta_source_directory, metadata_file), sep = "\t", header = T) %>% 
      filter(Relapse_Site != "Bone")
    
  } else {metadata_RNA <- read.table(paste0(meta_source_directory, metadata_file), sep = "\t", header = T)}
  
  table(metadata_RNA[[test_train_column]], metadata_RNA[[comparison_name]])
  
  ## Load combined feature data
  matrix_RNA <- read.table(file.path(source_directory, get(gene_file)), sep = "\t", header = TRUE, row.names = 1)
  matrix_RNA <- t(as.matrix(matrix_RNA))
  matrix_RNA <- matrix_RNA[metadata_RNA$Sample_ID,]
  
  # Read this in again and save if needed for any site testing comparisons
  og_matrix_RNA <- read.table(file.path(source_directory, get(gene_file)), sep = "\t", header = TRUE, row.names = 1)
  og_matrix_RNA <- t(as.matrix(og_matrix_RNA))
  
  metadata_RNA[[comparison_name]] <- factor(metadata_RNA[[comparison_name]], levels=c(level1,level2))
  table(metadata_RNA[[comparison_name]], metadata_RNA[[test_train_column]])
  
  meta_training <- metadata_RNA %>% filter(.data[[test_train_column]] == 'Train')
  meta_testing <- metadata_RNA %>% filter(.data[[test_train_column]] == 'Test')
  
  ## Define training and testing sets
  train_x <- matrix_RNA[meta_training$Sample_ID,]
  train_y <- meta_training[[comparison_name]]
  test_x <- matrix_RNA[meta_testing$Sample_ID,]
  test_y <- meta_testing[[comparison_name]]
  
} else {print("relapse_site1 is not No_Relapse")}


## Train Elastic Net Model -----------------------------------------------------

## Define alpha values to use in the model
# The alpha parameter in the Elastic Net combines L1 and L2 regularization
# methods, with alpha = 1 equivalent to Lasso regression and alpha = 0
# equivalent to Ridge regression.
### Alpha values are set in the input file

## Define lambda values to use in the model
# The range of lambda values was determined by fitting binomial generalized
# linear models via penalized maximum likelihood (R package glmnet) with each
# alpha value on the training set.

# Initialize the number of lambda values (another regularization parameter) to
# test and creates an empty list to store these values.
### n_lambda are set in the input file
lambda <- list()

# Runs a loop over the set alpha values, initializes the EN model with each
# alpha, and stores the sequence of lambda values that glmnet generates
# internally for each model.

for(i in 1:length(alphas)) {
  init <- glmnet(train_x, train_y, alpha = alphas[i], family = 'binomial')
  lambda[[i]] <- init$lambda
}

# Determines the minimum and maximum lambda values across all models.
lambda_min <- min(unlist(lapply(lambda, min)))
lambda_max <- max(unlist(lapply(lambda, max)))

# Creates a grid of tuning parameters to be used in cross-validation, combining
# every pair of alpha and lambda values, creating 900 possible combinations.
tune_grid <- expand.grid(alpha = alphas,
                         lambda = seq(lambda_min, lambda_max, length.out = n_lambdas))

## Select the optimal tuning parameters to use
# To select the optimal tuning parameters: Bootstrapping will sample with replacement
train_control <- trainControl(method = method, number = n_bootstraps, classProbs = T)

# Trains the Elastic Net model using the caret package's unified interface,
# specifying 'Accuracy' or 'Kappa' as the evaluation metric and providing the
# tuning grid and training control object.
# N.B. Kappa is a better option when your dataset is unbalanced between groups
glmnet_obj <- train(train_x, train_y, method = "glmnet", metric = train_metric,
                    trControl = train_control, tuneGrid = tune_grid)

saveRDS(glmnet_obj, paste0(output_directory, relapse_name, "/EN_", test_train_name, "_", data, "_glmnet_obj_s", seed_name, "_", relapse_name, "_", method_metric, ".rds"))
glmnet_obj <- readRDS(paste0(output_directory, relapse_name, "/EN_", test_train_name, "_", data, "_glmnet_obj_s", seed_name, "_", relapse_name, "_", method_metric, ".rds"))
glmnet_obj
mean(glmnet_obj$results$Accuracy)

## Evaluating the trained model ------------------------------------------------

## Make Predictions ##
# Using the trained model, make predictions on train and test sets. 
# These predictions are probabilities of relapse, derived from RNA-seq features
pred_train <- predict(glmnet_obj, newdata = train_x, type = 'prob')
pred_test <- predict(glmnet_obj, newdata = test_x, type = 'prob')
pred_train <- prediction(pred_train[[level2]], labels = train_y,label.ordering = c(level1, level2))
pred_test <- prediction(pred_test[[level2]], labels = test_y,label.ordering = c(level1, level2))

## Assess Performance ##
# The script calculates the True Positive Rate (TPR) and False Positive Rate
# (FPR), essential components of the ROC (Receiver Operating Characteristic)
# curve, a graphical plot that illustrates the diagnostic ability of a binary
# classifier.
perf_train <- performance(pred_train, measure = 'tpr', x.measure = 'fpr')
perf_test <- performance(pred_test, measure = 'tpr', x.measure = 'fpr')

# The AUC (Area Under the ROC Curve) is computed, which is a comprehensive
# measure of the model's performance, with 1 indicating a perfect model and 0.5
# indicating a worthless model / no features selected.
auc_train <- signif(performance(pred_train, measure = 'auc')@y.values[[1]][1], 4)
auc_test <- signif(performance(pred_test, measure = 'auc')@y.values[[1]][1], 4)

# Plots the receiving operating curve as a way of evaluating the performance
# of a given elastic net model (2 curves).
plot_ROC2 <- function(perf1, perf2, a1, a2, set1, set2, plot_name) {
  png(paste0(output_directory, relapse_name, "/ROC_plot/", plot_name, ".png"), width = 1.5, height = 1.5, units = "in", res = 350)
  par(mai = c(0.24, 0.24, 0.05, 0.05), cex.axis = 0.3)
  plot(perf1, col = "red", lwd = 1.2, xlab = "", ylab = "", box.lwd=0.8,
       xaxis.xaxt = "n", yaxis.yaxt = "n")
  plot(perf2, col = "darkblue", lwd = 1.2, add = TRUE)
  abline(a=0, b=1, lwd = 0.8)
  axis(1, tck=(-0.02), lwd = 0.8)
  axis(2, tck=(-0.02), lwd = 0.8)
  mtext(side = 1, text = seq(0, 1, by=0.2), at = seq(0, 1, by=0.2), cex = 0.35, line = (-0.3))
  mtext(side = 2, text = seq(0, 1, by=0.2), at = seq(0, 1, by=0.2), cex = 0.35, line = 0.1)
  mtext(side = 2, text = "True positive rate", at = 0.5, cex = 0.4, line = 0.6)
  mtext(side = 1, text = "False positive rate", at = 0.5, cex = 0.4, line = 0.2)
  legend("bottomright", legend = c(paste0(set1, " (AUC = ", round(a1, 2), ")"), 
                                   paste0(set2, " (AUC = ", round(a2, 2), ")")), 
         lty = c(1,1), lwd = c(1,1),
         col = c('red','darkblue'), cex = 0.35, bty = "n")
  dev.off()
}

plot_ROC2_pdf <- function(perf1, perf2, a1, a2, set1, set2, plot_name) {
  pdf(paste0(output_directory, relapse_name, "/ROC_plot/", plot_name, ".pdf"), width = 1.5, height = 1.5)
  par(mai = c(0.24, 0.24, 0.05, 0.05), cex.axis = 0.3)
  plot(perf1, col = "red", lwd = 1.2, xlab = "", ylab = "", box.lwd=0.8,
       xaxis.xaxt = "n", yaxis.yaxt = "n")
  plot(perf2, col = "darkblue", lwd = 1.2, add = TRUE)
  abline(a=0, b=1, lwd = 0.8)
  axis(1, tck=(-0.02), lwd = 0.8)
  axis(2, tck=(-0.02), lwd = 0.8)
  mtext(side = 1, text = seq(0, 1, by=0.2), at = seq(0, 1, by=0.2), cex = 0.35, line = (-0.3))
  mtext(side = 2, text = seq(0, 1, by=0.2), at = seq(0, 1, by=0.2), cex = 0.35, line = 0.1)
  mtext(side = 2, text = "True positive rate", at = 0.5, cex = 0.4, line = 0.6)
  mtext(side = 1, text = "False positive rate", at = 0.5, cex = 0.4, line = 0.2)
  legend("bottomright", legend = c(paste0(set1, " (AUC = ", round(a1, 2), ")"), 
                                   paste0(set2, " (AUC = ", round(a2, 2), ")")), 
         lty = c(1,1), lwd = c(1,1),
         col = c('red','darkblue'), cex = 0.35, bty = "n")
  dev.off()
}

plot_ROC2(perf_train, perf_test, auc_train, auc_test, 
          train_name, test_name, paste0("EN_", test_train_name, "_", data, "_glmnet_s", seed_name, "_", relapse_name, "_", method_metric))

plot_ROC2_pdf(perf_train, perf_test, auc_train, auc_test, 
              train_name, test_name, paste0("EN_", test_train_name, "_", data, "_glmnet_s", seed_name, "_", relapse_name, "_", method_metric))

## Extract Model Coefficients --------------------------------------------------
# The script extracts the model's coefficients at the best lambda (determined
# during bootstrapping), providing insights into which features are most
# important in predicting site specific relapse.

if (RNA == "modules") {
  # Feature B-coefficients must be scaled if using modules 
  
  # Extract best alpha and lambda values
  best_alpha <- glmnet_obj$bestTune$alpha
  best_lambda <- glmnet_obj$bestTune$lambda
  
  # Extract coefficients for the best lambda value
  coef_matrix <- as.matrix(coef(glmnet_obj$finalModel, s = best_lambda))
  coef_matrix <- data.frame(
    feature = rownames(coef_matrix),
    unscaled_coef = as.numeric(coef_matrix),
    alpha = best_alpha,
    lambda = best_lambda)
  
  coef_matrix_no_intercept <- subset(coef_matrix, feature != "(Intercept)")
  
  # Re-scale B-coefficient to account for original scaling of features
  # This is especially important when including features like ER status, which 
  # is a binary 0/1 variable instead of a continuous scale like an RNA feature.
  coef_matrix_no_intercept$scaled_coef <- coef_matrix_no_intercept$unscaled_coef * 
    apply(train_x, 2, function(x) sqrt(sum((x - mean(x)) ^ 2) / length(train_y)))
  
  # Re-order features based on B-coefficient
  coef_matrix_no_intercept <- coef_matrix_no_intercept[order(coef_matrix_no_intercept$scaled_coef, decreasing = T), ]
  
  # Remove features that were not chosen in the model
  coef_matrix_no_intercept <- coef_matrix_no_intercept %>% 
    filter(scaled_coef != 0)
  
  write.table(coef_matrix_no_intercept,
              paste0(output_directory, relapse_name, "/model_txt/EN_", test_train_name, "_", data, "_glmnet_obj_s", seed_name, "_", relapse_name, "_", method_metric, ".txt"),
              sep = "\t", row.names=F)
  
  coef_matrix_no_intercept$coef <- coef_matrix_no_intercept$scaled_coef
  
} else if (RNA == "genes") {
  
  #  Extract best alpha and lambda values
  best_alpha <- glmnet_obj$bestTune$alpha
  best_lambda <- glmnet_obj$bestTune$lambda
  
  # Extract coefficients for the best lambda value
  coef_matrix <- as.matrix(coef(glmnet_obj$finalModel, s = best_lambda))
  coef_matrix <- data.frame(
    EnsemblID = rownames(coef_matrix), 
    unscaled_coef = as.numeric(coef_matrix), 
    alpha = best_alpha, 
    lambda = best_lambda)
  
  coef_matrix_no_intercept <- subset(coef_matrix, EnsemblID != "(Intercept)")
  
  # Re-order features based on B-coefficient
  coef_matrix_no_intercept <- coef_matrix_no_intercept[order(coef_matrix_no_intercept$unscaled_coef, decreasing = T), ]
  
  # Remove features that were not chosen in the model
  coef_matrix_no_intercept <- coef_matrix_no_intercept %>% 
    filter(unscaled_coef != 0)
  
  if (nrow(coef_matrix_no_intercept) > 0) {
    
  # Convert features from Ensembl IDs to gene symbols
  coef_matrix_no_intercept$GeneSymbol <- mapIds(
    org.Hs.eg.db, 
    keys = coef_matrix_no_intercept$EnsemblID, 
    column = "SYMBOL", # output is gene symbols
    keytype = "ENSEMBL", # providing EnsemblIDs
    multiVals = "first")
  
  write.table(coef_matrix_no_intercept,
              paste0(output_directory, relapse_name, "/model_txt/EN_", test_train_name, "_", data, "_glmnet_obj_s", seed_name, "_", relapse_name, "_", method_metric, ".txt"),
              sep = "\t", row.names=F)
  
  coef_matrix_no_intercept$feature <- coef_matrix_no_intercept$GeneSymbol
  coef_matrix_no_intercept$coef <- coef_matrix_no_intercept$unscaled_coef
  } else {
    write.table(coef_matrix_no_intercept,
                paste0(output_directory, relapse_name, "/model_txt/EN_", test_train_name, "_", data, "_glmnet_obj_s", seed_name, "_", relapse_name, "_", method_metric, ".txt"),
                sep = "\t", row.names=F)
  }
}


## Visualize top features chosen in the model ----------------------------------
  # Only generate these plots if features are chosen in the final model
if (nrow(coef_matrix_no_intercept) > 0) {
  
plot_data <- as.data.frame(coef_matrix_no_intercept)
plot_data <- plot_data %>%
  mutate(feature = as.character(feature)) %>% 
  mutate(coef = as.numeric(coef)) %>% 
  dplyr::select(feature, coef)

top_up <- plot_data %>%
  filter(coef > 0) %>%
  arrange(desc(coef)) %>%
  slice_head(n = top_features)

top_down <- plot_data %>%
  filter(coef < 0) %>%
  arrange(coef) %>%
  slice_head(n = top_features)

plot_data_plot <- bind_rows(top_up, top_down) %>% 
  mutate(feature = factor(feature, levels = unique(feature[order(coef, decreasing = FALSE)])))

p <- ggplot(plot_data_plot, aes(x = feature, y = coef, fill = coef > 0)) + 
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick")) + 
  coord_flip() + 
  labs(title = paste("Top", nrow(plot_data_plot), "out of", nrow(coef_matrix_no_intercept), "features selected in", relapse_site1, "vs.", relapse_site2, "EN model"), 
       x = "Feature", 
       y = case_when(
         RNA == "modules" ~ "Scaled B-coefficient", 
         RNA == "genes" ~ "B-coefficient")) + 
  theme_minimal()

png(paste0(output_directory, relapse_name, "/top_features_plot/EN_", relapse_name, "_top_features_", file, "_", test_train_name, "_", seed_name, ".png"),
    width = 10, height = 12, units = "in", res = 350)
print(p)
dev.off()

pdf(paste0(output_directory, relapse_name, "/top_features_plot/EN_", relapse_name, "_top_features_", file, "_", test_train_name, "_", seed_name, ".pdf"),
    width = 10, height = 12)
print(p)
dev.off()

## Youden's J cut point --------------------------------------------------------
# Determine the optimal cut point for the logistic regression model to determine
# sensitivity and specificity for model evaluation.

### Training Set
true_calls <- train_y # true outcomes / labels
probabilities <- predict(glmnet_obj, newdata = train_x, type = "prob")[, level2] # probabilities based on model
  
train_opt_cut <- cutpointr(x = probabilities, class = true_calls, direction = ">=",
                           pos_class = level2, neg_class = level1, 
                           method = maximize_metric, metric = youden)
  
png(paste0(output_directory, relapse_name, "/Youden/EN_opt_cut_", relapse_name, "_", file, "_", test_train_name, "_", seed_name, "_train.png"),
    width = 8, height = 4, units = 'in', res = 200)
plot(train_opt_cut)
dev.off()
  
sink(paste0(output_directory, relapse_name, "/Youden/EN_opt_cut_", relapse_name, "_", file, "_", test_train_name, "_", seed_name, "_train.txt"))
summary(train_opt_cut)
sink()
  
### Testing Set
true_calls <- test_y # true outcomes / labels
probabilities <- predict(glmnet_obj, newdata = test_x, type = "prob")[, level2] # probabilities based on model
  
test_opt_cut <- cutpointr(x = probabilities, class = true_calls, direction = ">=",
                          pos_class = level2, neg_class = level1, 
                          method = maximize_metric, metric = youden)
  
png(paste0(output_directory, relapse_name, "/Youden/EN_opt_cut_", relapse_name, "_", file, "_", test_train_name, "_", seed_name, "_test.png"), width = 8, height = 4, units = 'in', res = 200)
plot(test_opt_cut)
dev.off()
  
sink(paste0(output_directory, relapse_name, "/Youden/EN_opt_cut_", relapse_name, "_", file, "_", test_train_name, "_", seed_name, "_test.txt"))
summary(test_opt_cut)
sink()
}
  
## Generate summary file with all model parameters -----------------------------
y_int_value <- if (RNA == "modules") {
  coef_matrix %>% 
  filter(feature == "(Intercept)") %>% 
  dplyr::select(unscaled_coef) %>% 
  pull() 
} else if (RNA == "genes") {
  coef_matrix %>% 
    filter(EnsemblID == "(Intercept)") %>% 
    dplyr::select(unscaled_coef) %>% 
    pull()
}

model_coef_features <- nrow(coef_matrix_no_intercept)
if (nrow(coef_matrix_no_intercept) > 0) {
feature_list <- coef_matrix_no_intercept %>%
  dplyr::select(feature) %>%
  pull()
} else {
  feature_list <- NA_character_
}
feature_list <- paste(feature_list, collapse = ", ")
pos_bcoef_features <- coef_matrix_no_intercept$feature[coef_matrix_no_intercept$coef > 0]
pos_feature_list <- paste0(pos_bcoef_features, collapse = ', ')
neg_bcoef_features <- coef_matrix_no_intercept$feature[coef_matrix_no_intercept$coef < 0]
neg_feature_list <- paste0(neg_bcoef_features, collapse = ', ')

# Save parameters from this model
results_summary <- data.frame(
  Dataset = file,
  Input = RNA, 
  Training = train_name,
  Testing = test_name,
  Training_Method = method,
  Metric = metric, 
  Seed = seed_name,
  Ref_Site = level1, 
  Met_Site = relapse_site2,
  Training_AUC = signif(performance(pred_train, measure = 'auc')@y.values[[1]][1], 4), 
  Testing_AUC = signif(performance(pred_test, measure = 'auc')@y.values[[1]][1], 4), 
  Opt_Cut_Train = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(train_opt_cut$optimal_cutpoint, 6), NA_character_),
  Sensitivity_Train = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(train_opt_cut$sensitivity, 4), NA_character_),
  Specificity_Train = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(train_opt_cut$specificity, 4), NA_character_),
  Accuracy_Train = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(train_opt_cut$acc, 4), NA_character_),
  Opt_Cut_Test = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(test_opt_cut$optimal_cutpoint, 6), NA_character_),
  Sensitivity_Test = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(test_opt_cut$sensitivity, 4), NA_character_),
  Specificity_Test = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(test_opt_cut$specificity, 4), NA_character_),
  Accuracy_Test = ifelse(nrow(coef_matrix_no_intercept) > 0, signif(test_opt_cut$acc, 4), NA_character_),
  alpha = glmnet_obj$bestTune$alpha, 
  lambda = glmnet_obj$bestTune$lambda,
  y_int = y_int_value, 
  n_genes_selected = nrow(coef_matrix_no_intercept),
  n_pos_bcoef_feat_selected = length(coef_matrix_no_intercept$feature[coef_matrix_no_intercept$coef > 0]),
  n_neg_bcoef_feat_selected = length(coef_matrix_no_intercept$feature[coef_matrix_no_intercept$coef < 0]), 
  feature_names = feature_list,
  features_with_positive_bcoef = pos_feature_list, 
  features_with_negative_bcoef = neg_feature_list, 
  stringsAsFactors = FALSE
)

# Write results file
write.table(results_summary, paste0(output_directory, relapse_name, '/feature_summary/EN_', relapse_name, '_model_parameters_', file, '_', test_train_name, '_', seed_name, '.txt'), sep = '\t', row.names = F)  

## Evaluate performance of site specific models --------------------------------
model <- glmnet_obj

if (nrow(coef_matrix_no_intercept) > 0) {

### No v distant relapse -------------------------------------------------------
# For models predictive of distant relapse, plot the performance of the model
# for each site specific group separately:
if (relapse_site1 == "NoRel" & relapse_site2 %in% c("Distant", "Distant_notBone")) {
  
  site_list <- c("Bone", "Brain", "Liver", "Lung", "Multimet")
  results_list <- list()
  
  # Add in 75% training results
  metadata_train <- metadata_RNA %>% filter(.data[[test_train_column]] == "Train")
  train_counts <- table(metadata_train[[comparison_name]])
  
  results_list[[train_name]] <- data.frame(
    Group = train_name, 
    AUC = auc_train, 
    n_level1 = train_counts[level1],
    n_level2 = train_counts[level2],
    Total = (train_counts[level1] + train_counts[level2]))
  
  # Add in 25% testing results
  metadata_test <- metadata_RNA %>% filter(.data[[test_train_column]] == "Test")
  test_counts <- table(metadata_test[[comparison_name]])
  
  results_list[[test_name]] <- data.frame(
    Group = test_name, 
    AUC = auc_test, 
    n_level1 = test_counts[level1],
    n_level2 = test_counts[level2],
    Total = (test_counts[level1] + test_counts[level2]))
  
  for(site in site_list) {
    # Subset metadata test set for a specific site
    metadata_RNA_site <- read.table(paste0(meta_source_directory, metadata_file), sep = "\t", header = T) %>% 
      filter(Relapse_Site == site) %>% 
      filter(.data[[test_train_column]] == "Test")
    
    # Subset metadata test set for no_relapse
    metadata_RNA_noR <- read.table(paste0(meta_source_directory, metadata_file), sep = "\t", header = T) %>% 
      filter(.data[[comparison_name]] == "No_Relapse") %>% 
      filter(.data[[test_train_column]] == "Test")
    
    # Combine all sample metadata
    sub_metadata_RNA <- rbind(metadata_RNA_noR, metadata_RNA_site)
    
    # Subset feature matrix
    sub_matrix_RNA <- og_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    sub_metadata_RNA[[comparison_name]] <- factor(sub_metadata_RNA[[comparison_name]], levels = c(level1, level2))
    table(sub_metadata_RNA[[comparison_name]])
    
    test_x <- sub_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    test_y <- sub_metadata_RNA[[comparison_name]]
    
    # Sample counts per group
    sample_counts <- table(sub_metadata_RNA[[comparison_name]])
    n_level1 <- sample_counts[level1]
    n_level2 <- sample_counts[level2]
    
    # Make predictions #
    # Using the trained model, make prediction on this sub-test set
    pred_test <- predict(model, newdata = test_x, type = 'prob')
    pred_test <- prediction(pred_test[[level2]], labels = test_y,label.ordering = c(level1,level2))
    perf_test <- performance(pred_test, measure = 'tpr', x.measure = 'fpr')
    auc <- signif(performance(pred_test, measure = 'auc')@y.values[[1]][1], 4)
    
    results_list[[site]] <- data.frame(
      Group = paste0(site, " Only"), 
      AUC = auc, 
      n_level1 = n_level1, 
      n_level2 = n_level2, 
      Total = n_level1 + n_level2)
  }
  
  results_df <- do.call(rbind, results_list)
  
  results_df <- results_df %>% 
    mutate(sample_size_flag = ifelse(n_level2 > 10, n_level2, NA_character_))
  
  write.table(results_df,
              file = paste0(output_directory, relapse_name, "/comparisons/EN_", 
                            relapse_name, "_other_test_sets_", file, "_",
                            test_train_name, "_", seed_name, ".txt"),
              sep = "\t", col.names = T, row.names = F)

  plot_df <- results_df %>% 
    mutate(Group = factor(Group, levels = c("75% Training", "25% Testing", "Bone Only", "Brain Only", "Liver Only", "Lung Only", "Multimet Only"))) %>% 
    filter(!is.na(sample_size_flag))
  
  p <- ggplot(plot_df, aes(x = Group, y = AUC)) +
    geom_bar(stat = "identity") +
    geom_col(color = "black", fill = "gray80", linewidth = 1) +
    geom_text(aes(label = round(AUC,3)), vjust = -0.3, size = 4) + 
    theme_minimal() + 
    labs(x = "Group",
         y = "AUC", 
         title = paste0(comparison_name, " elastic net model using ", RNA), 
         subtitle = paste0("(seed ", seed_name, ")")) +
   # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(limits = c(0, 1))
  ggsave(filename = paste0(output_directory, relapse_name, "/comparisons/EN_", 
                           relapse_name, "_other_test_sets_", file, "_",
                           test_train_name, "_", seed_name, ".png"),
         plot = p, width = 5, height = 4, units = "in", dpi = 350)
  
  pdf(paste0(output_directory, relapse_name, "/comparisons/EN_", 
             relapse_name, "_other_test_sets_", file, "_",
             test_train_name, "_", seed_name, ".pdf"), height = 4, width = 5)
  print(p)
  dev.off()
  
### No v only site relapse -----------------------------------------------------
  # For models predictive of site specific relapse, plot performance on other sites:
} else if (relapse_site1 == "NoRel" & relapse_site2 %in%
           c("OnlyBone", "OnlyBrain", "OnlyLiver", "OnlyLung", "MultiMet")) {
  
  results_list <- list()
  
  # Add in 75% training results
  metadata_train <- metadata_RNA %>% filter(.data[[test_train_column]] == "Train")
  train_counts <- table(metadata_train[[comparison_name]])
  
  results_list[[train_name]] <- data.frame(
    Group = train_name, 
    AUC = auc_train, 
    n_level1 = train_counts[level1],
    n_level2 = train_counts[level2],
    Total = (train_counts[level1] + train_counts[level2]))
  
  # Add in 25% testing results
  metadata_test <- metadata_RNA %>% filter(.data[[test_train_column]] == "Test")
  test_counts <- table(metadata_test[[comparison_name]])
  
  results_list[[test_name]] <- data.frame(
    Group = test_name, 
    AUC = auc_test, 
    n_level1 = test_counts[level1],
    n_level2 = test_counts[level2],
    Total = (test_counts[level1] + test_counts[level2]))
  
  if (relapse_site2 == "OnlyBone") {
    site_list <- c("Brain", "Liver", "Lung", "Multimet")
  } else if (relapse_site2 == "OnlyBrain") {
    site_list <- c("Bone", "Liver", "Lung", "Multimet")
  } else if (relapse_site2 == "OnlyLiver") {
    site_list <- c("Bone", "Brain", "Lung", "Multimet")
  } else if (relapse_site2 == "OnlyLung") {
    site_list <- c("Bone", "Brain", "Liver", "Multimet")
  } else if (relapse_site2 == "MultiMet") {
    site_list <- c("Bone", "Brain", "Liver", "Lung")}
  
  for(site in site_list) {
    # Subset metadata set for a specific site (all samples)
    metadata_RNA_site <- read.table(paste0(meta_source_directory, metadata_file), sep = "\t", header = T) %>% 
      filter(Relapse_Site == site)
    
    # Subset metadata test set for no_relapse
    metadata_RNA_noR <- read.table(paste0(meta_source_directory, metadata_file), sep = "\t", header = T) %>% 
      filter(.data[[comparison_name]] == "No_Relapse") %>% 
      filter(.data[[test_train_column]] == "Test")
    
    # Combine all sample metadata
    sub_metadata_RNA <- rbind(metadata_RNA_noR, metadata_RNA_site)
    
    # Rename samples so that they correspond to levels in the model
    sub_metadata_RNA <- sub_metadata_RNA %>% 
      mutate(!!sym(comparison_name) := case_when(
        !!sym(comparison_name) == "Other_Relapse" ~ level2, 
        TRUE ~ !!sym(comparison_name)))
    
    # Subset feature matrix
    sub_matrix_RNA <- og_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    sub_metadata_RNA[[comparison_name]] <- factor(sub_metadata_RNA[[comparison_name]], levels = c(level1, level2))
    table(sub_metadata_RNA[[comparison_name]])
    
    test_x <- sub_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    test_y <- sub_metadata_RNA[[comparison_name]]
    
    # Sample counts per group
    sample_counts <- table(sub_metadata_RNA[[comparison_name]])
    n_level1 <- sample_counts[level1]
    n_level2 <- sample_counts[level2]
    
    # Make predictions #
    # Using the trained model, make prediction on this sub-test set
    pred_test <- predict(model, newdata = test_x, type = 'prob')
    pred_test <- prediction(pred_test[[level2]], labels = test_y,label.ordering = c(level1,level2))
    perf_test <- performance(pred_test, measure = 'tpr', x.measure = 'fpr')
    auc <- signif(performance(pred_test, measure = 'auc')@y.values[[1]][1], 4)
    
    results_list[[site]] <- data.frame(
      Group = paste0(site, " Only"), 
      AUC = auc, 
      n_level1 = n_level1, 
      n_level2 = n_level2, 
      Total = n_level1 + n_level2)
  }
  
  if (relapse_site2 != "MultiMet") {
    ### MultiMet with site of interest ###
    single_site <- case_when(
      relapse_site2 == "OnlyBone" ~ "Bone", 
      relapse_site2 == "OnlyBrain" ~ "Brain", 
      relapse_site2 == "OnlyLiver" ~ "Liver", 
      relapse_site2 == "OnlyLung" ~ "Lung")
    
    # Subset metadata to include MultiMet patients with site of interest only 
    metadata_RNA_site <- read.table(paste0(meta_source_directory, metadata_file), 
                                    sep = "\t", header = T) %>% 
      filter(Relapse_Site == "Multimet") %>% 
      filter(grepl(single_site, MultiMet_sites, ignore.case = T))
    
    # Subset metadata test set for no_relapse
    metadata_RNA_noR <- read.table(paste0(meta_source_directory, metadata_file), 
                                   sep = "\t", header = T) %>% 
      filter(.data[[comparison_name]] == "No_Relapse") %>% 
      filter(.data[[test_train_column]] == "Test")
    
    # Combine all sample metadata
    sub_metadata_RNA <- rbind(metadata_RNA_noR, metadata_RNA_site)
    
    # Rename samples so that they correspond to levels in the model
    sub_metadata_RNA <- sub_metadata_RNA %>% 
      mutate(!!sym(comparison_name) := case_when(
        !!sym(comparison_name) == "Other_Relapse" ~ level2, 
        TRUE ~ !!sym(comparison_name)))
    
    # Subset feature matrix
    sub_matrix_RNA <- og_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    sub_metadata_RNA[[comparison_name]] <- factor(sub_metadata_RNA[[comparison_name]], levels = c(level1, level2))
    table(sub_metadata_RNA[[comparison_name]])
    
    test_x <- sub_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    test_y <- sub_metadata_RNA[[comparison_name]]
    
    # Sample counts per group
    sample_counts <- table(sub_metadata_RNA[[comparison_name]])
    n_level1 <- sample_counts[level1]
    n_level2 <- sample_counts[level2]
    
    # Make predictions #
    # Using the trained model, make prediction on this sub-test set
    pred_test <- predict(model, newdata = test_x, type = 'prob')
    pred_test <- prediction(pred_test[[level2]], labels = test_y,label.ordering = c(level1,level2))
    perf_test <- performance(pred_test, measure = 'tpr', x.measure = 'fpr')
    auc <- signif(performance(pred_test, measure = 'auc')@y.values[[1]][1], 4)
    
    results_list$MultiMet_site <- data.frame(
      Group = paste0("MultiMet w/ ", single_site), 
      AUC = auc, 
      n_level1 = n_level1, 
      n_level2 = n_level2, 
      Total = n_level1 + n_level2)
    
    ### MultiMet without site of interest ###
    single_site <- case_when(
      relapse_site2 == "OnlyBone" ~ "Bone", 
      relapse_site2 == "OnlyBrain" ~ "Brain", 
      relapse_site2 == "OnlyLiver" ~ "Liver", 
      relapse_site2 == "OnlyLung" ~ "Lung")
    
    # Subset metadata to include MultiMet patients without site of interest only 
    metadata_RNA_site <- read.table(paste0(meta_source_directory, metadata_file), 
                                    sep = "\t", header = T) %>% 
      filter(Relapse_Site == "Multimet") %>% 
      filter(!grepl(single_site, MultiMet_sites, ignore.case = T))
    
    # Subset metadata test set for no_relapse
    metadata_RNA_noR <- read.table(paste0(meta_source_directory, metadata_file), 
                                   sep = "\t", header = T) %>% 
      filter(.data[[comparison_name]] == "No_Relapse") %>% 
      filter(.data[[test_train_column]] == "Test")
    
    # Combine all sample metadata
    sub_metadata_RNA <- rbind(metadata_RNA_noR, metadata_RNA_site)
    
    # Rename samples so that they correspond to levels in the model
    sub_metadata_RNA <- sub_metadata_RNA %>% 
      mutate(!!sym(comparison_name) := case_when(
        !!sym(comparison_name) == "Other_Relapse" ~ level2, 
        TRUE ~ !!sym(comparison_name)))
    
    # Subset feature matrix
    sub_matrix_RNA <- og_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    sub_metadata_RNA[[comparison_name]] <- factor(sub_metadata_RNA[[comparison_name]], levels = c(level1, level2))
    table(sub_metadata_RNA[[comparison_name]])
    
    test_x <- sub_matrix_RNA[sub_metadata_RNA$Sample_ID,]
    test_y <- sub_metadata_RNA[[comparison_name]]
    
    # Sample counts per group
    sample_counts <- table(sub_metadata_RNA[[comparison_name]])
    n_level1 <- sample_counts[level1]
    n_level2 <- sample_counts[level2]
    
    # Make predictions #
    # Using the trained model, make prediction on this sub-test set
    pred_test <- predict(model, newdata = test_x, type = 'prob')
    pred_test <- prediction(pred_test[[level2]], labels = test_y,label.ordering = c(level1,level2))
    perf_test <- performance(pred_test, measure = 'tpr', x.measure = 'fpr')
    auc <- signif(performance(pred_test, measure = 'auc')@y.values[[1]][1], 4)
    
    results_list$MultiMet_wo_site <- data.frame(
      Group = paste0("MultiMet w/o ", single_site), 
      AUC = auc, 
      n_level1 = n_level1, 
      n_level2 = n_level2,
      Total = n_level1 + n_level2)
  }  
  
  results_df <- do.call(rbind, results_list)
  results_df <- results_df %>% 
    mutate(sample_size_flag = ifelse(n_level2 >= 10, n_level2, NA_character_))
  
  write.table(results_df,
              file = paste0(output_directory, relapse_name, "/comparisons/EN_", 
                            relapse_name, "_other_test_sets_", file, "_",
                            test_train_name, "_", seed_name, ".txt"),
              sep = "\t", col.names = T, row.names = F)
  
  if (relapse_site2 == "OnlyBone") {
    plot_df <- results_df %>% 
      mutate(Group = factor(Group, 
                            levels = c("75% Training", "25% Testing", "Brain Only", "Liver Only", "Lung Only", "Multimet Only", "MultiMet w/ Bone", "MultiMet w/o Bone"))) %>% 
      filter(!is.na(sample_size_flag))
    
  } else if (relapse_site2 == "OnlyBrain") {
    plot_df <- results_df %>% 
      mutate(Group = factor(Group,
                            levels = c("75% Training", "25% Testing", "Bone Only", "Liver Only", "Lung Only", "Multimet Only", "MultiMet w/ Brain", "MultiMet w/o Brain"))) %>% 
      filter(!is.na(sample_size_flag))
    
  } else if (relapse_site2 == "OnlyLiver") {
    plot_df <- results_df %>% 
      mutate(Group = factor(Group, 
                            levels = c("75% Training", "25% Testing", "Bone Only", "Brain Only", "Lung Only", "Multimet Only", "MultiMet w/ Liver", "MultiMet w/o Liver"))) %>% 
      filter(!is.na(sample_size_flag))
    
  } else if (relapse_site2 == "OnlyLung") {
    plot_df <- results_df %>% 
      mutate(Group = factor(Group,
                            levels = c("75% Training", "25% Testing", "Bone Only", "Brain Only", "Liver Only", "Multimet Only", "MultiMet w/ Lung", "MultiMet w/o Lung"))) %>% 
      filter(!is.na(sample_size_flag))
    
  } else if (relapse_site2 == "MultiMet") {
    plot_df <- results_df %>% 
      mutate(Group = factor(Group, 
                            levels = c("75% Training", "25% Testing", "Bone Only", "Brain Only", "Liver Only", "Lung Only"))) %>% 
      filter(!is.na(sample_size_flag))
  }
  
  p <- ggplot(plot_df, aes(x = Group, y = AUC)) +
    geom_bar(stat = "identity") + 
    geom_col(color = "black", fill = "gray80", linewidth = 1) +
    geom_text(aes(label = round(AUC,3)), vjust = -0.3, size = 4) + 
    theme_minimal() + 
    labs(x = "Group",
         y = "AUC", 
         title = paste0(comparison_name, " elastic net model using ", RNA),
         subtitle = paste0("(seed ", seed_name, ")")) +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(limits = c(0, 1))
  ggsave(filename = paste0(output_directory, relapse_name, "/comparisons/EN_", 
                           relapse_name, "_other_test_sets_", file, "_",
                           test_train_name, "_", seed_name, ".png"),
         plot = p, width = 5, height = 4, units = "in", dpi = 350)
  
  pdf(paste0(output_directory, relapse_name, "/comparisons/EN_", 
             relapse_name, "_other_test_sets_", file, "_",
             test_train_name, "_", seed_name, ".pdf"), height = 4, width = 5)
  print(p)
  dev.off()
  
}
}
