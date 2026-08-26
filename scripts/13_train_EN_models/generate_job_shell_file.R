# script to generate .sh files and submit jobs for 1 (or more) seeds for EN models

library(tidyverse)

set.seed(123)

seed_list <- 124022
# seed_list <- sample(1:1e6, 20)
percent_train <- "75"
percent_test <- "25"

# Options for feature input
file_input <- c("pub_modules", "dwd")

# Options for site1
relapse_site1 <- c("NoRel")

# Options for site2, dependent upon value of site1
relapse_site_NoRel <- c("Distant", "Distant_notBone",
                        "OnlyBone", "OnlyBrain", "OnlyLiver", "OnlyLung", "MultiMet")

# Map relapse_site1 values to their corresponding site2 options
relapse_site_map <- list(NoRel = relapse_site_NoRel)

# Map metadata files and locations for subgroup analyses
group_list <- c("n2486", "HR+HER2-", "HER2+", "TNBC",
                "basal", "lumA", "lumB", "lumA_lumB", "HER2E")

metadata_map <- tibble(
  group = group_list,
  metadata_directory = c(
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/n2486/",
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/HR+HER2-/",
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/HER2+/",
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/TNBC/",
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/basal_only/",
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/lumA_only/",
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/lumB_only/",
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/lumA_lumB_only/", 
    "/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/balanced_stratification/HER2E_only/"),
  metadata_file = c(
    "n2486_complete_meta_w_onlysite_stratification.txt",
    "HR+HER2-_only_complete_meta_w_onlysite_stratification.txt",
    "HER2+_only_complete_meta_w_onlysite_stratification.txt",
    "TNBC_only_complete_meta_w_onlysite_stratification.txt",
    "basal_only_complete_meta_w_onlysite_stratification.txt",
    "lumA_only_complete_meta_w_onlysite_stratification.txt",
    "lumB_only_complete_meta_w_onlysite_stratification.txt",
    "lumA_lumB_only_complete_meta_w_onlysite_stratification.txt", 
    "HER2E_only_complete_meta_w_onlysite_stratification.txt"))

# Generate combinations of relapse site(s) for models
combinations <- expand.grid(features = file_input,
                            relapse_site1 = relapse_site1, 
                            relapse_site2 = relapse_site_NoRel,
                            seed = seed_list,
                            group = group_list,
                            stringsAsFactors = FALSE)

combinations <- combinations %>%
  left_join(metadata_map, by = "group")

for (i in 1:nrow(combinations)) {
  
  # Extract values for the current job:
  features <- combinations$features[i]
  feature_folder <- case_when(
    features == "dwd" ~ "genes_10814", 
    features == "pub_modules" ~ "modules_935+3")
  
  seed <- combinations$seed[i]
  relapse_site1 <- combinations$relapse_site1[i]
  relapse_site2 <- combinations$relapse_site2[i]
  group <- combinations$group[i]
  metadata_directory <- combinations$metadata_directory[i]
  metadata_file <- combinations$metadata_file[i]
  
  group_label <- case_when(
    group == "n2486" ~ "n2486", 
    group == "HR+HER2-" ~ "HR+HER2_only", 
    group == "HER2+" ~ "HER2+_only", 
    group == "TNBC" ~ "TNBC_only", 
    group == "basal" ~ "basal_only", 
    group == "lumA" ~ "lumA_only", 
    group == "lumB" ~ "lumB_only", 
    group == "HER2E" ~ "HER2E_only", 
    group == "lumA_lumB" ~ "lumA_lumB_only"
  )
  
  job_name <- paste0(relapse_site1, "-v-", relapse_site2, "_", seed, "_", features)
  file_name <- paste0("/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/scripts/13_train_EN_models/jobs/", group_label, "/", feature_folder, "/", relapse_site1, "_v_", relapse_site2, "_n2486_run_", seed, ".sh")
  dir.create(dirname(file_name), recursive = TRUE, showWarnings = FALSE)
  
  # Generate shell script for the job:
  sink(file_name)
  cat(
    "#!/bin/bash\n",
    "#SBATCH --job-name=", job_name, "\n",
    "#SBATCH --time=48:00:00\n",
    "#SBATCH --cpus-per-task=1\n",
    "#SBATCH --mem=16g\n",
    "#SBATCH --partition=allnodes\n",
    "#SBATCH --output=/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/logs/EN_models/", group_label, "/", job_name, "_%j.out\n",
    "#SBATCH --error=/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/logs/EN_models/", group_label, "/", job_name, "_%j.err\n",
    "#SBATCH --mail-user=anna_michmerhuizen@med.unc.edu\n",
    "#SBATCH --mail-type=FAIL\n\n",
    
    "# Load conda and activate environment\n",
    "source /home/amich/conda/etc/profile.d/conda.sh\n",
    "conda activate metastatic_site_EN\n\n",
    
    "# Set environment variables\n",
    "export R_HOME=/home/amich/conda/envs/metastatic_site_EN/lib/R\n",
    "export LD_LIBRARY_PATH=/home/amich/conda/envs/metastatic_site_EN/lib/R/lib:/home/amich/conda/envs/metastatic_site_EN/lib:$LD_LIBRARY_PATH\n\n",
    
    "# Print debug info\n",
    "echo \"LD_LIBRARY_PATH=$LD_LIBRARY_PATH\"\n",
    "echo \"R_HOME=$R_HOME\"\n",
    "which Rscript\n",
    "Rscript --version\n\n",
    
    "# Run script\n",
    "Rscript /datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/scripts/13_train_EN_models/train_EN_site_AM.R \\\n",
    "  --seed ", seed, " \\\n",
    "  --percent_train ", percent_train, " \\\n",
    "  --percent_test ", percent_test, " \\\n",
    "  --relapse_site_1 ", relapse_site1, " \\\n",
    "  --relapse_site_2 ", relapse_site2, " \\\n",
    "  --input_file ", features, " \\\n",
    "  --resampling_method boot \\\n",
    "  --metric Accuracy \\\n",
    "  --n_bootstraps 1000 \\\n",
    "  --alphas 0.1,0.9,0.1 \\\n",
    "  --n_lambdas 100 \\\n",
    "  --n_features 30 \\\n",
    "  --output_directory /datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/output/train_EN_models/", group_label, "/", feature_folder,"/ \\\n",
    "  --metadata_directory ", metadata_directory, " \\\n",
    "  --metadata_file ", metadata_file, " \n",
    sep = ""
  )
  
  sink()
}
