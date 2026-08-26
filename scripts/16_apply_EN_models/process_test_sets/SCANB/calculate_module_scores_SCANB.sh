#!/bin/bash
#SBATCH --job-name=SCANB_modules
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128g
#SBATCH --partition=allnodes
#SBATCH --output=/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/logs/SCANB_modules_%j.out
#SBATCH --error=/datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/logs/SCANB_modules_%j.err
#SBATCH --mail-user=anna_michmerhuizen@med.unc.edu
#SBATCH --mail-type=END,FAIL

# Load conda and activate environment
source /home/amich/conda/etc/profile.d/conda.sh
conda activate metastatic_site_EN

# Set environment variables
export R_HOME=/home/amich/conda/envs/metastatic_site_EN/lib/R
export LD_LIBRARY_PATH=/home/amich/conda/envs/metastatic_site_EN/lib/R/lib:/home/amich/conda/envs/metastatic_site_EN/lib:$LD_LIBRARY_PATH

# Print debug info
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "R_HOME=$R_HOME"
which Rscript
Rscript --version

# Run script
Rscript /datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/scripts/16_apply_EN_models/process_test_sets/SCANB/calculate_module_scores_SCANB.R