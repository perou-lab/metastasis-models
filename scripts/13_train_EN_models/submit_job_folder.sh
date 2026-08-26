for job in $(ls /datastore/nextgenout4/share/labs/peroulab/Anna/metastatic_site_EN_github/scripts/13_train_EN_models/jobs/HR+HER2_only/genes_10814/*); do
  if [ -f "$job" ]; then
  echo "Submitting: $job"
  sbatch "$job"
  fi
done