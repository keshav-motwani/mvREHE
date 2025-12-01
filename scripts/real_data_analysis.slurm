#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --array=0
#SBATCH --mail-type=END
#SBATCH --job-name=ncp
#SBATCH --cpus-per-task=110
#SBATCH --mem-per-cpu=2gb
#SBATCH --output=logs/real_data_analysis.log
#SBATCH --partition=gpu-l40
#SBATCH --account=stf

export OMP_NUM_THREADS=1

PER_TASK=1

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1))
END_NUM=$(( ($SLURM_ARRAY_TASK_ID) * $PER_TASK))

for (( PARAMETER_ID=$START_NUM; PARAMETER_ID<=END_NUM; PARAMETER_ID++ )); do
  echo "Running PARAMETER_ID=$PARAMETER_ID"
  Rscript scripts/real_data_analysis.R $PARAMETER_ID || continue
done
