#!/bin/bash

#SBATCH --array=1-625
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=mvrehe
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2gb
#SBATCH --output=logs/mvrehe_simulation_%a.log
#SBATCH --partition=gpu-l40s
#SBATCH --account=stf

export OMP_NUM_THREADS=1

PER_TASK=2

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1))
END_NUM=$(( ($SLURM_ARRAY_TASK_ID) * $PER_TASK))

for (( PARAMETER_ID=$START_NUM; PARAMETER_ID<=END_NUM; PARAMETER_ID++ )); do
  echo "Running PARAMETER_ID=$PARAMETER_ID"
  Rscript scripts/simulation_runner.R $1 $2 $PARAMETER_ID || continue
done

