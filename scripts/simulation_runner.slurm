#!/bin/bash

#SBATCH --constraint=intel
#SBATCH --array=1-1680%500
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad
#SBATCH --job-name=test
#SBATCH --mem-per-cpu=6gb
#SBATCH --time 500
#SBATCH --output=logs/simulation_%a.log

export OMP_NUM_THREADS=1

ml R/4.1

PER_TASK=5

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

for (( PARAMETER_ID=$START_NUM; PARAMETER_ID<=END_NUM; PARAMETER_ID++ )); do
  timeout 1h Rscript --vanilla mvREHE.R $PARAMETER_ID
done

