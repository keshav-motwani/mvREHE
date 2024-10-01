#!/bin/bash

#SBATCH --array=1-465
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=test
#SBATCH --mem-per-cpu=16gb
#SBATCH --output=logs/simulation_%a.log

export OMP_NUM_THREADS=1

ml R
ml gemma

PER_TASK=10

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

for (( PARAMETER_ID=$START_NUM; PARAMETER_ID<=END_NUM; PARAMETER_ID++ )); do
  Rscript --vanilla scripts/simulation_runner.R $1 $2 $PARAMETER_ID || continue
done

