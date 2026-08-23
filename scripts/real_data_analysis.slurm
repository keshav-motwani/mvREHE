#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=nc
#SBATCH --cpus-per-task=180
#SBATCH --mem-per-cpu=4gb
#SBATCH --array=1
#SBATCH --output=logs/222matrixmvrehe_ica_data_%A_%a.log
#SBATCH --partition=ckpt-all
#SBATCH --account=shojaiegpu
#SBATCH --requeue

export OMP_NUM_THREADS=1

/gscratch/stf/kmotwani/R/miniconda3/envs/R/bin/Rscript scripts/real_data_analysis.R
