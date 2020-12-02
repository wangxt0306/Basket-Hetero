#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mem 5000
#SBATCH --output=./cluster-out/all-%a.out
#SBATCH --error=./cluster-err/all-%a.err
#SBATCH --array=1-200

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./programs/powers.R ./cluster-logs/power-$SLURM_ARRAY_TASK_ID.Rout
