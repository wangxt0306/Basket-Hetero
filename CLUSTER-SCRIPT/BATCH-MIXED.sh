#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 36:00:00
#SBATCH --mem 4000
#SBATCH --output=./../CLUSTER-OUT/bhm-mixed-%a.out
#SBATCH --error=./../CLUSTER-ERR/bhm-mixed-%a.err
#SBATCH --array=1-648

## add R module
module add r/3.3.1

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../PROGRAMS/bhm-mixed.R ./../CLUSTER-LOG/bhm-mixed-$SLURM_ARRAY_TASK_ID.Rout