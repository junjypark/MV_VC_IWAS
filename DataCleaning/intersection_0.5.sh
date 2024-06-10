#!/bin/bash
#SBATCH --account=rrg-junpark  # replace this with your own account
#SBATCH --time=1:30:00
#SBATCH --job-name="other"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80


module load intel/2019u4 r/4.1.2

cd /gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5

export R_LIBS=~/local/R_libs/

R CMD BATCH --no-save "--args $set" ${name}.R ${name}_set${set}.Rout


