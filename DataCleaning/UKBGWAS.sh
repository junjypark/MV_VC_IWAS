#!/bin/bash
#SBATCH --account=rrg-junpark  # replace this with your own account
#SBATCH --time=1:30:00
#SBATCH --job-name="GWAS"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80


module load intel/2019u4 r/4.1.2

cd /gpfs/fs0/scratch/j/junpark/tianyu47/UKB_GWAS

export R_LIBS=~/local/R_libs/

R CMD BATCH --no-save "--args $set" ${name}.R ${name}_set${set}.Rout


