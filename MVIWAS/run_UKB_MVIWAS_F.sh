#!/bin/bash
#SBATCH --account=rrg-junpark 
#SBATCH --time=1:00:00
#SBATCH --job-name="run_IGAP_MVIWAS_F"
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=yuant.tian@mail.utoronto.ca
#SBATCH --mail-type=END
#SBATCH --output=%x.%a.out

THISJOBVALUE=$(( 0 + $SLURM_ARRAY_TASK_ID ))

module load intel/2019u4 r/4.1.2
cd /gpfs/fs0/scratch/j/junpark/tianyu47/MV_VC_IWAS/MVIWAS

export R_LIBS=~/local/R_libs/

srun Rscript UKB_MVIWAS_F.R $SLURM_ARRAY_TASK_ID