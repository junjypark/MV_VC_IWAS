#!/bin/bash
#SBATCH --account=rrg-junpark 
#SBATCH --time=1:00:00
#SBATCH --job-name="run_UKB_part2"
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=junjy.park@utoronto.ca
#SBATCH --mail-type=END
#SBATCH --output=%x.%a.out

THISJOBVALUE=$(( 0 + $SLURM_ARRAY_TASK_ID ))

module load CCEnv StdEnv/2020 gcc/9.3.0 r/4.1.0
cd /gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/files

export R_LIBS=~/local/R_libs/

srun Rscript run_IGAP_part2_F.R $SLURM_JOB_NAME $THISJOBVALUE 80 $SLURM_ARRAY_TASK_COUNT