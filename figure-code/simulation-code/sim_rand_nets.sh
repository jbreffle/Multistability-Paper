#!/bin/bash
#SBATCH -J genRandDist
#SBATCH --account=<your account here>
#SBATCH --partition=<your partition here>
#SBATCH --qos=low
#SBATCH --time=23:59:59
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<your email here>
#SBATCH --output=jobLogs/log_genRandDist-%A-%a.out
#SBATCH --error=jobLogs/err_genRandDist-%A-%a.out
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 4000 # Default units is MB
#SBATCH --exclude=compute-9-[0-7],compute-5-[1-3]

echo Starting job $SLURM_JOB_ID
echo Using node $SLURM_JOB_NODELIST

module load share_modules/MATLAB/R2019a

matlab -nodisplay -nodesktop -nosplash -nojvm -singleCompThread -r "sim_rand_nets; exit"
