#!/bin/bash
#SBATCH --job-name=dada2_run1%j
#SBATCH --output=dada2_run1%j.log
#SBATCH --error=dada2_run1%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user acauvin@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=60gb
#SBATCH --time=4-00:00:00
#SBATCH --account=juliemeyer
#SBATCH --qos=juliemeyer-b

#
pwd; hostname; date

cd $SLURM_SUBMIT_DIR

module load R

Rscript dada2_all.R