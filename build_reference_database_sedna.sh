#!/bin/bash
#SBATCH --partition=node # Queue selection
#SBATCH --job-name=build_reference_database_sedna# Job name
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user= [EMAIL HERE] # Where to send mail\
#SBATCH --ntasks=1 # Run a single task\
#SBATCH --cpus-per-task=1  # Number of CPU cores per task\
#SBATCH --mem=10000 # Job memory request\
#SBATCH --time=24:00:00 # Time limit hrs:min:sec\
#SBATCH --output=build_reference_database.log # Standard output/error\
#export OMP_NUM_THREADS=8

module load R

#run dada2 script
Rscript build_reference_database_sedna.R
