#!/bin/bash -l

#SBATCH -J basilisk_job
#SBATCH -o slurm.out
#SBATCH -e slurm.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=test

# Load all needed modules for Basilisk
module purge > /dev/null 2>&1
source load_modules.sh

# Run the file passed onto these parameters.
."$1"
