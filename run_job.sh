#!/bin/bash

#SBATCH -J basilisk_job
#SBATCH -o slurm_out/slurm.out
#SBATCH -e slurm_out/slurm.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --time=00:20:00

# Load all needed modules for Basilisk
echo "Loading modules."
module purge > /dev/null 2>&1
source load_modules.sh
echo "Modules loaded."

# Run the file passed onto these parameters.
echo "Passing $1"
./"$1"
