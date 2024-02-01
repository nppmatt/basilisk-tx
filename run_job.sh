#!/bin/bash -l

#SBATCH -J basilisk_job
#SBATCH -o out/slurm.out
#SBATCH -e slurm.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=general

echo "test"

# Load all needed modules for Basilisk
module purge > /dev/null 2>&1
source load_modules.sh
echo "Modules loaded."

# Run the file passed onto these parameters.
echo "Passing $1"
./"$1"
