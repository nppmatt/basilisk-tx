#!/usr/bin/env bash

#SBATCH -J basilisk_job
#SBATCH -o slurm-out/slurm.out
#SBATCH -e slurm-out/slurm.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=shahriar
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=4G

# Load all needed modules for Basilisk
echo "Loading modules."
module purge > /dev/null 2>&1
source scripts/load_modules.sh
echo "Modules loaded."

# Run the file passed onto these parameters.
echo "Passing $2 to $1"
./"$1" "$2"
#mpirun -n $SLURM_NTASKS ./"$1"

