#!/bin/bash

#SBATCH -J basilisk_job
#SBATCH --output=%x.%j.out # %x.%j expands to slurm JobName.JobID
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=16384M

# Load all needed modules for Basilisk
echo "Loading modules."
module purge > /dev/null 2>&1
source load_modules.sh
echo "Modules loaded."

# Run the file passed onto these parameters.
echo "Running drop.c with arg vel: $1"
./out/drop $1

