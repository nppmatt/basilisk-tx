#!/bin/bash

#SBATCH -J basilisk_job
#SBATCH --partition=gpu
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:2
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=4000M

# Load all needed modules for Basilisk
echo "Loading modules."
module purge > /dev/null 2>&1
source load_modules.sh
echo "Modules loaded."

# Run the file passed onto these parameters.
echo "Running drop.c with args sigma: $1 | angle: $2 | vel: $3"
./out/drop $1 $2 $3

