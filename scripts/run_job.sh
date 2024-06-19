#!/usr/bin/env bash

#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=shahriar
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-socket=64
#SBATCH --threads-per-core=1

# Load all needed modules for Basilisk
echo "Loading modules."
module purge > /dev/null 2>&1
source scripts/load_modules.sh
echo "Modules loaded."

# Run the file passed onto these parameters.
program_hash=$(cat $1 | md5sum | cut -c1-8)
config_hash=$(cat $2 | md5sum | cut -c1-8)

cpus="$3"
time_start=$(date +%s)
echo -e "[$(date +%D-%H:%M:%S)] START\t- $1 ($program_hash) $2 ($config_hash)\t- $cpus CPUs" >> log/job.log
#./"$1" "$2"
mpirun -n $SLURM_NTASKS ./"$1" "$2"

time_end=$(date +%s)
time_diff=$(echo "scale=2; ($time_end - $time_start) / 3600" | bc)
cost=$(echo "scale=1; $time_diff * $cpus" | bc)
echo -e "[$(date +%D-%H:%M:%S)] END  \t- $1 ($program_hash) $2 ($config_hash)\t- $cpus CPUs | $time_diff Hours | $cost SUs used" >> log/job.log

