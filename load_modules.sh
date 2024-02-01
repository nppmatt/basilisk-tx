#!/bin/bash -l

echo "Loading easybuild, site, slurm/wulver"
module load easybuild site slurm/wulver

echo "Loading GCCcore, binutils, Bison"
module load GCCcore binutils Bison

echo "Loading GCC libunwind FFmpeg OpenMPI Mesa"
module load GCC libunwind FFmpeg OpenMPI Mesa

module list
