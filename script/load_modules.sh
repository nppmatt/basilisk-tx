#!/bin/bash -l

# 2024-05-19: Dependencies last pinned. Decided to keep it simple.

echo "Loading easybuild, site, slurm/wulver"
module load easybuild site slurm/wulver

echo "Loading GCCcore, binutils, Bison"
module load GCCcore/12.2.0 binutils/2.39 Bison/3.8.2

echo "Loading GCC libunwind FFmpeg OpenMPI Mesa intel"
module load GCC/12.2.0 libunwind/1.6.2 OpenMPI/4.1.6 Mesa/22.2.4

module list
