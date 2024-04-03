#!/bin/bash -l

# 2024-04-03: GCC toolchain pinned to 12.2.0, libunwind 1.6.2
# binutils set to 2.39 to avoid slurm warning
# TODO: Need to push version pinning to config file
# TODO: Make build system check source whether Mesa or equivalent is needed

echo "Loading easybuild, site, slurm/wulver"
module load easybuild site slurm/wulver

echo "Loading GCCcore, binutils, Bison"
module load GCCcore/12.2.0 binutils/2.39 Bison

echo "Loading GCC libunwind FFmpeg OpenMPI Mesa intel"
module load GCC/12.2.0 libunwind/1.6.2 FFmpeg OpenMPI Mesa intel

module list
