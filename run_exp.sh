#!/bin/bash

vel=$1
sbatch --error="slurm_out/drop-$vel.err" --output="slurm_out/drop-$vel.out" vel_exp.sh $vel

