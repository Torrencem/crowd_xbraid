#!/bin/bash
#SBATCH -N 1
#SBATCH -t 0:30:00
#SBATCH -p pbatch
#SBATCH -o schur_%j.out
#SBATCH -e schur_%j.err
srun -n36 ./crowd_horesh_xbraid_schur
