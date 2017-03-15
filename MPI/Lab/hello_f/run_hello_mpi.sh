#!/bin/bash
#SBATCH -p normal
#SBATCH -t 0:02:00
#SBATCH -N 2
#SBATCH --job-name=mpi_hello
#SBATCH --output=out.mpi_hello
#SBATCH --ntasks-per-node=20

module load intel/psxe-2015

srun --mpi=pmi2 ./mpi_hello

