#!/bin/bash
#SBATCH -p UV                               # Partition
#SBATCH -A training                         # Account
#SBATCH -t 0:01:00                          # Max walltime
#SBATCH -n 15                               # Number of requested cores
#SBATCH --job-name=hello                    # Job submission name
#SBATCH --output=hello.out                  # Output file name
#SBATCH --reservation=training_UV           # Reservation - will only work during workshop

module load intel

ifort -qopenmp omp_hello.f -o hello

echo "Running on the default number of cores"
./hello

echo "Running on 5 cores"
export OMP_NUM_THREADS=5
./hello
