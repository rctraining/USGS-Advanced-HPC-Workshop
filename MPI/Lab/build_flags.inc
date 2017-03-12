
## MPI+Fortran compiler
FC = mpif90   

## Optimization flags for gfortran
FCFLAGS = -O3 -mavx -ffixed-line-length-132

## Debuggin flags for gfortran
DBGFLAGS = -O0 -fbounds-check -fbacktrace -ffixed-line-length-132 
