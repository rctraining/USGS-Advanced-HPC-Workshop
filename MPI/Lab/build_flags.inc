####################################
#  GFORTRAN FLAGS
## MPI+Fortran compiler
FC = mpif90   

## Optimization flags for gfortran
FCFLAGS = -O3 -mavx -ffixed-line-length-132

## Debugging flags for gfortran
DBGFLAGS = -O0 -fbounds-check -fbacktrace -ffixed-line-length-132 

####################################
#  IFORRT FLAGS
## MPI+Fortran compiler
##FC = mpiifort

## Optimization flags for Intel
FCFLAGS = -O3 ###-xCORE-AVX2 ##-ffixed-line-length-132

## Debugging flags for Intel
DBGFLAGS = -O0 -CB -g -traceback ##-ffixed-line-length-132 

