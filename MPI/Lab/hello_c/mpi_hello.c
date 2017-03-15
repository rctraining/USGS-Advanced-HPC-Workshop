
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char *argv[]) 
{
  int ierr, comm_size, comm_rank, length;
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  double result;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  ierr = MPI_Get_processor_name(proc_name, &length);

  printf("Hello World from process = %d on processor %s\n", comm_rank, proc_name);
  result = exp((double) comm_rank);
  printf("Exp(%d) = %f\n", comm_rank, result);

  ierr = MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank == 0) {
    printf("Number of mpi processes = %d\n", comm_size);
  }
  ierr = MPI_Finalize();
} 
