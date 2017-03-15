
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char *argv[]) 
{
  int comm_size, comm_rank, length;
  char proc_name[100] = "test";
  double result;

  comm_rank = 0;
  comm_size = 1;
  
  printf("Hello World from process = %d on processor %s\n", comm_rank, proc_name);
  result = exp((double) comm_rank);
  printf("Exp(%d) = %f\n", comm_rank, result);

  printf("Number of mpi processes = %d\n", comm_size);
} 
