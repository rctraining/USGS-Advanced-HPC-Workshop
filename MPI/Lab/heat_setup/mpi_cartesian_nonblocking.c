#include "mpi.h"
#include <stdio.h>
#define LEFT    0
#define RIGHT  1
#define UP 2
#define DOWN 3

int main(argc,argv)
int argc;
char *argv[];  {
int numtasks, rank, source, dest, outbuf, i;
int tag=1; 
int inbuf[4]={MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL};
int nbrs[4];
int dims[2] = {0,0};
int periods[2]={0,0}, reorder=0, coords[2];

MPI_Request reqs[8];
MPI_Status stats[8];
MPI_Comm cartcomm;

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

MPI_Dims_create(numtasks, 2, dims);
MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
MPI_Comm_rank(cartcomm, &rank);
MPI_Cart_coords(cartcomm, rank, 2, coords);
MPI_Cart_shift(cartcomm, 0, 1, &nbrs[LEFT], &nbrs[RIGHT]);
MPI_Cart_shift(cartcomm, 1, 1, &nbrs[UP], &nbrs[DOWN]);

printf("rank= %d coords= %d %d  neighbors(u,d,l,r)= %d %d %d %d\n",
        rank,coords[0],coords[1],nbrs[UP],nbrs[DOWN],nbrs[LEFT],
        nbrs[RIGHT]);

outbuf = rank;

for (i=0; i<4; i++) {
   dest = nbrs[i];
   source = nbrs[i];
   MPI_Isend(&outbuf, 1, MPI_INT, dest, tag, 
             MPI_COMM_WORLD, &reqs[i]);
   MPI_Irecv(&inbuf[i], 1, MPI_INT, source, tag, 
             MPI_COMM_WORLD, &reqs[i+4]);
   }

MPI_Waitall(8, reqs, stats);
   
printf("rank= %d                  inbuf(u,d,l,r)= %d %d %d %d\n",
       rank,inbuf[UP],inbuf[DOWN],inbuf[LEFT],inbuf[RIGHT]);
   
MPI_Finalize();
}