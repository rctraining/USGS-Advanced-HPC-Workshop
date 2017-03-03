/*
 *  Copyright (C) 2015  Timothy Brown
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <mpi.h>

#ifdef __BIGGEST_ALIGNMENT__
#  define ALIGNMENT __BIGGEST_ALIGNMENT__
#else
#  define ALIGNMENT 64
#endif

int
main(int argc, char **argv)
{
	int ierr   = 0;			/* Error status                     */
	int rank   = 0;			/* MPI rank                         */
	int nprocs = 0;			/* MPI number of processors         */
	int i      = 0;			/* Indexer                          */
	int N      = 0;			/* Array size                       */
	int *data  = NULL;		/* Data array                       */
	int total  = 0;			/* Total of the data array          */

	/* Initialize MPI */
	ierr = MPI_Init(&argc, &argv);

	/* Find out our rank and total number of processors */
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	/*  Only the root rank reads the command line arguments */
	if (rank == 0) {
		N = strtol(argv[1], NULL, 10);
		if (N <= 0) {
			warnx("Will not run with negative array size.");
			MPI_Abort(MPI_COMM_WORLD, ierr);
		}
		printf("Running on %d processors.\nWith %d elements.\n", nprocs,
		       N);
	}

	/* Broadcast the data array size */
	ierr = MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* Allocate the data array */
	ierr = posix_memalign((void **)&data, ALIGNMENT, N * sizeof(int));
	if (ierr) {
		warnx("Unable to allocate: %lu", N * sizeof(int));
		MPI_Abort(MPI_COMM_WORLD, ierr);
	}

	/* Initialize the data array to contain my rank number */
	for (i = 0; i < N; ++i) {
		data[i] = rank;
	}

	/* All Reduce in-place with summation */
	ierr = MPI_Allreduce(MPI_IN_PLACE, &(data[0]), N, MPI_INT, MPI_SUM,
			     MPI_COMM_WORLD);

	/* Calculate our sum */
	for (i = 0; i < N; ++i) {
		total += data[i];
	}
	printf("[%02d] Total is: %d\n", rank, total);

	/* Free the data array */
	if (data) {
		free(data);
		data = NULL;
	}

	/* Finalize MPI and exit */
	ierr = MPI_Finalize();
	return(EXIT_SUCCESS);
}
