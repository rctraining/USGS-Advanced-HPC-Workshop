/*
 * Copyright (C) 2015  Timothy Brown
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <err.h>

#include "atts.h"


int
plus(double * restrict a, double * restrict b, double ** restrict c, int n)
{
	double *x = ASSUME_ALIGNED(a, ALIGNMENT);
	double *y = ASSUME_ALIGNED(b, ALIGNMENT);
	double *z = ASSUME_ALIGNED(*c, ALIGNMENT);
	int i = 0;

	for (i = 0; i < n; ++i) {
		z[i] = x[i] + y[i];
	}

	return(EXIT_SUCCESS);
}

int
main(int argc, char **argv)
{
	int ierr = 0;
	int i = 0;
	int n = 32;
	double *a = NULL;
	double *b = NULL;
	double *c = NULL;

	ierr = posix_memalign((void **)&a, ALIGNMENT,
			      n * sizeof(double));
	if (ierr) {
		err(ierr, "Unable to allocate: %ld", n * sizeof(double));
	}
	ierr = posix_memalign((void **)&b, ALIGNMENT,
			      n * sizeof(double));
	if (ierr) {
		err(ierr, "Unable to allocate: %ld", n * sizeof(double));
	}
	ierr = posix_memalign((void **)&c, ALIGNMENT,
			      n * sizeof(double));
	if (ierr) {
		err(ierr, "Unable to allocate: %ld", n * sizeof(double));
	}

	for (i = 0; i < n; ++i) {
		a[i] = 1;
		b[i] = 2;
		c[i] = 0;
	}

	ierr = plus(a, b, &c, n);
	if (ierr) {
		err(ierr, "Unable add vectors");
	}
	for (i = 0; i < n; ++i) {
		printf("%2d: a: %g\tb: %g\tc: %g\n",i, a[i], b[i], c[i]);
	}

	if (a) {
		free(a);
		a = NULL;
	}
	if (b) {
		free(b);
		b = NULL;
	}
	if (c) {
		free(c);
		c = NULL;
	}
	return(EXIT_SUCCESS);
}
