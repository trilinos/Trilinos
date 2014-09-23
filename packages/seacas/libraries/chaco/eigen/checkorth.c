/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Governement retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>                       // for fabs, log10
#include <stdio.h>                      // for printf
#include "defs.h"                       // for min

/* Check orthogonality of vector set */
void 
checkorth (double **mat, int n, int dim)
{
    int       i, j;		/* loop idices */
    double    measure;		/* Froebenius norm */
    double    prod;		/* value of dot product */
    double    worst;		/* greatest off-diagonal dot product */
    int       lim;		/* index of last vec to check against */
    int       screenlim;	/* value of lim that will fit on screen */
    int       option;		/* which option to use */

    double    dot();		/* standard dot product routine */

    /* The T/F argument in the conditionals is just a convenient option: */

    screenlim = 20;
    option = 3;

    /* Check orthogonality over whole set. */
    if (option == 1) {
	printf("Orthogonality check:\n");
	for (i = 1; i <= dim; i++) {
	    printf("%2d)", i);
	    for (j = 1; j <= i; j++) {
		prod = dot(mat[i], 1, n, mat[j]);
		/* printf(" %g ",prod); */
		/* printf(" %4.2e ",prod); */
		/* printf(" %4.2e ",fabs(prod)); */
		printf(" %2d", -(int) log10(prod));
	    }
	    printf("\n");
	}
    }

    if (option == 2) {
	printf("Frobenius orthogonality measure:");
	measure = 0;
	for (i = 1; i <= dim; i++) {
	    for (j = i; j <= dim; j++) {
		prod = dot(mat[i], 1, n, mat[j]);
		if (i == j) {
		    measure += fabs(1.0 - prod);
		}
		else {
		    measure += 2.0 * fabs(prod);
		}
	    }
	}
	printf("%g \n", measure);
    }

    /* Check orthogonality against last vector. Allows you to build up orthogonality
       matrix much faster if previous columns stay the same when add a new column,
       but may interact with other debug output to give a confusing presentation. */
    if (option == 3) {
	printf("%3d) ", dim);
	lim = min(dim, screenlim);
	worst = 0;
	for (i = 1; i <= dim; i++) {
	    prod = dot(mat[i], 1, n, mat[dim]);
	    if (i <= lim) {
		printf(" %2d", -(int) log10(fabs(prod)));
	    }
	    if ((i != dim) && (fabs(prod) > fabs(worst))) {
		worst = prod;
	    }
	}
	printf(" worst %4.2e\n", worst);
    }
}


/* Check orthogonality of vector set */
void 
checkorth_float (float **mat, int n, int dim)
{
    int       i, j;		/* loop idices */
    double    measure;		/* Froebenius norm */
    double    prod;		/* value of dot product */
    double    worst;		/* greatest off-diagonal dot product */
    int       lim;		/* index of last vec to check against */
    int       screenlim;	/* value of lim that will fit on screen */
    int       option;		/* which option to use */

    double    dot_float();	/* standard dot product routine */

    /* The T/F argument in the conditionals is just a convenient option: */

    screenlim = 20;
    option = 3;

    /* Check orthogonality over whole set. */
    if (option == 1) {
	printf("Orthogonality check:\n");
	for (i = 1; i <= dim; i++) {
	    printf("%2d)", i);
	    for (j = 1; j <= i; j++) {
		prod = dot_float(mat[i], 1, n, mat[j]);
		/* printf(" %g ",prod); */
		/* printf(" %4.2e ",prod); */
		/* printf(" %4.2e ",fabs(prod)); */
		printf(" %2d", -(int) log10(prod));
	    }
	    printf("\n");
	}
    }

    if (option == 2) {
	printf("Frobenius orthogonality measure:");
	measure = 0;
	for (i = 1; i <= dim; i++) {
	    for (j = i; j <= dim; j++) {
		prod = dot_float(mat[i], 1, n, mat[j]);
		if (i == j) {
		    measure += fabs(1.0 - prod);
		}
		else {
		    measure += 2.0 * fabs(prod);
		}
	    }
	}
	printf("%g \n", measure);
    }

    /* Check orthogonality against last vector. Allows you to build up orthogonality
       matrix much faster if previous columns stay the same when add a new column,
       but may interact with other debug output to give a confusing presentation. */
    if (option == 3) {
	printf("%3d) ", dim);
	lim = min(dim, screenlim);
	worst = 0;
	for (i = 1; i <= dim; i++) {
	    prod = dot_float(mat[i], 1, n, mat[dim]);
	    if (i <= lim) {
		printf(" %2d", -(int) log10(fabs(prod)));
	    }
	    if ((i != dim) && (fabs(prod) > fabs(worst))) {
		worst = prod;
	    }
	}
	printf(" worst %4.2e\n", worst);
    }
}

