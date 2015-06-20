/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
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

#include	<math.h>
#include	"defs.h"


/* Find eigenvalues of 2x2 symmetric system by solving quadratic. */
void 
evals2 (
    double H[2][2],		/* symmetric matrix for eigenvalues */
    double *eval1,		/* smallest eigenvalue */
    double *eval2		/* middle eigenvalue */
)
{
    double    M[2][2];		/* normalized version of matrix */
    double    b, c;		/* coefficents of cubic equation */
    double    root1, root2;	/* roots of quadratic */
    double    xmax;		/* largest matrix element */
    int       i, j;		/* loop counters */

    M[0][0]=M[1][0]=M[0][1]=M[1][1]=0.0;
    
    xmax = 0.0;
    for (i = 0; i < 2; i++) {
	for (j = i; j < 2; j++) {
	    if (fabs(H[i][j]) > xmax)
		xmax = fabs(H[i][j]);
	}
    }
    if (xmax != 0) {
	for (i = 0; i < 2; i++) {
	    for (j = 0; j < 2; j++)
		M[i][j] = H[i][j] / xmax;
	}
    }


    b = -M[0][0] - M[1][1];
    c = M[0][0] * M[1][1] - M[1][0] * M[1][0];
    root1 = -.5 * (b + sign(b) * sqrt(b * b - 4 * c));
    root2 = c / root1;

    root1 *= xmax;
    root2 *= xmax;
    *eval1 = min(root1, root2);
    *eval2 = max(root1, root2);
}


/* Solve for eigenvector of SPD 2x2 matrix, with given eigenvalue. */
void 
eigenvec2 (
    double A[2][2],		/* matrix */
    double eval,			/* eigenvalue */
    double evec[2],		/* eigenvector returned */
    double *res			/* normalized residual */
)
{
    double    norm;		/* norm of eigenvector */
    double    res1, res2;	/* components of residual vector */
    int       i;		/* loop counter */

    if (fabs(A[0][0] - eval) > fabs(A[1][1] - eval)) {
	evec[0] = -A[1][0];
	evec[1] = A[0][0] - eval;
    }
    else {
	evec[0] = A[1][1] - eval;
	evec[1] = -A[1][0];
    }

    /* Normalize eigenvector and calculate a normalized eigen-residual. */
    norm = sqrt(evec[0] * evec[0] + evec[1] * evec[1]);
    if (norm == 0) {
	evec[0] = 1;
	evec[1] = 0;
	norm = 1;
    }
    for (i = 0; i < 2; i++)
	evec[i] /= norm;
    res1 = (A[0][0] - eval) * evec[0] + A[1][0] * evec[1];
    res2 = A[1][0] * evec[0] + (A[1][1] - eval) * evec[1];
    *res = sqrt(res1 * res1 + res2 * res2);

    res1 = fabs(A[0][0]) + fabs(A[1][0]);
    res2 = fabs(A[1][1]) + fabs(A[1][0]);
    *res /= max(res1, res2);
}
