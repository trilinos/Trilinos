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

#include <stdio.h>                      // for NULL
#include "structs.h"                    // for scanlink

/* Assemble eigenvectors, return bounds, etc. */
void 
mkeigvecs (
    struct scanlink *scanlist,      /* linked list of fields to do with min ritz vals */
    double *lambda,               /* ritz approximation to eigenvals of A */
    double *bound,                /* on ritz pair approximations to eig pairs of A */
    int *index,                /* the Ritz index of an eigenpair */
    double *bj,                   /* beta(j)*(last el. of corr. eigvec s of T) */
    int d,                    /* problem dimension = number of eigvecs to find */
    double *Sres_max,             /* Max value of Sres */
    double *alpha,		/* vector of Lanczos scalars */
    double *beta,			/* vector of Lanczos scalars */
    int j,			/* number of Lanczos iterations taken */
    double *s,			/* approximate eigenvector of T */
    double **y,                    /* columns of y are eigenvectors of A  */
    int n,                    /* problem size */
    double **q                    /* columns of q are Lanczos basis vectors */
)
{

    int       i,k;		/* indcies */
    double    Sres;             /* how well Tevec calculated eigvec s */
    struct scanlink *curlnk;    /* for traversing the scanlist */
    void      setvec();         /* initialize a vector */
    void      scadd();          /* add scalar multiple of vector to another */
    double    Tevec();          /* calc eigenvector of T by linear recurrence */

    /* Scan for some data which is used here or later */ 
    i = d;
    curlnk = scanlist;
    while (curlnk != NULL) {
        lambda[i] = curlnk->val;
        bound[i] = bj[curlnk->indx];
        index[i] = curlnk->indx;
        curlnk = curlnk->pntr;
        i--;
    }
 
    /* Assemble the evecs from the Lanczos basis */
    for (i = 1; i <= d; i++) {
        Sres = Tevec(alpha, beta - 1, j, lambda[i], s);
        if (Sres > *Sres_max) {
            *Sres_max = Sres;
        }
        setvec(y[i], 1, n, 0.0);
        for (k = 1; k <= j; k++) {
            scadd(y[i], 1, n, s[k], q[k]);
        }
    }
}
