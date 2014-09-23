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

#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for NULL
#include "smalloc.h"                    // for smalloc_ret
#include "structs.h"                    // for bilist


int 
kl_init (
    struct bilist *****bucket_ptrs,	/* space for multiple bucket sorts */
    struct bilist ***listspace,	/* space for all elements of linked lists */
    int ***dvals,		/* change in cross edges for each move */
    int ***tops,			/* top dval for each type of move */
    int nvtxs,		/* number of vertices in the graph */
    int nsets,		/* number of sets created at each step */
    int maxchange		/* maximum change by moving a vertex */
)
{
    struct bilist *spacel;	/* space for all listspace entries */
    struct bilist **spaceb;	/* space for all buckets entries */
    size_t    sizeb;		/* size of set of buckets */
    size_t    sizel;		/* size of set of pointers for all vertices */
    int       i, j;		/* loop counters */
    double   *array_alloc_2D_ret(size_t dim1, size_t dim2, size_t size);

    /* Allocate appropriate data structures for buckets, and listspace. */

    *bucket_ptrs = (struct bilist ****)
       array_alloc_2D_ret(nsets, nsets, sizeof(struct bilist *));

    *dvals = (int **) array_alloc_2D_ret(nvtxs + 1, nsets - 1, sizeof(int));

    *tops = (int **) array_alloc_2D_ret(nsets, nsets, sizeof(int));

    /* By using '-1' in the next line, I save space, but I need to */
    /* be careful to get the right element in listspace each time. */
    *listspace = smalloc_ret((nsets - 1) * sizeof(struct bilist *));

    sizeb = (2 * maxchange + 1) * sizeof(struct bilist *);
    sizel = (nvtxs + 1) * sizeof(struct bilist);
    spacel = smalloc_ret((nsets - 1) * sizel);
    spaceb = smalloc_ret(nsets * (nsets - 1) * sizeb);

    if (*bucket_ptrs == NULL || *dvals == NULL || *tops == NULL ||
	*listspace == NULL || spacel == NULL || spaceb == NULL) {
	return(1);
    }

    for (i = 0; i < nsets; i++) {
	if (i != nsets - 1) {
	    (*listspace)[i] = spacel;
	    spacel += nvtxs + 1;
	}

	for (j = 0; j < nsets; j++) {
	    if (i != j) {
		(*bucket_ptrs)[i][j] = spaceb;
		spaceb += 2 * maxchange + 1;
	    }
	}
    }

    return(0);
}
