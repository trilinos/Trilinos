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

#include	<stdio.h>
#include	"defs.h"


static void recursesort(double *vals, int nvals, int *indices, int *space);
static void merge(double *vals, int *indices, int length1, int length2, int *space);

void      mergesort(double *vals, int nvals, int *indices, int *space)
/* Merge sort values in vals, returning sorted indices. */
               			/* values to be sorted */
                		/* number of values */
                  		/* indices of sorted list */
                		/* space for nvals integers */
{
    extern int DEBUG_BPMATCH;	/* debugging flag for bipartite matching */
    int       flag;		/* has sorting screwed up? */
    int       i;

    for (i = 0; i < nvals; i++)
	indices[i] = i;

    recursesort(vals, nvals, indices, space);

    if (DEBUG_BPMATCH > 0) {
	flag = FALSE;
	for (i = 1; i < nvals; i++) {
	    if (vals[indices[i - 1]] > vals[indices[i]])
		flag = TRUE;
	}
	if (flag) {
	    printf("List improperly sorted in mergesort\n");
	    if (DEBUG_BPMATCH > 1) {
		for (i = 1; i < nvals; i++) {
		    printf("%d  %f\n", indices[i], vals[indices[i]]);
		}
	    }
	}
    }
}


static void recursesort(double *vals, int nvals, int *indices, int *space)
/* Recursive implementation of mergesort. */
               			/* values to be sorted */
                		/* number of values */
                  		/* indices of sorted list */
                		/* space for nvals integers */
{
    int       length1, length2;	/* lengths of two sublists */
    int       i;		/* temporary value */

    /* First consider base cases */
    if (nvals <= 1)
	return;

    else if (nvals == 2) {
	if (vals[indices[0]] <= vals[indices[1]])
	    return;
	else {
	    i = indices[0];
	    indices[0] = indices[1];
	    indices[1] = i;
	    return;
	}
    }

    else {
	length1 = nvals / 2;
	length2 = nvals - length1;
	recursesort(vals, length1, indices, space);
	recursesort(vals, length2, indices + length1, space);
	merge(vals, indices, length1, length2, space);
    }
}


static void merge(double *vals, int *indices, int length1, int length2, int *space)
/* Merge two sorted lists to create longer sorted list. */
               			/* values to be sorted */
                  		/* start of first sorted list */
                  		/* number of values in first list */
                  		/* number of values in second list */
                		/* sorted answer */
{
    int       n1, n2;		/* number of values seen in each list */
    int      *spaceptr;		/* loops through space array */
    int      *index1, *index2;	/* pointers into two lists */

    n1 = n2 = 0;
    spaceptr = space;
    index1 = indices;
    index2 = indices + length1;

    while (n1 < length1 && n2 < length2) {
	if (vals[*index1] <= vals[*index2]) {
	    *spaceptr++ = *index1++;
	    ++n1;
	}
	else {
	    *spaceptr++ = *index2++;
	    ++n2;
	}
    }

    /* Now add on remaining elements */
    while (n1 < length1) {
	*spaceptr++ = *index1++;
	++n1;
    }

    while (n2 < length2) {
	*spaceptr++ = *index2++;
	++n2;
    }

    spaceptr = space;
    index1 = indices;
    for (n1=length1+length2; n1; n1--) *index1++ = *spaceptr++;
}
