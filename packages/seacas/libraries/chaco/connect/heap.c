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
#include <stdlib.h>
#include <stdio.h>
#include "structs.h"

#define left(i)		(2 * i)
#define right(i)	(2 * (i) + 1 )
#define parent(i)	((int) ((i)/2))


/* NOTE: heap assumes indices are 1-based. */

/* Note further: All the map arguments and manipulations allow me to */
/* quickly find a particular tag value.  This assumes that tags */
/* are integers between 0 and nvals. */

void 
heapify (
    struct heap *heap,		/* array of vals/tag to make into heap */
    int index,		/* root of subtree to heapify */
    int nvals,		/* number of values in array */
    int *map			/* maps from tag values to heap indices */
)
{
    double    swap_val;		/* temporary storage for swapping values */
    double    swap_tag;		/* temporary storage for swapping values */
    int       l, r;		/* indices of left & right children */
    int       largest;		/* index of largest value */

    l = left(index);
    r = right(index);

    if (l <= nvals && heap[l].val > heap[index].val)
	largest = l;
    else
	largest = index;

    if (r <= nvals && heap[r].val > heap[largest].val)
	largest = r;

    if (largest != index) {	/* swap index with largest and recurse */
	swap_val = heap[index].val;
	swap_tag = heap[index].tag;

	heap[index].val = heap[largest].val;
	heap[index].tag = heap[largest].tag;

	heap[largest].val = swap_val;
	heap[largest].tag = swap_tag;

	if (map != NULL) {	/* update pointers to heap tags */
	    map[heap[index].tag] = index;
	    map[heap[largest].tag] = largest;
	}

	heapify(heap, largest, nvals, map);
    }
}

/* Construct a heap from an unordered set of values. */
void 
heap_build (
    struct heap *heap,		/* array of vals/tag to make into heap */
    int nvals,		/* number of values in array */
    int *map			/* maps from tag values to heap indices */
)
{
    int       i;		/* loop counter */

    for (i = nvals / 2; i; i--) {
	heapify(heap, i, nvals, (int *) NULL);
    }

    if (map != NULL) {
	for (i = 1; i <= nvals; i++)
	    map[heap[i].tag] = i;
    }
}



double 
heap_extract_max (
    struct heap *heap,		/* array of vals/tag in a heap */
    int nvals,		/* number of values in array */
    int *ptag,			/* tag associated with return value */
    int *map			/* maps from tag values to heap indices */
)
{
    double    maxval;		/* return value */

    if (nvals < 1) {
	printf("Heap underflow\n");
	exit(0);
    }

    if (map != NULL) {		/* turn off map value for extracted tag */
	map[heap[1].tag] = 0;
    }

    maxval = heap[1].val;
    *ptag = heap[1].tag;

    heap[1].val = heap[nvals].val;
    heap[1].tag = heap[nvals].tag;

    if (map != NULL) {		/* update map value for root */
	map[heap[1].tag] = 1;
    }

    heapify(heap, 1, nvals - 1, map);

    return (maxval);
}


void 
heap_update_val (
    struct heap *heap,		/* array of vals/tag in a heap */
    int index,		/* index of value to update */
    double newval,		/* new value to insert */
    int *map			/* maps from tag values to heap indices */
)
{
    int       tag;		/* tag value associated with updated val */
    int       dad;		/* parent of a tree node */

    tag = heap[index].tag;

    dad = parent(index);
    while (index > 1 && heap[dad].val < newval) {
	heap[index].val = heap[dad].val;
	heap[index].tag = heap[dad].tag;
	if (map != NULL) {
	    map[heap[index].tag] = index;
	}
	index = dad;
	dad = parent(index);
    }

    heap[index].val = newval;
    heap[index].tag = tag;
    if (map != NULL) {
	map[tag] = index;
    }
}
