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
#include "smalloc.h"                    // for smalloc_ret
#include "structs.h"                    // for bilist


int 
klv_init (
    struct bilist ***lbucket_ptr,	/* space for left bucket sorts */
    struct bilist ***rbucket_ptr,	/* space for right bucket sorts */
    struct bilist **llistspace,	/* space for elements of linked lists */
    struct bilist **rlistspace,	/* space for elements of linked lists */
    int **ldvals,		/* change in separator for left moves */
    int **rdvals,		/* change in separator for right moves */
    int nvtxs,		/* number of vertices in the graph */
    int maxchange		/* maximum change by moving a vertex */
)
{
    int       sizeb;		/* size of set of buckets */
    int       sizel;		/* size of set of pointers for all vertices */
    int       flag;		/* return code */

    /* Allocate appropriate data structures for buckets, and listspace. */

    sizeb = (2 * maxchange + 1) * sizeof(struct bilist *);
    *lbucket_ptr = smalloc_ret(sizeb);
    *rbucket_ptr = smalloc_ret(sizeb);

    *ldvals = smalloc_ret((nvtxs + 1) * sizeof(int));
    *rdvals = smalloc_ret((nvtxs + 1) * sizeof(int));

    sizel = (nvtxs + 1) * sizeof(struct bilist);
    *llistspace = smalloc_ret(sizel);
    *rlistspace = smalloc_ret(sizel);

    if (*lbucket_ptr == NULL || *rbucket_ptr == NULL || *ldvals == NULL ||
	*rdvals == NULL || *llistspace == NULL || *rlistspace == NULL) {
	flag = 1;
    }

    else {
        flag = 0;
    }

    return(flag);
}
