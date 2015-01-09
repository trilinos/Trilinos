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

#include <stdio.h>                      // for NULL
#include "structs.h"                    // for bilist


/* Note: bi-directional lists aren't assumed to be sorted. */

void 
add2bilist (/* add val to unsorted list */
    struct bilist *lptr,		/* element to add */
    struct bilist **list		/* list added to */
)
{
    lptr->next = *list;
    if (*list != NULL)
	(*list)->prev = lptr;
    lptr->prev = NULL;
    *list = lptr;
}


void 
removebilist (
    struct bilist *lptr,		/* ptr to element to remove */
    struct bilist **list		/* head of list to remove it from */
)

/* Remove an element from a bidirectional list. */
{
    if (lptr->next != NULL)
	lptr->next->prev = lptr->prev;
    if (lptr->prev != NULL)
	lptr->prev->next = lptr->next;
    else
	*list = lptr->next;
}


void 
movebilist (
    struct bilist *lptr,		/* ptr to element to move */
    struct bilist **oldlist,	/* head of list to remove it from */
    struct bilist **newlist	/* head of list to add it to */
)

/* Move an element from a old bidirectional list to new one. */
{
    removebilist(lptr, oldlist);

    add2bilist(lptr, newlist);
}
