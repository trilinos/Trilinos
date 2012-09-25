/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#include <string.h>
#include "comm.h"
#include "zoltan_mem.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Given a list of processors & message data, sort them by proc ID            */
/* This routine will ensure that the ordering produced by the invert_map      */
/* routines is deterministic.  This should make bugs more reproducible.       */
/* This is a distribution count sort (Knuth)                                  */
/* An easily fixed assumption is that the smallest integer is zero!           */

int  Zoltan_Comm_Sort_Ints(
int *vals_sort,     /* values to be sorted */
int *vals_other,    /* other array to be reordered w/ sort */
int  nvals)         /* length of these two arrays */
{
    int *store=NULL, *copy_sort=NULL, *copy_other=NULL, *p;
    int i;
    int top;         /* largest integer to sort, smallest is 0 by assumption */
    int err = ZOLTAN_OK;
    int already_sorted = 1;  /* flag indicating whether vals_sort is
                                already sorted; can exit early and skip
                                memory allocations if it is.  */

    if (nvals < 1 || vals_sort == NULL  || vals_other == NULL)
       return ZOLTAN_FATAL;
    if (nvals == 1)
       return ZOLTAN_OK;           /* fastest way to sort 1 item is to return */
       
    /* find largest value (sort sometimes used for non processor lists) */   
    top = vals_sort[0];
    for (i = 1; i < nvals; i++) {
       if (vals_sort[i-1] > vals_sort[i]) already_sorted = 0;
       if (top < vals_sort[i]) top = vals_sort[i];
    }

    if (already_sorted)
       return ZOLTAN_OK;

    store      = (int*) ZOLTAN_CALLOC (top+2,  sizeof(int));
    copy_sort  = (int*) ZOLTAN_MALLOC (nvals * sizeof(int));
    copy_other = (int*) ZOLTAN_MALLOC (nvals * sizeof(int));

    if (store  &&  copy_sort  &&  copy_other)  {
       memcpy (copy_sort,  vals_sort,  nvals * sizeof(int));
       memcpy (copy_other, vals_other, nvals * sizeof(int));

       p = store+1;
       for (i = 0; i < nvals; i++)
          p[copy_sort[i]]++;                /* count number of occurances */

       for (i = 1; i < top+1; i++)
          p[i] += p[i-1];                   /* compute partial sums */
                                            /* assert: p[top] = nvals */

       p = store;                           /* effectively shifts down by one */
       for (i = 0; i < nvals; i++)  {
          vals_sort  [p[copy_sort [i]]] = copy_sort [i];
          vals_other [p[copy_sort [i]]] = copy_other[i];
          ++p[copy_sort [i]];
       }
    }
    else
       err =  ZOLTAN_MEMERR;
       
    Zoltan_Multifree (__FILE__, __LINE__, 3, &copy_sort, &copy_other, &store);
    return err;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
