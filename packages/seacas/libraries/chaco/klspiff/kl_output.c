/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
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
 *     * Neither the name of NTESS nor the names of its
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

#include "structs.h" // for bilist
#include <stdio.h>   // for printf, NULL

/*static void p1bucket();*/

void pbuckets(struct bilist ****buckets,   /* pointers to bucket lists */
              struct bilist **  listspace, /* elements within buckets */
              int               maxdeg,    /* maximum degree of a vertex */
              int               nsets      /* number of sets being divided into */
              )
{
  struct bilist *lptr; /* points to correct listspace */
  int            i, j; /* loop counter */
  void           p1bucket();

  printf("\n");
  for (i = 0; i < nsets; i++) {
    for (j = 0; j < nsets; j++) {
      if (i != j) {
        printf("For transition %d -> %d\n", i, j);
        if (j > i) {
          lptr = listspace[j - 1];
        }
        else {
          lptr = listspace[j];
        }
        p1bucket(buckets[i][j], lptr, maxdeg);
        printf("\n");
      }
    }
  }
  printf("\n");
}

/*static*/ void p1bucket(struct bilist **bucket, /* buckets holding bucket list */
                         struct bilist * lptr,   /* elements within bucket */
                         int             maxdeg  /* maximum degree of a vertex */
                         )
{
  struct bilist *bptr; /* loops through list at a bucket */
  int            val;  /* element in a bucket */
  int            size; /* array spacing */
  int            i;    /* loop counter */

  size = (int)(&(lptr[1]) - &(lptr[0]));
  for (i = 2 * maxdeg; i >= 0; i--) {
    if (bucket[i] != NULL) {
      printf("  Bucket %d:", i - maxdeg);
      for (bptr = bucket[i]; bptr != NULL; bptr = bptr->next) {
        val = ((int)(bptr - lptr)) / size;
        printf(" %d", val);
      }
      printf("\n");
    }
  }
}
