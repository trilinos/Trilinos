/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
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
