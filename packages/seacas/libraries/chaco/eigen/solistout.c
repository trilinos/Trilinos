/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h" // for orthlink, orthlink_float
#include <stdio.h>   // for printf

/* Print out the orthogonalization set, double version */
void solistout(struct orthlink **solist, /* vector of pntrs to orthlnks */
               int               n,      /* length of vecs to orth. against */
               int               ngood,  /* number of good vecs on list */
               int               j       /* current number of Lanczos steps */
)
{
  int        i;           /* index */
  extern int DEBUG_EVECS; /* debugging output level for eigen computations */

  for (i = 1; i <= ngood; i++) {
    if ((solist[i])->index <= (j / 2)) {
      printf(".");
    }
    else {
      printf("+");
    }
    /* Really detailed output: printf("\n"); printf("depth
       %d\n",(solist[i])->depth); printf("index %d\n",(solist[i])->index);
       printf("ritzval %g\n",(solist[i])->ritzval); printf("betaji
       %g\n",(solist[i])->betaji); printf("tau %g\n",(solist[i])->tau);
       printf("prevtau %g\n",(solist[i])->prevtau);
       vecout((solist[i])->vec,1,n,"vec", NULL); */
  }
  printf("%d\n", ngood);

  if (DEBUG_EVECS > 2) {
    printf("  actual indices: ");
    for (i = 1; i <= ngood; i++) {
      printf(" %2d", solist[i]->index);
    }
    printf("\n");
  }
}

/* Print out the orthogonalization set, float version */
void solistout_float(struct orthlink_float **solist, /* vector of pntrs to orthlnks */
                     int                     n,      /* length of vecs to orth. against */
                     int                     ngood,  /* number of good vecs on list */
                     int                     j       /* current number of Lanczos steps */
)
{
  int i; /* index */

  for (i = 1; i <= ngood; i++) {
    if ((solist[i])->index <= (j / 2)) {
      printf(".");
    }
    else {
      printf("+");
    }
  }
  printf("%d\n", ngood);
}
