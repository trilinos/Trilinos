/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc_ret
#include "structs.h" // for bilist
#include <stddef.h>  // for size_t
#include <stdio.h>   // for NULL

int kl_init(struct bilist *****bucket_ptrs, /* space for multiple bucket sorts */
            struct bilist ***  listspace,   /* space for all elements of linked lists */
            int ***            dvals,       /* change in cross edges for each move */
            int ***            tops,        /* top dval for each type of move */
            int                nvtxs,       /* number of vertices in the graph */
            int                nsets,       /* number of sets created at each step */
            int                maxchange    /* maximum change by moving a vertex */
)
{
  struct bilist * spacel; /* space for all listspace entries */
  struct bilist **spaceb; /* space for all buckets entries */
  size_t          sizeb;  /* size of set of buckets */
  size_t          sizel;  /* size of set of pointers for all vertices */
  int             i, j;   /* loop counters */
  double *        array_alloc_2D_ret(size_t dim1, size_t dim2, size_t size);

  /* Allocate appropriate data structures for buckets, and listspace. */

  *bucket_ptrs = (struct bilist ****)array_alloc_2D_ret(nsets, nsets, sizeof(struct bilist *));

  *dvals = (int **)array_alloc_2D_ret(nvtxs + 1, nsets - 1, sizeof(int));

  *tops = (int **)array_alloc_2D_ret(nsets, nsets, sizeof(int));

  /* By using '-1' in the next line, I save space, but I need to */
  /* be careful to get the right element in listspace each time. */
  *listspace = smalloc_ret((nsets - 1) * sizeof(struct bilist *));

  sizeb  = (2 * maxchange + 1) * sizeof(struct bilist *);
  sizel  = (nvtxs + 1) * sizeof(struct bilist);
  spacel = smalloc_ret((nsets - 1) * sizel);
  spaceb = smalloc_ret(nsets * (nsets - 1) * sizeb);

  if (*bucket_ptrs == NULL || *dvals == NULL || *tops == NULL || *listspace == NULL ||
      spacel == NULL || spaceb == NULL) {
    sfree(spacel);
    sfree(spaceb);
    return (1);
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

  return (0);
}
