/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc_ret
#include "structs.h" // for bilist
#include <stdio.h>   // for NULL

int klv_init(struct bilist ***lbucket_ptr, /* space for left bucket sorts */
             struct bilist ***rbucket_ptr, /* space for right bucket sorts */
             struct bilist  **llistspace,  /* space for elements of linked lists */
             struct bilist  **rlistspace,  /* space for elements of linked lists */
             int            **ldvals,      /* change in separator for left moves */
             int            **rdvals,      /* change in separator for right moves */
             int              nvtxs,       /* number of vertices in the graph */
             int              maxchange    /* maximum change by moving a vertex */
)
{
  int sizeb; /* size of set of buckets */
  int sizel; /* size of set of pointers for all vertices */
  int flag;  /* return code */

  /* Allocate appropriate data structures for buckets, and listspace. */

  sizeb        = (2 * maxchange + 1) * sizeof(struct bilist *);
  *lbucket_ptr = smalloc_ret(sizeb);
  *rbucket_ptr = smalloc_ret(sizeb);

  *ldvals = smalloc_ret((nvtxs + 1) * sizeof(int));
  *rdvals = smalloc_ret((nvtxs + 1) * sizeof(int));

  sizel       = (nvtxs + 1) * sizeof(struct bilist);
  *llistspace = smalloc_ret(sizel);
  *rlistspace = smalloc_ret(sizel);

  if (*lbucket_ptr == NULL || *rbucket_ptr == NULL || *ldvals == NULL || *rdvals == NULL ||
      *llistspace == NULL || *rlistspace == NULL) {
    flag = 1;
  }

  else {
    flag = 0;
  }

  return (flag);
}
