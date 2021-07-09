/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"
#include <stdio.h>

void orthogonalize(double *         vec,     /* vector to be orthogonalized */
                   int              n,       /* length of the columns of orth */
                   struct orthlink *orthlist /* set of vectors to orthogonalize against */
)
{
  struct orthlink *curlnk;
  void             orthogvec();

  curlnk = orthlist;
  while (curlnk != NULL) {
    orthogvec(vec, 1, n, curlnk->vec);
    curlnk = curlnk->pntr;
  }
}
