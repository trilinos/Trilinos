/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h" // for orthlink, orthlink_float

void sorthog(double *          vec,    /* vector to be orthogonalized */
             int               n,      /* length of the columns of orth */
             struct orthlink **solist, /* set of vecs to orth. against */
             int               ngood   /* number of vecs in solist */
)
{
  double  alpha;
  double *dir;
  double  dot();
  void    scadd();
  int     i;

  for (i = 1; i <= ngood; i++) {
    dir   = (solist[i])->vec;
    alpha = -dot(vec, 1, n, dir) / dot(dir, 1, n, dir);
    scadd(vec, 1, n, alpha, dir);
  }
}

void sorthog_float(float *                 vec,    /* vector to be orthogonalized */
                   int                     n,      /* length of the columns of orth */
                   struct orthlink_float **solist, /* set of vecs to orth. against */
                   int                     ngood   /* number of vecs in solist */
)
{
  float  alpha;
  float *dir;
  double dot_float();
  void   scadd_float();
  int    i;

  for (i = 1; i <= ngood; i++) {
    dir   = (solist[i])->vec;
    alpha = -dot_float(vec, 1, n, dir) / dot_float(dir, 1, n, dir);
    scadd_float(vec, 1, n, alpha, dir);
  }
}
