/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h" // for scanlink
#include <stdio.h>   // for NULL

/* Assemble eigenvectors, return bounds, etc. */
void mkeigvecs(struct scanlink *scanlist, /* linked list of fields to do with min ritz vals */
               double          *lambda,   /* ritz approximation to eigenvals of A */
               double          *bound,    /* on ritz pair approximations to eig pairs of A */
               int             *index,    /* the Ritz index of an eigenpair */
               double          *bj,       /* beta(j)*(last el. of corr. eigvec s of T) */
               int              d,        /* problem dimension = number of eigvecs to find */
               double          *Sres_max, /* Max value of Sres */
               double          *alpha,    /* vector of Lanczos scalars */
               double          *beta,     /* vector of Lanczos scalars */
               int              j,        /* number of Lanczos iterations taken */
               double          *s,        /* approximate eigenvector of T */
               double         **y,        /* columns of y are eigenvectors of A  */
               int              n,        /* problem size */
               double         **q         /* columns of q are Lanczos basis vectors */
)
{

  int              i, k;     /* indcies */
  double           Sres;     /* how well Tevec calculated eigvec s */
  struct scanlink *curlnk;   /* for traversing the scanlist */
  void             setvec(); /* initialize a vector */
  void             scadd();  /* add scalar multiple of vector to another */
  double           Tevec(double *, double *, int, double,
                         double *); /* calc eigenvector of T by linear recurrence */

  /* Scan for some data which is used here or later */
  i      = d;
  curlnk = scanlist;
  while (curlnk != NULL) {
    lambda[i] = curlnk->val;
    bound[i]  = bj[curlnk->indx];
    index[i]  = curlnk->indx;
    curlnk    = curlnk->pntr;
    i--;
  }

  /* Assemble the evecs from the Lanczos basis */
  for (i = 1; i <= d; i++) {
    Sres = Tevec(alpha, beta - 1, j, lambda[i], s);
    if (Sres > *Sres_max) {
      *Sres_max = Sres;
    }
    setvec(y[i], 1, n, 0.0);
    for (k = 1; k <= j; k++) {
      scadd(y[i], 1, n, s[k], q[k]);
    }
  }
}
