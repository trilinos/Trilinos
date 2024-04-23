/*
 * Copyright(C) 1999-2020, 2022, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for ipairs
#include <stdio.h>   // for NULL

static struct ipairs *pedges; /* perturbed edges */
static double        *pvals;  /* perturbed values */

/* Initialize the perturbation */
void perturb_init(int n /* graph size at this level */
)
{
  extern int    NPERTURB;    /* number of edges to perturb */
  extern double PERTURB_MAX; /* maximum perturbation */
  int           i, j;        /* loop counter */
  double        drandom(void);

  /* Initialize the diagonal perturbation weights */
  pedges = smalloc(NPERTURB * sizeof(struct ipairs));
  pvals  = smalloc(NPERTURB * sizeof(double));

  if (n <= 1) {
    for (i = 0; i < NPERTURB; i++) {
      pedges[i].val1 = pedges[i].val2 = 0;
      pvals[i]                        = 0;
    }
    return;
  }

  for (i = 0; i < NPERTURB; i++) {
    pedges[i].val1 = 1 + (n * drandom());

    /* Find another vertex to define an edge. */
    j = 1 + (n * drandom());
    while (j == i) {
      j = 1 + (n * drandom());
    }
    pedges[i].val2 = 1 + (n * drandom());

    pvals[i] = PERTURB_MAX * drandom();
  }
}

void perturb_clear(void)
{

  sfree(pedges);
  sfree(pvals);
  pedges = NULL;
  pvals  = NULL;
}

/* Modify the result of splarax to break any graph symmetry */
void perturb(double *result, /* result of matrix-vector multiply */
             double *vec     /* vector matrix multiplies */
)
{
  extern int NPERTURB; /* number of edges to perturb */
  int        i;        /* loop counter */

  for (i = 0; i < NPERTURB; i++) {
    result[pedges[i].val1] += pvals[i] * (vec[pedges[i].val2] - vec[pedges[i].val1]);
    result[pedges[i].val2] += pvals[i] * (vec[pedges[i].val1] - vec[pedges[i].val2]);
  }
}

/* Modify the result of splarax to break any graph symmetry, float version */
void perturb_float(float *result, /* result of matrix-vector multiply */
                   float *vec     /* vector matrix multiplies */
)
{
  extern int NPERTURB; /* number of edges to perturb */
  int        i;        /* loop counter */

  for (i = 0; i < NPERTURB; i++) {
    result[pedges[i].val1] += pvals[i] * (vec[pedges[i].val2] - vec[pedges[i].val1]);
    result[pedges[i].val2] += pvals[i] * (vec[pedges[i].val1] - vec[pedges[i].val2]);
  }
}
