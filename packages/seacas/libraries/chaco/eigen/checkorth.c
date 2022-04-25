/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"  // for min
#include <math.h>  // for fabs, log10
#include <stdio.h> // for printf

/* Check orthogonality of vector set */
void checkorth(double **mat, int n, int dim)
{
  int    i, j;      /* loop indices */
  double measure;   /* Froebenius norm */
  double prod;      /* value of dot product */
  double worst;     /* greatest off-diagonal dot product */
  int    lim;       /* index of last vec to check against */
  int    screenlim; /* value of lim that will fit on screen */
  int    option;    /* which option to use */

  double dot(double *vec1, int beg, int end, double *vec2); /* standard dot product routine */

  /* The T/F argument in the conditionals is just a convenient option: */

  screenlim = 20;
  option    = 3;

  /* Check orthogonality over whole set. */
  if (option == 1) {
    printf("Orthogonality check:\n");
    for (i = 1; i <= dim; i++) {
      printf("%2d)", i);
      for (j = 1; j <= i; j++) {
        prod = dot(mat[i], 1, n, mat[j]);
        /* printf(" %g ",prod); */
        /* printf(" %4.2e ",prod); */
        /* printf(" %4.2e ",fabs(prod)); */
        printf(" %2d", -(int)log10(prod));
      }
      printf("\n");
    }
  }

  if (option == 2) {
    printf("Frobenius orthogonality measure:");
    measure = 0;
    for (i = 1; i <= dim; i++) {
      for (j = i; j <= dim; j++) {
        prod = dot(mat[i], 1, n, mat[j]);
        if (i == j) {
          measure += fabs(1.0 - prod);
        }
        else {
          measure += 2.0 * fabs(prod);
        }
      }
    }
    printf("%g \n", measure);
  }

  /* Check orthogonality against last vector. Allows you to build up orthogonality
     matrix much faster if previous columns stay the same when add a new column,
     but may interact with other debug output to give a confusing presentation. */
  if (option == 3) {
    printf("%3d) ", dim);
    lim   = min(dim, screenlim);
    worst = 0;
    for (i = 1; i <= dim; i++) {
      prod = dot(mat[i], 1, n, mat[dim]);
      if (i <= lim) {
        printf(" %2d", -(int)log10(fabs(prod)));
      }
      if ((i != dim) && (fabs(prod) > fabs(worst))) {
        worst = prod;
      }
    }
    printf(" worst %4.2e\n", worst);
  }
}

/* Check orthogonality of vector set */
void checkorth_float(float **mat, int n, int dim)
{
  int    i, j;      /* loop indices */
  double measure;   /* Froebenius norm */
  double prod;      /* value of dot product */
  double worst;     /* greatest off-diagonal dot product */
  int    lim;       /* index of last vec to check against */
  int    screenlim; /* value of lim that will fit on screen */
  int    option;    /* which option to use */

  double dot_float(float *vec1, int beg, int end, float *vec2); /* standard dot product routine */

  /* The T/F argument in the conditionals is just a convenient option: */

  screenlim = 20;
  option    = 3;

  /* Check orthogonality over whole set. */
  if (option == 1) {
    printf("Orthogonality check:\n");
    for (i = 1; i <= dim; i++) {
      printf("%2d)", i);
      for (j = 1; j <= i; j++) {
        prod = dot_float(mat[i], 1, n, mat[j]);
        /* printf(" %g ",prod); */
        /* printf(" %4.2e ",prod); */
        /* printf(" %4.2e ",fabs(prod)); */
        printf(" %2d", -(int)log10(prod));
      }
      printf("\n");
    }
  }

  if (option == 2) {
    printf("Frobenius orthogonality measure:");
    measure = 0;
    for (i = 1; i <= dim; i++) {
      for (j = i; j <= dim; j++) {
        prod = dot_float(mat[i], 1, n, mat[j]);
        if (i == j) {
          measure += fabs(1.0 - prod);
        }
        else {
          measure += 2.0 * fabs(prod);
        }
      }
    }
    printf("%g \n", measure);
  }

  /* Check orthogonality against last vector. Allows you to build up orthogonality
     matrix much faster if previous columns stay the same when add a new column,
     but may interact with other debug output to give a confusing presentation. */
  if (option == 3) {
    printf("%3d) ", dim);
    lim   = min(dim, screenlim);
    worst = 0;
    for (i = 1; i <= dim; i++) {
      prod = dot_float(mat[i], 1, n, mat[dim]);
      if (i <= lim) {
        printf(" %2d", -(int)log10(fabs(prod)));
      }
      if ((i != dim) && (fabs(prod) > fabs(worst))) {
        worst = prod;
      }
    }
    printf(" worst %4.2e\n", worst);
  }
}
