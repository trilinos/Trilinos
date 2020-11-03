/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "chaco_random.h"

/* Randomly permute elements of an array. */
void randomize(int *array, int n)
/* array of integer values */
/* number of values */
{
  int    i; /* loop counter */
  double drandom(void);

  for (i = 1; i <= n; i++) {
    double value = drandom();
    int    index = n * value + 1;
    int    temp  = array[i];
    array[i]     = array[index];
    array[index] = temp;
  }
}

double drandom(void) { return rand_rect_port(); }

void setrandom(long int seed) { init_rand_port(seed); }
