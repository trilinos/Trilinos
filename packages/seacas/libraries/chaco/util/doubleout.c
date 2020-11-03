/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <math.h>
#include <stdio.h>

/* Print a double precision number with filtering format to screen. */
void doubleout(double number, int mode)
/* argument to print */
/* currently just one */
{
  if (mode == 1) {
    if (fabs(number) < 100) {
      printf("  %19.16f", number);
    }
    else {
      printf("  %19g", number);
    }
  }
}

/* Print a double precision number with filtering format to file. */
void doubleout_file(FILE *outfile, double number, int mode)
/* output file if not NULL */
/* argument to print */
/* currently just one */
{
  if (outfile == NULL) {
    return;
  }

  if (mode == 1) {
    if (fabs(number) < 100) {
      fprintf(outfile, "  %19.16f", number);
    }
    else {
      fprintf(outfile, "  %19g", number);
    }
  }
}
