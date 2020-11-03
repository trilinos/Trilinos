/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h> // for fprintf, printf, FILE, NULL

/* Wrapper for a printf statement with a string as only arg.
   Prints to screen and to output file if there is one. */
void strout(char *msg)
{
  extern FILE *Output_File; /* output file or null */

  printf("%s\n", msg);
  if (Output_File != NULL) {
    fprintf(Output_File, "%s\n", msg);
  }
}
