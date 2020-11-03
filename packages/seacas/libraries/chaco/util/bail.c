/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h>
#include <string.h>

/* Wrapper for exit() - print message and exit with status code. Exit code
   of 0 indicates normal termination. Exit code of 1 indicates early
   termination following detection of some problem. Call with bail(NULL,status)
   to suppress message. */
void bail(char *msg, int status)
{
  extern FILE *Output_File; /* Output file or NULL */
  void         exit(int);

  if (msg != NULL && (int)strlen(msg) > 0) {
    printf("%s\n", msg);
    if (Output_File != NULL) {
      fprintf(Output_File, "%s\n", msg);
    }
  }
  exit(status);
}
