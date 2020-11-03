/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h>
#include <string.h>

/* Debug break point. */
void checkpnt(char *tag)
{
  int  affirm(char *prompt);
  void bail(char *msg, int status);

  if (tag != NULL && (int)strlen(tag) > 0) {
    printf("%s: ", tag);
  }
  if (!affirm("continue")) {
    bail(NULL, 0);
  }
}
