/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"   // for FALSE, TRUE
#include "params.h" // for LINE_LENGTH
#include <ctype.h>
#include <stdio.h> // for getchar, sscanf

/* Robust routine to read an integer */
int input_int(void)
{
  char line[LINE_LENGTH]; /* space to read input line */
  int  done;              /* flag for end of integer */
  int  val;               /* value returned */
  int  i;                 /* loop counter */

  for (i = 0; i < LINE_LENGTH; i++) {
    line[i] = '\0';
  }

  i    = 0;
  done = FALSE;
  while (!done && i < LINE_LENGTH) {
    int c = getchar();
    if (c >= 0 && c <= 127) {
      line[i] = (char)c;
      if (isdigit(line[i]) || line[i] == '-') {
        i++;
      }
      else if (i != 0) {
        done = TRUE;
      }
    }
    else {
      done = TRUE;
    }
  }

  (void)sscanf(line, "%d", &val);
  return (val);
}
