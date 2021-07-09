/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>

/* Record a return TRUE if answer is yes, FALSE if no. */
int affirm(char *prompt)
{
  int  reply; /* character typed in */
  int  done;  /* loop control */
  void bail(char *msg, int status);

  if (prompt != NULL && (int)strlen(prompt) > 0) {
    printf("%s? ", prompt);
  }
  done = 0;
  while (!done) {
    reply = getchar();
    /* while (reply == ' ' || reply== '\n') reply= getchar(); */
    while (isspace(reply)) {
      reply = getchar();
    }

    if (reply == 'y' || reply == 'Y') {
      done = 1;
    }
    else if (reply == 'n' || reply == 'N') {
      done = 2;
    }
    else if (reply == 'q' || reply == 'Q') {
      done = 3;
    }
    else if (reply == 'x' || reply == 'X') {
      done = 3;
    }
    else {
      printf("Valid responses begin with: y Y n N q Q x X\n");
      if (prompt != NULL) {
        printf("%s? ", prompt);
      }
      /* Flush rest of input line. */
      while (reply != '\n') {
        reply = getchar();
      }
    }
  }
  if (done > 2) {
    bail(NULL, 0);
  }
  else if (done == 2) {
    return (FALSE);
  }
  return (TRUE);
}
