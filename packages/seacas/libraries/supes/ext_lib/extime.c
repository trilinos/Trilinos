/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 */

/*
************************************************************************

C     DESCRIPTION:
C     This routine returns the current time in a CHARACTER string. The
C     format is as follows:

C       HH:MM:SS

C     HH is a two digit hour
C     MM is a two digit minute
C     SS is a two digit second

C     FORMAL PARAMETERS:
C     STRING    CHARACTER       String to receive the time

************************************************************************

*/

#define STRLEN 8
#include <stdio.h>
#include <string.h>
#include <time.h>

static char *copy_string(char *dest, char const *source, long int elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

#if defined(ADDC_)
void extime_(char *string, long int len)
#else
void extime(char *string, long int len)
#endif
{
  struct tm *t;
  time_t     tim;

  char Temp[STRLEN + 1]; /* My temporary string storage slot. */

  tim = time(0);
  t   = localtime(&tim);

  sprintf(Temp, "%02d:%02d:%02d", t->tm_hour, t->tm_min, t->tm_sec);
  copy_string(string, Temp, STRLEN + 1);
}
