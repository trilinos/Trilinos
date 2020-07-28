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
C     DESCRIPTION:
C     This routine returns the current date in a character string. The
C     format is as follows:

C       YYYYMMDD

C     MM is a two digit month
C     DD is a two digit day
C     YYYY is a four digit year

C     This is known as the "Compact ISO 8601 format"

C     FORMAL PARAMETERS:
C     STRING    CHARACTER       String to receive the date

************************************************************************

*/

#define STRLEN 8
#include <stdio.h> // for sprintf
#include <time.h>  // for tm, localtime, time, time_t

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
void exdate_(char *string, long int len)
#else
void exdate(char *string, long int len)
#endif
{
  struct tm *t;
  time_t     tim;

  char Temp[STRLEN + 1]; /* Temporary string storage slot. */

  tim = time((time_t *)0);
  t   = localtime(&tim);
  t->tm_year += 1900;

  sprintf(Temp, "%04d%02d%02d", t->tm_year, t->tm_mon + 1, t->tm_mday);
  copy_string(string, Temp, STRLEN + 1);
}
