/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
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
#include <time.h>

#if defined(ADDC_)
void exdate_(char *string, long int len)
#else
void exdate(char *string, long int len)
#endif
{
  time_t     tim = time((time_t *)0);
  struct tm *t   = localtime(&tim);
  strftime(string, STRLEN + 1, "%Y%m%d", t);
}
