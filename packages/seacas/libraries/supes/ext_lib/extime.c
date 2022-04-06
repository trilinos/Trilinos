/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
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

#include <time.h>

#if defined(ADDC_)
void extime_(char *string, long int len)
#else
void extime(char *string, long int len)
#endif
{
  time_t     tim = time(0);
  struct tm *t   = localtime(&tim);
  strftime(string, 9, "%H:%M:%S", t);
}
