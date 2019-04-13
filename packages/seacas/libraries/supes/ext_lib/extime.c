/*
 * Copyright(C) 2008-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *
 * * Neither the name of NTESS nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*
 * $Id: extime.c,v 1.16 2008/03/14 13:22:38 gdsjaar Exp $
 */

/*
************************************************************************
C
C     DESCRIPTION:
C     This routine returns the current time in a CHARACTER string. The
C     format is as follows:
C
C       HH:MM:SS
C
C     HH is a two digit hour
C     MM is a two digit minute
C     SS is a two digit second
C
C     FORMAL PARAMETERS:
C     STRING    CHARACTER       String to receive the time
C
************************************************************************
C
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
  copy_string(string, Temp, STRLEN);
}
