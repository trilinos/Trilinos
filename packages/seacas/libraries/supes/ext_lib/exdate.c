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
 * $Id: exdate.c,v 1.19 2008/03/14 13:22:36 gdsjaar Exp $
 */

/*
C     DESCRIPTION:
C     This routine returns the current date in a character string. The
C     format is as follows:
C
C       YYYYMMDD
C
C     MM is a two digit month
C     DD is a two digit day
C     YYYY is a four digit year
C
C     This is known as the "Compact ISO 8601 format"
C
C     FORMAL PARAMETERS:
C     STRING    CHARACTER       String to receive the date
C
************************************************************************
C
*/

#define STRLEN 8
#include <stdio.h>  // for sprintf
#include <string.h> // for strncpy
#include <time.h>   // for tm, localtime, time, time_t

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
  strncpy(string, Temp, STRLEN);
}
