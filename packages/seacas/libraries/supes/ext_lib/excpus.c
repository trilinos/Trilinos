/*
 * Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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
 * * Neither the name of Sandia Corporation nor the names of its
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
 * $Id: excpus.c,v 1.26 2008/03/14 13:22:35 gdsjaar Exp $
 */

/*
 *   DESCRIPTION:
 *   This routine returns an accumulated CPU time. The base time is
 *   undefined; only relative times are valid.
 *
 *   FORMAL PARAMETERS:
 *   CPUSEC    REAL            CPU time
 *
 */

#include "fortranc.h"

#if defined (interix)
#include <time.h>
     void excpus_(FTNREAL *cpusec)
#undef linux
#endif             /* interix */

#if defined(aix) || defined(__VACPP__) || defined(hpux) || defined (sun) || defined (sgi) || defined (__osf__) || defined(__linux__) || defined (__APPLE__) || defined(__CYGWIN__) || defined(p6)

#if defined(__NO_CYGWIN_OPTION__)
#include <windows.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif
#if defined ADDC_       
     void excpus_(FTNREAL *cpusec)
#else
     void excpus(FTNREAL *cpusec)
#endif       
#endif

#if defined (paragon) || defined (pumagon) || defined(cougar)

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <nx.h>
     void excpus_(FTNREAL *cpusec)
#endif

{
#if defined(interix)
  *cpusec = ((FTNREAL) clock()) / CLOCKS_PER_SEC;
#endif              /* interix */

#if defined(__NO_CYGWIN_OPTION__)
  FILETIME createTime;
  FILETIME exitTime;
  FILETIME kernelTimeProcess;
  FILETIME userTimeProcess;
  HANDLE hProcess = GetCurrentProcess();
  LONGLONG tUser, tKernel;

  // Retrieve the file times for the process.
  if (hProcess)
  {
    // Get process times
    if (!GetProcessTimes(hProcess, &createTime, &exitTime,
                         &kernelTimeProcess, &userTimeProcess))
        *cpusec = 0.0;
    else
    {
      // Convert to 64-bit LONGLONG
      tUser = *(LONGLONG *)&userTimeProcess;
      tKernel = *(LONGLONG *)&kernelTimeProcess;
      // Time interval 100 nanoseconds, divide by 10,000,000 
      // to convert to seconds
      *cpusec = ((double)(tUser)+(double)(tKernel))/10000000.0;
    }
  }
  else
  {
    *cpusec = 0.0;
  }
#endif

#if defined (sun) || defined (sgi) || defined (__osf__) || defined(linux) || defined(aix) || defined(__VACPP__) || defined(paragon) || defined(hpux) || defined(__APPLE__)
  struct rusage rusage;
  int secs,mics;

  getrusage(RUSAGE_SELF,&rusage);
  secs = rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec;
  mics = rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec;
  *cpusec = (FTNREAL)secs + ((FTNREAL)mics / 1.e6);
#endif

#if defined (pumagon)
     *cpusec = dclock();
# endif     

#if defined (cougar)
     double tickval;
     tickval = dclock();
     *cpusec = tickval;
#  endif     

#  if defined (p6)
     *cpusec = 0.0;
#  endif     
}

