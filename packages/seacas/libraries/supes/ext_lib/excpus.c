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
 *   DESCRIPTION:
 *   This routine returns an accumulated CPU time. The base time is
 *   undefined; only relative times are valid.
 *
 *   FORMAL PARAMETERS:
 *   CPUSEC    REAL            CPU time
 *
 */

#include "fortranc.h"

#if defined(interix)
#include <time.h>
void excpus_(FTNREAL *cpusec)
#endif /* interix */

#if defined(aix) || defined(__VACPP__) || defined(hpux) || defined(sun) || defined(sgi) ||         \
    defined(__osf__) || defined(__linux__) || defined(__APPLE__) || defined(__CYGWIN__)

#include <sys/resource.h>
#include <sys/time.h>

#if defined(__NO_CYGWIN_OPTION__)
#include <windows.h>
#else
#include <sys/resource.h>
#include <sys/time.h>
#endif

#if defined ADDC_
    void excpus_(FTNREAL *cpusec)
#else
    void excpus(FTNREAL *cpusec)
#endif
#endif

{
#if defined(interix)
  *cpusec = ((FTNREAL)clock()) / CLOCKS_PER_SEC;
#endif /* interix */

#if defined(__NO_CYGWIN_OPTION__)
  FILETIME createTime;
  FILETIME exitTime;
  FILETIME kernelTimeProcess;
  FILETIME userTimeProcess;
  HANDLE   hProcess = GetCurrentProcess();
  LONGLONG tUser, tKernel;

  // Retrieve the file times for the process.
  if (hProcess) {
    // Get process times
    if (!GetProcessTimes(hProcess, &createTime, &exitTime, &kernelTimeProcess, &userTimeProcess))
      *cpusec = 0.0;
    else {
      // Convert to 64-bit LONGLONG
      tUser   = *(LONGLONG *)&userTimeProcess;
      tKernel = *(LONGLONG *)&kernelTimeProcess;
      // Time interval 100 nanoseconds, divide by 10,000,000
      // to convert to seconds
      *cpusec = ((double)(tUser) + (double)(tKernel)) / 10000000.0;
    }
  }
  else {
    *cpusec = 0.0;
  }
#endif

#if defined(sun) || defined(sgi) || defined(__osf__) || defined(__linux__) || defined(aix) ||      \
    defined(__VACPP__) || defined(paragon) || defined(hpux) || defined(__APPLE__)
  struct rusage rusage;
  int           secs, mics;

  getrusage(RUSAGE_SELF, &rusage);
  secs    = rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec;
  mics    = rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec;
  *cpusec = (FTNREAL)secs + ((FTNREAL)mics / 1.e6);
#endif
}
