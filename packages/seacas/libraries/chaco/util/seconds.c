/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef _MSC_VER
#include <sys/time.h>
#else
#include <time.h>
#endif
#if !defined(__CYGWIN__) && !defined(_MSC_VER)
#include <sys/resource.h>
#endif

double seconds(void)
{
  double curtime;

#ifdef RUSAGE_SELF

  /* This timer is faster and more robust (if it exists). */
  struct rusage rusage;
  int           getrusage(int, struct rusage *);

  getrusage(RUSAGE_SELF, &rusage);
  curtime = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
             1.0e-6 * (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));
#else
  /* ANSI timer, but lower resolution & wraps around after ~36 minutes. */
  curtime = clock() / ((double)CLOCKS_PER_SEC);
#endif

  return (curtime);
}
