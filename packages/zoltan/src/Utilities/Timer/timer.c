// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "timer.h"
#include "zoltan_util.h"
#ifdef __PUMAGON__
#include <nx.h>
#endif

/****************************************
 * Machine independent timing utilities.
 * ANSI C and MPI are required.
 ****************************************/


/* Timer routine that returns either CPU or wall-clock time.
   The ANSI function clock() may roll over at approx 71.5 minutes,
   on some systems so we try to determine the number of rollovers.
*/

double Zoltan_Time(int timer)
{
  double t = -1.;

#if !defined(__PUMAGON__) && !defined(__LIBCATAMOUNT__)
  clock_t num_ticks;
  static clock_t last_num_ticks = 0;
  static int     clock_rollovers = 0;
  static double  secs_per_clock = (double) 1./((double) CLOCKS_PER_SEC);
  static double  clock_width = 
    ((double)(1L<<((int)sizeof(clock_t)*8-2)))*4./((double) CLOCKS_PER_SEC);
  static double  secs_per_tick  = 0.; /* Not necessarily the same as 
         secs_per_clock; system-dependent; get value from sysconf(). */
#endif

  if (timer==ZOLTAN_TIME_WALL)
    /* Wall clock */
    t = MPI_Wtime();
  else if (timer==ZOLTAN_TIME_CPU) {
#if defined(__PUMAGON__) || defined(__LIBCATAMOUNT__)
    /* CPU time on ASCI Red and Red Storm. */
    t = dclock();
#else
    /* CPU time */
    num_ticks = clock();
    if (num_ticks < last_num_ticks) clock_rollovers++;
    t = num_ticks * secs_per_clock;
    if (clock_rollovers) t += clock_rollovers * clock_width;
    last_num_ticks = num_ticks;
#endif
  }
#ifndef NO_TIMES
  else if (timer==ZOLTAN_TIME_USER) {
    struct tms tm;
    if (secs_per_tick == 0.)
      secs_per_tick = (double) 1. / ((double) sysconf(_SC_CLK_TCK));
    times(&tm);
    t = tm.tms_utime * secs_per_tick;
  }
#endif

  return t;
}


/* Resolution (precision) of timer. 
 * This is really a lower bound, the actual resolution may be worse.
 * If the precision is unknown, -1 is returned.  
 */

double Zoltan_Time_Resolution(int timer)
{
  double t = -1.;

  if (timer==ZOLTAN_TIME_WALL)
    t = MPI_Wtick();
  else if (timer==ZOLTAN_TIME_CPU)
    t = (double) 1. / ((double) CLOCKS_PER_SEC);
#ifndef NO_TIMES
  else if (timer==ZOLTAN_TIME_USER)
    t = (double) 1. / ((double) sysconf(_SC_CLK_TCK));
#endif 

  return t;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
