/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __TIMER_H
#define __TIMER_H

#include <time.h> /* ANSI C; defines clock_t and clock() */
#include "timer_const.h"

/* Skip advanced timers. They give us more trouble than
   they are useful.  */
#ifdef USE_ADVANCED_TIMERS 

/* POSIX compliant systems should use times() for user+system timing */
#if (defined(_POSIX_) || defined(POSIX) || \
     defined(_POSIX_SOURCE) || defined(_POSIX_C_SOURCE))
#define HAVE_TIMES
#endif

/* BSD-like systems should use getrusage() for user+system timing */
#if (defined(BSD) || defined(BSD43) || defined(sun))
#define HAVE_RUSAGE
#endif

#endif /* USE_ADVANCED_TIMERS */

/* Include more header files depending on HAVE_* */
#if defined(HAVE_TIMES)
#include <sys/types.h>
#include <sys/times.h>
/* extern clock_t times(); */ /* Should be defined in sys/times.h */
#elif defined(HAVE_RUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
extern int getrusage(); /* Should be in sys/resource.h, but isn't always */
#endif

/* Some non-ANSI systems return clock() in microseconds
 * and do not define CLOCKS_PER_SEC.
 */
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000 /* Correct for SunOs 4.1 */
#endif
/* Similarly, guess the value of CLK_TCK if it is not defined */
#ifndef CLK_TCK
#define CLK_TCK 60 /* Correct for SunOs 4.1 */
#endif

#endif
