// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TIMER_H
#define __TIMER_H

#include "zoltan_timer.h"
#include <time.h> /* ANSI C; defines clock_t and clock() */

#ifndef CLOCKS_PER_SEC /* Should have been defined in time.h */
#define CLOCKS_PER_SEC 1000000 /* To prevent compile errors, not always the correct value. */
#endif

/*
 * POSIX compliant systems should use times() for user timing. 
 * This is the default in Zoltan. Make Zoltan with -DNO_TIMES if
 * your system does not have sys/times.h and times().
 * Note: BSD-like systems may use getrusage() instead for user timing,
 * but that has not been implemented here. 
 */

#if defined(__PUMAGON__) || defined(__LIBCATAMOUNT__) || defined(_WIN32)
/* Tflops with Cougar & Red Storm w/Catamount does not have sysconf() or times() */
/* Microsoft Visual Studio does not have times either */
#define NO_TIMES
#endif /* __PUMAGON__ */

#ifndef NO_TIMES
/* #include <sys/types.h> -- Included by sys/times.h on most systems. */
#include <sys/times.h>
#include <unistd.h> /* Needed for sysconf() and _SC_CLK_TCK */
#endif

#endif
