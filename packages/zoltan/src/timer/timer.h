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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <time.h> /* ANSI C; defines clock_t and clock() */
#include "timer_const.h"

#ifdef SMOS
/* Tflops is special. Only use MPI's wall time. */
#define NO_TIMES
#else /* !SMOS */
#include <unistd.h> /* Needed for sysconf() */
#ifndef CLOCKS_PER_SEC /* Should have been defined in time.h */
#define CLOCKS_PER_SEC 1000000 /* To prevent compile errors, not always the correct value. */
#endif
#endif /* !SMOS */

/*
 * POSIX compliant systems should use times() for user timing. 
 * This is the default in Zoltan. Make Zoltan with -DNO_TIMES if
 * your system does not have sys/times.h and times().
 * Note: BSD-like systems may use getrusage() instead for user timing.
 */

#ifndef NO_TIMES
/* #include <sys/types.h> -- Not required on most systems. */
#include <sys/times.h>
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
