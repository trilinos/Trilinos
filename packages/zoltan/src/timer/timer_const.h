/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#ifndef TIMER_CONST_H
#define TIMER_CONST_H

#ifndef lint
static char *cvs_timerconsth_id = "$Id$";
#endif

#include <time.h> /* ANSI C; defines clock_t and clock() */

/* POSIX.1 compliant systems should use times() for user+system timing */
#if (defined(_POSIX_) || defined(POSIX) || defined(POSIX1))
#define HAVE_TIMES
#endif

/* BSD-like systems should use getrusage() for user+system timing */
#if (defined(BSD) || defined(BSD_COMP) || defined(sun))
#define HAVE_RUSAGE
#endif

/* Include more header files depending on HAVE_* */
#if defined(HAVE_TIMES)
#include <sys/types.h>
#include <sys/times.h>
extern clock_t times();
#elif defined(HAVE_RUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
extern int getrusage();
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

/* Constants used in LB timer routines */
#define LB_TIME_WALL 1
#define LB_TIME_CPU  2
#define LB_TIME_USER 3
#define LB_TIME_USERSYS 4

/* Function prototypes */
double LB_Time();
double LB_Time_Resolution();
void LB_Print_Time (LB *, double, char *);
int LB_Set_Timer_Param(char *, char *);

#endif
