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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <mpi.h>
#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */

#include "timer.h"
#include "params_const.h"
#include "zoltan_util.h"

/*
 * Machine independent timing utilities.
 * ANSI C and MPI are required.
 */

/* Interpret and set timer parameters */

int Zoltan_Set_Timer_Param(
char *name,                     /* input:  name of variable */
char *val,                      /* input:  value of variable */
int *timer)                     /* output: timer type */
{
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    int status;
    PARAM_VARS Timer_params[] = {
        { "TIMER", NULL, "STRING" },
        { NULL, NULL, NULL }
    };
    char *yo = "Zoltan_Set_Timer_Param";

    (*timer) = ZOLTAN_TIME_WALL;  /* default timer value */

    status = Zoltan_Check_Param(name, val, Timer_params, &result, &index);

    if (status == 0 && index == 0) {
        if (!strcmp(result.sval, "WALL"))
          (*timer) = ZOLTAN_TIME_WALL;
        else if (strcmp(result.sval, "CPU")==0) {
#if (! defined(SMOS))
          (*timer) = ZOLTAN_TIME_CPU;
#else  /* SMOS */
          ZOLTAN_PRINT_WARN(-1, yo, "CPU time not available;"
                        " Wall clock time will be used.");
#endif /* SMOS */
        }
        else if (strcmp(result.sval, "USER")==0){
#if (! defined(NO_TIMES))
          (*timer) = ZOLTAN_TIME_USER;
#else
          ZOLTAN_PRINT_WARN(-1, yo, "User time not available;"
                          " Wall clock time will be used.");
#endif
        }
        else{
          char msg[256];
          sprintf(msg, "Unknown timer option %s.", result.sval);
          ZOLTAN_PRINT_WARN(-1, yo, msg);
          status = 2; /* Illegal parameter */
        }
    }

    return(status);
}


/* Timer routine that returns either CPU or wall-clock time.
   The ANSI function clock() might roll over at approx 71.5 minutes,
   so we try to determine the number of rollovers.
*/

double Zoltan_Time(int timer)
{
  double t = -1.;

#if (! defined(SMOS))
  clock_t num_ticks;
  static clock_t last_num_ticks = 0;
  static int     clock_rollovers = 0;
  static double  secs_per_clock = 1./((double) CLOCKS_PER_SEC);
  static double  clock_width = 
    ((double)(1L<<((int)sizeof(clock_t)*8-2)))*4./((double) CLOCKS_PER_SEC);
  static double  secs_per_tick  = 0.; /* Not necessarily the same as secs_per_clock; system-dependent; get value from sysconf(). */
#endif /* !SMOS */

  if (timer==ZOLTAN_TIME_WALL)
    /* Wall clock */
    t = MPI_Wtime();
  else if (timer==ZOLTAN_TIME_CPU) {
    /* CPU time */
#if (! defined(SMOS))
    num_ticks = clock();
    if (num_ticks < last_num_ticks) clock_rollovers++;
    t = num_ticks * secs_per_clock;
    if (clock_rollovers) t += clock_rollovers * clock_width;
    last_num_ticks = num_ticks;
#endif /* !SMOS */
  }
#if (! defined(NO_TIMES))
  else if (timer==ZOLTAN_TIME_USER) {
    struct tms tm;
    if (secs_per_tick == 0.)
      secs_per_tick = 1. / ((double) sysconf(_SC_CLK_TCK));
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
#if (! defined(SMOS))
  else if (timer==ZOLTAN_TIME_CPU)
    t = 1. / ((double) CLOCKS_PER_SEC);
#endif /* !SMOS */
#if (! defined(NO_TIMES))
  else if (timer==ZOLTAN_TIME_USER){
    t = 1. / ((double) sysconf(_SC_CLK_TCK));
  }
#endif 

  return t;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
