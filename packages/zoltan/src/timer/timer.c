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

/*
 * Symbols being used in this file:
 *
 * _CLOCK_T    : if defined, clock() is available
 * HAVE_TIMES  : if defined, times() is available
 * HAVE_RUSAGE : if defined, getrusage() is available
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
#if (defined(_CLOCK_T) && ! defined(SMOS))
          (*timer) = ZOLTAN_TIME_CPU;
#else  /* SMOS or !_CLOCK_T */
          ZOLTAN_PRINT_WARN(-1, yo, "CPU time not available;"
                        " Wall clock time will be used.");
#endif /* SMOS or !_CLOCK_T */
        }
        else if (strcmp(result.sval, "USER")==0){
#if (defined(HAVE_TIMES) || defined(HAVE_RUSAGE))
          (*timer) = ZOLTAN_TIME_USER;
#else
          ZOLTAN_PRINT_WARN(-1, yo, "User time not available;"
                          " Wall clock time will be used.");
#endif
        }
        else if (strcmp(result.sval, "USERSYS")==0){
#if (defined(HAVE_TIMES) || defined(HAVE_RUSAGE))
          (*timer) = ZOLTAN_TIME_USERSYS;
#else
          ZOLTAN_PRINT_WARN(-1, yo, "Usersys time not "
                        "available; Wall clock time will be used.");
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
   The ANSI function clock() might roll over at approx 36 minutes,
   so we try to determine the number of rollovers.
*/

double Zoltan_Time(int timer)
{
  double t = -1.;

#if (defined(_CLOCK_T) && ! defined(SMOS))
  clock_t num_ticks;
  static clock_t last_num_ticks = 0;
  static double  inv_clocks_per_sec = 1./(double)CLOCKS_PER_SEC;
  static double  clock_width =
    (double)(1L<<((int)sizeof(clock_t)*8-2))*4./(double)CLOCKS_PER_SEC;
  static int     clock_rollovers = 0;
#endif /* !SMOS */
#ifdef HAVE_TIMES
  static double  inv_clk_tck = 1./(double)CLK_TCK;
#endif

  if (timer==ZOLTAN_TIME_WALL)
    /* Wall clock */
    t = MPI_Wtime();
  else if (timer==ZOLTAN_TIME_CPU) {
    /* CPU time */
#if (defined(_CLOCK_T) && ! defined(SMOS))
    num_ticks = clock();
    if (num_ticks < last_num_ticks) clock_rollovers++;
    t = num_ticks * inv_clocks_per_sec;
    if (clock_rollovers) t += clock_rollovers * clock_width;
    last_num_ticks = num_ticks;
#endif /* !SMOS */
  }
#if defined(HAVE_TIMES)
  else if (timer==ZOLTAN_TIME_USER) {
    struct tms tms;
    times(&tms);
    t = tms.tms_utime * inv_clk_tck;
  }
  else if (timer==ZOLTAN_TIME_USERSYS) {
    struct tms tms;
    times(&tms);
    t = (tms.tms_utime + tms.tms_stime) * inv_clk_tck;
  }
#elif defined(HAVE_RUSAGE)
  else if (timer==ZOLTAN_TIME_USER) {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    t = (rusage.ru_utime.tv_sec + 1.0e-6 * rusage.ru_utime.tv_usec);
  }
  else if (timer==ZOLTAN_TIME_USERSYS) {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    t = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
        1.0e-6 * (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));
  }
#endif /* HAVE_TIMES or HAVE_RUSAGE */

  return t;
}

/* Resolution (precision) of timer. If the precision is unknown, -1 is returned.  */

double Zoltan_Time_Resolution(int timer)
{
  double t = -1.;

  if (timer==ZOLTAN_TIME_WALL)
    t = MPI_Wtick();
  else if (timer==ZOLTAN_TIME_CPU)
    t = 1.0/(double)CLOCKS_PER_SEC;
#if defined(HAVE_TIMES)
  else if ((timer==ZOLTAN_TIME_USER)||(timer==ZOLTAN_TIME_USERSYS)) {
    t = 1.0/(double)CLK_TCK;
  }
#elif defined(HAVE_RUSAGE)
  else if ((timer==ZOLTAN_TIME_USER)||(timer==ZOLTAN_TIME_USERSYS)) {
#ifdef SUN
    /* Use Sun-specific variable */
    extern int usec_per_tick; 
    t = 1.0e-6 * usec_per_tick;
#endif /* !SUN */
  }
#endif /* HAVE_TIMES or HAVE_RUSAGE */

  return t;
}


