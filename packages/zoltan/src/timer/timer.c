/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * Zoltan is distributed under the GNU Lesser General Public License 2.1.    * 
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"

/*
 * Symbols being used in this file:
 *
 * _CLOCK_T    : if defined, clock() is available
 * HAVE_TIMES  : if defined, times() is available
 * HAVE_RUSAGE : if defined, getrusage() is available
 */
 
/* Static variables */
static int Timer = LB_TIME_WALL;
static PARAM_VARS Timer_params[] = {
        { "TIMER", &Timer, "STRING" },
        { NULL, NULL, NULL }
};

/* Interpret and set timer parameters */

int LB_Set_Timer_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    int status;

    status = LB_Check_Param(name, val, Timer_params, &result, &index);

    if (status == 0 && index == 0) {
        if (result.def || (!strcmp(result.sval, "WALL")))
          Timer = LB_TIME_WALL;
        else if (strcmp(result.sval, "CPU")==0) {
#if (defined(_CLOCK_T) && ! defined(SMOS))
          Timer = LB_TIME_CPU;
#else  /* SMOS or !_CLOCK_T */
          fprintf(stderr, "LB_Set_Timer_Param warning: CPU time not available;"
                          " Wall clock time will be used.\n");
          Timer = LB_TIME_WALL;
#endif /* SMOS or !_CLOCK_T */
        }
        else if (strcmp(result.sval, "USER")==0){
#if (defined(HAVE_TIMES) || defined(HAVE_RUSAGE))
          Timer = LB_TIME_USER;
#else
          fprintf(stderr, "LB_Set_Timer_Param warning: user time not available;"
                          " Wall clock time will be used.\n");
          Timer = LB_TIME_WALL;
#endif
        }
        else if (strcmp(result.sval, "USERSYS")==0){
#if (defined(HAVE_TIMES) || defined(HAVE_RUSAGE))
          Timer = LB_TIME_USERSYS;
#else
          fprintf(stderr, "LB_Set_Timer_Param warning: usersys time not "
                          "available; Wall clock time will be used.\n");
          Timer = LB_TIME_WALL;
#endif
        }
        else{
          fprintf(stderr, "Zoltan warning: Unknown timer option "
                          "%s\n", result.sval);
          status = 2; /* Illegal parameter */
        }
    }

    return(status);
}

/* Timer routine that returns either CPU or wall-clock time.
   The ANSI function clock() might roll over at approx 36 minutes,
   so we try to determine the number of rollovers.
*/

double LB_Time()
{
  double t = 0.0;
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

  if (Timer==LB_TIME_WALL)
    /* Wall clock */
    t = MPI_Wtime();
  else if (Timer==LB_TIME_CPU) {
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
  else if (Timer==LB_TIME_USER) {
    struct tms tms;
    times(&tms);
    t = tms.tms_utime * inv_clk_tck;
  }
  else if (Timer==LB_TIME_USERSYS) {
    struct tms tms;
    times(&tms);
    t = (tms.tms_utime + tms.tms_stime) * inv_clk_tck;
  }
#elif defined(HAVE_RUSAGE)
  else if (Timer==LB_TIME_USER) {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    t = (rusage.ru_utime.tv_sec + 1.0e-6 * rusage.ru_utime.tv_usec);
  }
  else if (Timer==LB_TIME_USERSYS) {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    t = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
        1.0e-6 * (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));
  }
#endif /* HAVE_TIMES or HAVE_RUSAGE */

  return t;
}

/* Precision of timer. If the precision is unknown, zero is returned.  */
double LB_Time_Resolution()
{
  double t = 0.0;

  if (Timer==LB_TIME_WALL)
    t = MPI_Wtick();
  else if (Timer==LB_TIME_CPU)
    t = 1.0/(double)CLOCKS_PER_SEC;
#if defined(HAVE_TIMES)
  else if ((Timer==LB_TIME_USER)||(Timer==LB_TIME_USERSYS)) {
    t = 1.0/(double)CLK_TCK;
  }
#elif defined(HAVE_RUSAGE)
  else if ((Timer==LB_TIME_USER)||(Timer==LB_TIME_USERSYS)) {
#ifdef SUN
    /* Use Sun-specific variable */
    extern int usec_per_tick; 
    t = 1.0e-6 * usec_per_tick;
#endif /* !SUN */
  }
#endif /* HAVE_TIMES or HAVE_RUSAGE */

  return t;
}


