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
#ifndef lint
static char *cvs_timerc_id = "$Id$";
#endif

#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"

/* Static variables */
static int Timer = LB_TIME_WALL;

/* Timer parameters */

int LB_Set_Timer_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    PARAM_VARS timer_params[] = {
        { "TIMER", &Timer, "STRING" },
        { NULL, NULL, NULL }
    };

    status = LB_Check_Param(name, val, timer_params, &result, &index);
/*
    printf("Debug: In LB_Set_Timer_Param, name=%s, val=%s, status=%d, index=%d\n", 
      name, val, status, index);
*/

    if (status == 0 && index == 0) {
        status = 3; /* Don't add to params list */
        if (strcasecmp(result.sval, "wall")==0)
          Timer = LB_TIME_WALL;
        else if (strcasecmp(result.sval, "cpu")==0) {
#if (defined(_CLOCK_T) && ! defined(SMOS))
          Timer = LB_TIME_CPU;
#else  /* SMOS or !_CLOCK_T */
          fprintf(stderr, "LB_Set_Timer_Param warning: CPU time not available;"
                          " Wall clock time will be used.\n");
          Timer = LB_TIME_WALL;
#endif /* SMOS or !_CLOCK_T */
        }
        else if (strcasecmp(result.sval, "user")==0){
#if (defined(HAVE_TIMES) || defined(HAVE_RUSAGE))
          Timer = LB_TIME_USER;
#else
          fprintf(stderr, "LB_Set_Timer_Param warning: user time not available;"
                          " Wall clock time will be used.\n");
          Timer = LB_TIME_WALL;
#endif
        }
        else if (strcasecmp(result.sval, "usersys")==0){
#if (defined(HAVE_TIMES) || defined(HAVE_RUSAGE))
          Timer = LB_TIME_USERSYS;
#else
          fprintf(stderr, "LB_Set_Timer_Param warning: usersys time not "
                          "available; Wall clock time will be used.\n");
          Timer = LB_TIME_WALL;
#endif
        }
        else{
          fprintf(stderr, "LB_Set_Timer_Param warning: Unknown timer option "
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



/* Print out a time with a message. Print max, sum, and imbalance. */

void LB_Print_Time (LB *lb, double time, char *msg)
{
  double sum, max;

  MPI_Reduce((void *)&time, (void *)&sum, 1, MPI_DOUBLE, MPI_SUM, 0, lb->Communicator);

  MPI_Reduce((void *)&time, (void *)&max, 1, MPI_DOUBLE, MPI_MAX, 0, lb->Communicator);

  if (lb->Proc == 0 && sum != 0.0)
    printf("%s: Max: %g, Sum: %g, Imbal.: %g\n", 
            msg, max, sum, max*(lb->Num_Proc)/sum-1.);

}
