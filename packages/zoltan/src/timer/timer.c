#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"

/* Static variables */
static int TIMER = TIME_WALL;

/* Timer parameters */

int LB_Set_Timer_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    PARAM_VARS timer_params[] = {
        { "TIMER", &TIMER, "INT" },
        { NULL, NULL, NULL }
    };

    status = LB_Check_Param(name, val, timer_params, &result, &index);
    if (status == 0 && index == 0) {
        status = 3; /* Don't add to params list */
        if (strcasecmp(result.sval, "wall")==0)
          TIMER = TIME_WALL;
        else if (strcasecmp(result.sval, "cpu")==0)
          TIMER = TIME_CPU;
        else
          status = 2; /* Unknown timer */
    }

    return(status);
}

/* Timer routine that returns either CPU or wall-clock time */

#ifndef CLOCKS_PER_SEC
/* Some non-ANSI systems return clock() in microseconds 
 * and do not define CLOCKS_PER_SEC.
 */
#define CLOCKS_PER_SEC 1000000
#endif

double LB_Time()
{
  double t;

  if (TIMER==TIME_WALL)
    /* Wall clock */
    t = MPI_Wtime();
  else if (TIMER==TIME_CPU)
    /* CPU time */
    t = ((double) clock())/CLOCKS_PER_SEC;
  else
    t = 0.0;  /* Error */

  return t;
}

/* Precision of timer */
double LB_Time_Resolution()
{
  double t;

  if (TIMER==TIME_WALL)
    t = MPI_Wtick();
  else if (TIMER==TIME_CPU)
    t = 1.0/CLOCKS_PER_SEC;
  else
    t = 0.0;  /* Error */

  return t;
}



/* Print out a time with a message. Print max, sum, and imbalance. */

void LB_Print_Time (LB *lb, double time, char *msg)
{
  double sum, max;

  MPI_Reduce((void *)&time, (void *)&sum, 1, MPI_DOUBLE, MPI_SUM, 0, lb->Communicator);

  MPI_Reduce((void *)&time, (void *)&max, 1, MPI_DOUBLE, MPI_MAX, 0, lb->Communicator);

  if (lb->Proc == 0 && sum != 0.0)
    printf("%s: Max: %7.3f, Sum: %7.3f, Imbal.: %7.3f\n", 
            msg, max, sum, max*(lb->Num_Proc)/sum-1.);

}
