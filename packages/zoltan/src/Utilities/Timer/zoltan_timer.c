/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#include "zoltan_timer.h"
#include "zoltan_types.h"
#include "zoltan_util.h"
#include "zoltan_mem.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/****************************************************************************/
/*
 * Functions that implement a Timer "class," creating a Timer object
 * with start, stop and print functions.
 *
 * This code was designed to be a stand-alone utility, relying only
 * on Zoltan_Time, Zoltan error codes, Zoltan utilities, and libzoltan_mem.a.
 * Some rearranging of the Zoltan_Time files is necessary to truly make
 * this utility a standalone one.
 */

/****************************************************************************/
/* Number of timers initially in Timer object. */
#define INITLENGTH 30   

/* Length of character strings naming each timer. */
#define MAXNAMELEN 30   

/* Flag indicating whether a timer is in use. */
#define INUSE 1

/* Flag indicating whether a timer is running. */
#define RUNNING 2


/* Macro to ensure that a given timer index is valid. */
#define TESTINDEX(zt, ts_idx, proc, yo) \
  if ((ts_idx) >= (zt)->NextTimeStruct) { \
    ZOLTAN_PRINT_ERROR(proc, yo, "Invalid Timer Index"); \
    return ZOLTAN_FATAL; \
  }

/****************************************************************************/
/* Structure that implements an individual timer. */
typedef struct TimeStruct {
  double Start_Time;      /* Most recent start time; 
                             set by Zoltan_Timer_Start */
  double Stop_Time;       /* Most recent end time;
                             set by Zoltan_Timer_Stop */
  double My_Tot_Time;     /* Sum of stop_time-start_time over all invocations
                             of this timer */
  double Global_Tot_Time; /* Sum of max stop_time-start_time over all procs
                             (in communicator) over all invocations of this 
                             timer */
  int Use_Barrier;        /* Flag indicating whether to perform a barrier
                             operation before starting the timer. */
  int Status;             /* Flag indicating status of TimeStruct:
                                > 0  -->  In Use
                                > 2  -->  Running */
  char Name[MAXNAMELEN];  /* String associated (and printed) with timer info */
} ZTIMER_TS;

/* Timer object consisting of many related timers. 
 * Applications access this structure. */
typedef struct Zoltan_Timer {
  int Timer_Flag;         /* Zoltan Timer_Flag flag passed to Zoltan_Time */
  int Length;             /* # of entries allocated in Times */
  int NextTimeStruct;     /* Index of next unused TimeStruct */
  ZTIMER_TS *Times;       /* Array of actual timing data -- individual timers */
} ZTIMER;

/****************************************************************************/
ZTIMER *Zoltan_Timer_Create(
  int timer_flag
)
{
/* Allocates a Timer object for the application; returns a pointer to it. 
 * Does not start any timers.
 */

ZTIMER *zt;
int i;

  zt = (ZTIMER *) ZOLTAN_MALLOC(sizeof(ZTIMER));
  zt->Times = (ZTIMER_TS *) ZOLTAN_MALLOC(sizeof(ZTIMER_TS) * INITLENGTH);
  zt->Timer_Flag = timer_flag;
  zt->Length = INITLENGTH;
  zt->NextTimeStruct = 0;

  for (i = 0; i < zt->Length; i++) 
    zt->Times[i].Status = 0;

  return zt;
}

/****************************************************************************/
int Zoltan_Timer_Init(
  ZTIMER *zt,           /* Ptr to Timer object */
  int proc,             /* Processor calling the Init */
  int use_barrier,      /* Flag indicating whether to perform a 
                           barrier operation before starting the
                           timer. */
  char *name            /* Name of this timer */
)
{
/* Function that returns the index of the next available Timer timer. */
int ret = zt->NextTimeStruct++;

  if (ret >= zt->Length) {
    /* Realloc -- need more individual timers */
    zt->Length += INITLENGTH;
    zt->Times = (ZTIMER_TS *) ZOLTAN_REALLOC(zt->Times, 
                                             zt->Length * sizeof(ZTIMER_TS));
  }

  Zoltan_Timer_Reset(zt, ret, proc, use_barrier, name);

  return ret;
}


/****************************************************************************/
int Zoltan_Timer_Reset(
  ZTIMER *zt,
  int ts_idx,            /* Index of the timer to reset */
  int proc,              /* Proc calling reset */
  int use_barrier,       /* Flag indicating whether to perform a 
                            barrier operation before starting the
                            timer. */
  char *name             /* Name of this timer */
)
{
/* Initialize a timer for INUSE; reset its values to zero. */
static char *yo = "Zoltan_Timer_Reset";
ZTIMER_TS *ts;

  TESTINDEX(zt, ts_idx, proc, yo);

  ts = &(zt->Times[ts_idx]);

  ts->Status = INUSE;
  ts->Start_Time = 0.;
  ts->Stop_Time = 0.;
  ts->My_Tot_Time = 0.;
  ts->Global_Tot_Time = 0.;
  ts->Use_Barrier = use_barrier;
  strcpy(ts->Name, name);

  return ZOLTAN_OK;
}

/****************************************************************************/
int Zoltan_Timer_Start(
  ZTIMER *zt,            /* Ptr to Timer object */
  int ts_idx,            /* Index of the timer to use */
  int proc,              /* Processor calling Start */
  MPI_Comm comm          /* Communicator to use for synchronization, 
                            if requested */
)
{
ZTIMER_TS *ts;
static char *yo = "Zoltan_Timer_Start";

  TESTINDEX(zt, ts_idx, proc, yo);

  ts = &(zt->Times[ts_idx]);
  if (ts->Status > 2) {
    ZOLTAN_PRINT_ERROR(proc, yo,
                       "Cannot start timer; timer is already running.");
    return ZOLTAN_FATAL;
  }

  ts->Status += RUNNING;
  if (ts->Use_Barrier)
    MPI_Barrier(comm);

  ts->Start_Time = Zoltan_Time(zt->Timer_Flag);
  return ZOLTAN_OK;
}

/****************************************************************************/

int Zoltan_Timer_Stop(
  ZTIMER *zt,            /* Ptr to Timer object */
  int ts_idx,            /* Index of the timer to use */
  int proc,              /* Processor calling Stop */
  MPI_Comm comm          /* Communicator to use for AllReduce */
)
{
/* Function to stop a timer and accrue its information */
ZTIMER_TS *ts;
static char *yo = "Zoltan_Timer_Stop";
double my_time, global_time;

  TESTINDEX(zt, ts_idx, proc, yo);

  ts = &(zt->Times[ts_idx]);
  if (ts->Status < 2) {
    ZOLTAN_PRINT_ERROR(proc, yo,
                       "Cannot stop timer; timer is not running.");
    return ZOLTAN_FATAL;
  }

  ts->Status -= RUNNING;
  ts->Stop_Time = Zoltan_Time(zt->Timer_Flag);
  my_time = ts->Stop_Time - ts->Start_Time;

  ts->My_Tot_Time += my_time;

  MPI_Allreduce(&my_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, comm);
  ts->Global_Tot_Time += global_time;

  return ZOLTAN_OK;
}


/****************************************************************************/
int Zoltan_Timer_Print(
  ZTIMER *zt,
  int ts_idx,
  int proc
)
{
/* Print a single timer's information */
static char *yo = "Zoltan_Timer_Print";
ZTIMER_TS *ts;

  TESTINDEX(zt, ts_idx, proc, yo);
  ts = &(zt->Times[ts_idx]);

  printf("%d ZOLTAN_TIMER %d %s:   ProcTime %e   GlobalTime %e\n", 
         proc, ts_idx, ts->Name, ts->My_Tot_Time, ts->Global_Tot_Time);
  return ZOLTAN_OK;
}

/****************************************************************************/
int Zoltan_Timer_PrintAll(
  ZTIMER *zt,
  int proc
)
{
/* Function to print all timer information */
int i, ierr = ZOLTAN_OK;
ZTIMER_TS *ts;

  for (i = 0; i < zt->NextTimeStruct; i++) 
    if ((ierr = Zoltan_Timer_Print(zt, i, proc)) != ZOLTAN_OK)
      break;

  return ierr;
}

/****************************************************************************/
void Zoltan_Timer_Destroy(
  ZTIMER **zt
)
{
/* Destroy a Timer object */
  ZOLTAN_FREE(&((*zt)->Times));
  ZOLTAN_FREE(zt);
}

/****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
