// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include "zoltan_timer.h"
#include "zoltan_types.h"
#include "zoltan_util.h"
#include "zoltan_mem.h"
#include "zoltan_comm.h"

#ifdef VAMPIR
#include <VT.h>
#endif

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
/* If you change this constant, change the string format  */
/* in Zoltan_Timer_Print, too. */
#define MAXNAMELEN 31   

/* Flag indicating whether a timer is in use. */
#define INUSE 1

/* Flag indicating whether a timer is running. */
#define RUNNING 2

#define FATALERROR(yo, str) \
  { \
    int ppproc; \
    MPI_Comm_rank(zoltan_get_global_comm(), &ppproc); \
    ZOLTAN_PRINT_ERROR(ppproc, yo, str); \
    return ZOLTAN_FATAL; \
  }
  

/* Macro to ensure that a Timer object is non-NULL */
#define TESTTIMER(zt, yo) \
  if ((zt) == NULL) FATALERROR(yo, "NULL Zoltan_Timer")

/* Macro to ensure that a given timer index is valid. */
#define TESTINDEX(zt, ts_idx, yo) \
  if ((ts_idx) >= (zt)->NextTimeStruct) FATALERROR(yo, "Invalid Timer Index")


/****************************************************************************/
/* Structure that implements an individual timer. */
typedef struct TimeStruct {
  double Start_Time;      /* Most recent start time; 
                             set by Zoltan_Timer_Start */
  double Stop_Time;       /* Most recent end time;
                             set by Zoltan_Timer_Stop */
  char Start_File[MAXNAMELEN+1];  /* Filename for most recent Start */
  char Stop_File[MAXNAMELEN+1];   /* Filename for most recent Stop */
  int Start_Line;         /* Line # in Start_File for most recent Start */
  int Stop_Line;          /* Line # in Stop_File for most recent Stop */
  double My_Tot_Time;     /* Sum of stop_time-start_time over all invocations
                             of this timer */
  int Use_Barrier;        /* Flag indicating whether to perform a barrier
                             operation before starting the timer. */
  int Status;             /* Flag indicating status of TimeStruct:
                                > 0  -->  In Use
                                > 2  -->  Running */
  char Name[MAXNAMELEN+1];/* String associated (and printed) with timer info */

#ifdef VAMPIR
  int vt_handle;          /* state handle for vampir traces */
#endif
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

ZTIMER *Zoltan_Timer_Copy(ZTIMER *from)
{
  ZTIMER *to = NULL;

  Zoltan_Timer_Copy_To(&to, from);

  return to;
}

/****************************************************************************/
int Zoltan_Timer_Copy_To(ZTIMER **to, ZTIMER *from)
{
  ZTIMER *toptr = NULL;
  if (!to){
    return ZOLTAN_FATAL;
  }

  if (*to){
    Zoltan_Timer_Destroy(to);
  }

  if (from){
    *to = (ZTIMER *)ZOLTAN_MALLOC(sizeof(ZTIMER));
    toptr = *to;

    toptr->Timer_Flag = from->Timer_Flag;
    toptr->Length = from->Length;
    toptr->NextTimeStruct = from->NextTimeStruct;
    if (toptr->Length > 0){
      toptr->Times = (ZTIMER_TS *)ZOLTAN_MALLOC(sizeof(ZTIMER_TS) * toptr->Length);
      memcpy(toptr->Times, from->Times, sizeof(ZTIMER_TS) * toptr->Length);
    }
    else{
      toptr->Times = NULL;
    }
  }
  
  return ZOLTAN_OK;
}

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
  int use_barrier,      /* Flag indicating whether to perform a 
                           barrier operation before starting the
                           timer. */
  const char *name            /* Name of this timer */
)
{
/* Function that returns the index of the next available Timer timer. */
int ret;
static char *yo = "Zoltan_Timer_Init";

  TESTTIMER(zt, yo);
  
  ret = zt->NextTimeStruct++;

  if (ret >= zt->Length) {
    /* Realloc -- need more individual timers */
    zt->Length += INITLENGTH;
    zt->Times = (ZTIMER_TS *) ZOLTAN_REALLOC(zt->Times, zt->Length * sizeof(ZTIMER_TS));
  }

  Zoltan_Timer_Reset(zt, ret, use_barrier, name);

#ifdef VAMPIR
  if (VT_funcdef(name, VT_NOCLASS, &((zt->Times[ret]).vt_handle)) != VT_OK)
      FATALERROR(yo, "VT_funcdef failed.");
#endif

  return ret;
}


/****************************************************************************/
int Zoltan_Timer_Reset(
  ZTIMER *zt,
  int ts_idx,            /* Index of the timer to reset */
  int use_barrier,       /* Flag indicating whether to perform a 
                            barrier operation before starting the
                            timer. */
  const char *name             /* Name of this timer */
)
{
/* Initialize a timer for INUSE; reset its values to zero. */
static char *yo = "Zoltan_Timer_Reset";
ZTIMER_TS *ts;

  TESTTIMER(zt, yo);
  TESTINDEX(zt, ts_idx, yo);

  ts = &(zt->Times[ts_idx]);

  ts->Status = INUSE;
  ts->Start_Time = 0.;
  ts->Stop_Time = 0.;
  ts->My_Tot_Time = 0.;
  ts->Use_Barrier = use_barrier;
  strncpy(ts->Name, name, MAXNAMELEN);
  ts->Name[MAXNAMELEN] = '\0';
  ts->Start_File[0] = '\0';
  ts->Start_Line = -1;
  ts->Stop_File[0] = '\0';
  ts->Stop_Line = -1;

  return ZOLTAN_OK;
}

/****************************************************************************/
int Zoltan_Timer_ChangeFlag(
  ZTIMER *zt,
  int timer
)
{
static char *yo = "Zoltan_Timer_ChangeFlag";

  TESTTIMER(zt, yo);
  zt->Timer_Flag = timer;
  return ZOLTAN_OK;
}

/****************************************************************************/
int Zoltan_Timer_Start(
  ZTIMER *zt,            /* Ptr to Timer object */
  int ts_idx,            /* Index of the timer to use */
  MPI_Comm comm,         /* Communicator to use for synchronization, 
                            if requested */
  char *filename,        /* Filename of file calling the Start */
  int lineno             /* Line number where Start was called */
)
{
ZTIMER_TS *ts;
static char *yo = "Zoltan_Timer_Start";

  TESTTIMER(zt, yo);
  TESTINDEX(zt, ts_idx, yo);

  ts = &(zt->Times[ts_idx]);
  if (ts->Status > RUNNING)  {
    char msg[256];
    sprintf(msg, 
            "Cannot start timer %d at %s:%d; timer already running from %s:%d.",
            ts_idx, filename, lineno, ts->Start_File, ts->Start_Line);
    FATALERROR(yo, msg)
  }

  ts->Status += RUNNING;
  strncpy(ts->Start_File, filename, MAXNAMELEN);
  ts->Start_Line = lineno;
  if (ts->Use_Barrier)
    MPI_Barrier(comm);

  ts->Start_Time = Zoltan_Time(zt->Timer_Flag);

#ifdef VAMPIR
  if (VT_begin(ts->vt_handle) != VT_OK)
      FATALERROR(yo, "VT_begin failed.");
#endif

  return ZOLTAN_OK;
}

/****************************************************************************/

int Zoltan_Timer_Stop(
  ZTIMER *zt,            /* Ptr to Timer object */
  int ts_idx,            /* Index of the timer to use */
  MPI_Comm comm,         /* Communicator to use for synchronization, 
                            if requested */
  char *filename,        /* Filename of file calling the Stop */
  int lineno             /* Line number where Stop was called */
)
{
/* Function to stop a timer and accrue its information */
ZTIMER_TS *ts;
static char *yo = "Zoltan_Timer_Stop";
double my_time;

  TESTTIMER(zt, yo);
  TESTINDEX(zt, ts_idx, yo);

  ts = &(zt->Times[ts_idx]);
  if (ts->Status < RUNNING) {
    if (ts->Stop_Line == -1)
      FATALERROR(yo, "Cannot stop timer; timer never started.")
    else {
      char msg[256];
      sprintf(msg, 
              "Cannot stop timer %d at %s:%d; "
              "timer already stopped from %s:%d.",
              ts_idx, filename, lineno, ts->Stop_File, ts->Stop_Line);
      FATALERROR(yo, msg)
    }
  }

#ifdef VAMPIR
  if (VT_end(ts->vt_handle) != VT_OK)
      FATALERROR(yo, "VT_end failed.");
#endif

  if (ts->Use_Barrier)
    MPI_Barrier(comm);
  ts->Stop_Time = Zoltan_Time(zt->Timer_Flag);
  ts->Status -= RUNNING;
  ts->Stop_Line = lineno;
  strncpy(ts->Stop_File, filename, MAXNAMELEN);
  my_time = ts->Stop_Time - ts->Start_Time;

  ts->My_Tot_Time += my_time;

  return ZOLTAN_OK;
}


/****************************************************************************/
int Zoltan_Timer_Print(
  ZTIMER *zt,
  int ts_idx,
  int proc,  /* Rank of the processor (in comm) that should print the data. */
  MPI_Comm comm,
  FILE *fp
)
{
/* Accrues a single timer's values across a communicator and prints 
 * its information.  This function must be called by all processors
 * within the communicator.  
 */
static char *yo = "Zoltan_Timer_Print";
ZTIMER_TS *ts;
int my_proc, nproc;
int restart = 0;
double max_time;
double min_time;
double sum_time;

  TESTTIMER(zt, yo);
  TESTINDEX(zt, ts_idx, yo);
  MPI_Comm_rank(comm, &my_proc);
  MPI_Comm_size(comm, &nproc);

  ts = &(zt->Times[ts_idx]);
  if (ts->Status > RUNNING)  {
    /* Timer is running; stop it before printing the times.
     * Don't want to include print times in timer.
     */
    restart = 1;
    ZOLTAN_TIMER_STOP(zt, ts_idx, comm);
  }

  MPI_Allreduce(&(ts->My_Tot_Time), &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&(ts->My_Tot_Time), &min_time, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&(ts->My_Tot_Time), &sum_time, 1, MPI_DOUBLE, MPI_SUM, comm);

  if (proc == my_proc) 
    fprintf(fp,
            "%3d ZOLTAN_TIMER %3d %23s:  MyTime %7.4f  "
            "MaxTime %7.4f  MinTime %7.4f  AvgTime %7.4f\n",
            proc, ts_idx, ts->Name, ts->My_Tot_Time, 
            max_time, min_time, sum_time/nproc);

  if (restart) {
    /* We stopped the timer for printing; restart it now. */
    ZOLTAN_TIMER_START(zt, ts_idx, comm);
  }

  return ZOLTAN_OK;
}

/****************************************************************************/
int Zoltan_Timer_PrintAll(
  ZTIMER *zt,
  int proc,    /* Rank of the processor (in comm) that should print the data. */
  MPI_Comm comm, 
  FILE *fp
)
{
/* Function to print all timer information */
static char *yo = "Zoltan_Timer_PrintAll";
int i, ierr = ZOLTAN_OK;

  TESTTIMER(zt, yo);
  for (i = 0; i < zt->NextTimeStruct; i++) 
    if ((ierr = Zoltan_Timer_Print(zt, i, proc, comm, fp)) != ZOLTAN_OK)
      break;

  return ierr;
}

/****************************************************************************/
void Zoltan_Timer_Destroy(
  ZTIMER **zt
)
{
/* Destroy a Timer object */
  if (*zt != NULL) {
    ZOLTAN_FREE(&((*zt)->Times));
    ZOLTAN_FREE(zt);
  }
}

/****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
