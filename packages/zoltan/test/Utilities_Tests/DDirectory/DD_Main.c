/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


/*  Program tests the Zoltan Distributed Directory software (stand-alone mode).
*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "zoltan_dd.h"
#include "zoltan_mem.h"
#include "zoltan_util.h"
#include "zoltan_id.h"
#include "zoltan_align.h"
#include "zoltan_timer.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#define DD_TEST_NAME_NUMERIC 100  /* arbitrary starting value for defines */
#define NO_PROC -1


typedef struct {
   int old_owner ;
   int new_owner ;      /* old != new ==> simultated data move    */
   int partition ;
   ZOLTAN_ID_TYPE id[1] ;   /* malloc'd beyond end of struct              */
                            /* ZOLTAN_ID_TYPE lid[]  at gid + glen        */
                            /* ZOLTAN_ID_TYPE user[] at gid + glen + llen */
} Data ;


typedef struct {
   int count ;          /* number of global IDs to simulate            */
   int pdelete ;        /* fraction of GIDs removed per cycle x100     */
   int pmove ;          /* fraction of GIDs moved per cycle x100       */
   int pscatter ;       /* fraction of GIDs scattered per cycle x 100  */
   int nloops ;         /* number of times to loop though simulation   */
   int rseed ;          /* seed for random number generator            */

   int glen ;           /* global ID length    */
   int llen ;           /* local ID length     */
   int ulen ;           /* user data ID length */
   int tlen ;           /* hash table length   */
   int slen ;           /* length of Data * id stuff */

   int name_scheme ;
   int debug_level ;
} Param ;


static int get_params (Param *p) ;

static int initialize_data (Param *p, char *store, int nproc) ;

/****************************************************************************/
#define MAXNAMELEN 31
#define INITLENGTH 30

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

/***********************************************************************************/

#define MACRO_TIMER_START(arg, message, sync)\
  {if (timer[arg] == -1)\
     timer[arg] = Zoltan_Timer_Init(zt, sync, message);\
   ZOLTAN_TIMER_START(zt, timer[arg], MPI_COMM_WORLD);}

#define MACRO_TIMER_STOP(arg)\
   {ZOLTAN_TIMER_STOP(zt, timer[arg], MPI_COMM_WORLD);}

int main (int argc, char *argv[])
{
   Zoltan_DD_Directory *dd;
   ZTIMER *zt;

   int        myproc;         /* MPI rank, my processor's MPI ID      */
   int        nproc;          /* MPI size, number of processors       */
   ZOLTAN_ID_PTR  glist = NULL;   /* list of GIDs                     */
   ZOLTAN_ID_PTR  llist = NULL;   /* list of LIDs                     */
   ZOLTAN_ID_PTR  ulist = NULL;   /* list of user data of type 
                                      ZOLTAN_ID_TYPE                   */
   int       *plist = NULL;    /* list of partitions                   */
   int       *olist = NULL;    /* list of owners                       */
   Param      param ;          /* program's adjustable parameters      */
   char      *store = NULL;    /* non directory storage of test data   */
   Data      *data  = NULL;    /* pointer into store                   */

   /* these are all temporary variables */
   int   new_owner;
   int   count;
   int   i;
   char *p;
   int   err;
   int   errcount;
   int   loop;
   char str[100];      /* for building output messages */
   static int timer[7] = {-1,-1,-1,-1,-1,-1,-1};
   char *yo = "DD_Main";


   /* initialize MPI communications, ignore errors */
   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &myproc);
   MPI_Comm_size (MPI_COMM_WORLD, &nproc);

   get_params (&param);     /* read input parameters */

   zt = Zoltan_Timer_Create(1);
   MACRO_TIMER_START(5, "Total time", 0);

   MACRO_TIMER_START(0, "DD_Create time", 0);
   err = Zoltan_DD_Create (&dd, MPI_COMM_WORLD, param.glen, param.llen,
    param.ulen*sizeof(ZOLTAN_ID_TYPE), param.tlen, param.debug_level);
   if (err != ZOLTAN_OK)
      ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Create");
   MACRO_TIMER_STOP(0);

   param.slen = sizeof (Data) + sizeof(ZOLTAN_ID_TYPE)
              * (param.glen + param.llen + param.ulen);
   param.slen = Zoltan_Align(param.slen);
   store = (char*) ZOLTAN_MALLOC (param.count * param.slen);

   /* allocate storage for various lists */
   glist = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC(sizeof(ZOLTAN_ID_TYPE) * param.count
         * param.glen);
   plist = (int*) ZOLTAN_MALLOC (sizeof(int) * param.count);
   olist = (int*) ZOLTAN_MALLOC (sizeof(int) * param.count);
   if (param.llen)
      llist = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC (sizeof (ZOLTAN_ID_TYPE)
            * param.count * param.llen);
   if (param.ulen)
      ulist = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC (sizeof (ZOLTAN_ID_TYPE)
            * param.count * param.ulen) ;


   if (store == NULL || glist == NULL || (param.llen != 0 && llist == NULL)
    || (param.ulen != 0 && ulist == NULL) || plist == NULL || olist == NULL)  {
       ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to malloc storage lists");
       return ZOLTAN_MEMERR;
   }

   initialize_data (&param, store, nproc);

   /* create & update directory with myproc's initial simulated GIDs */
   count = 0;
   for (p = store; p < store + param.count * param.slen; p += param.slen)
      if (((Data *)p)->new_owner == myproc) {
         ZOLTAN_SET_ID (param.glen, glist + count * param.glen, ((Data*)p)->id);
         if (param.llen)
            ZOLTAN_SET_ID (param.llen, llist + count * param.llen, 
             ((Data *)p)->id + param.glen);
         if (param.ulen)
            ZOLTAN_SET_ID (param.ulen, ulist + count * param.ulen, 
             ((Data *)p)->id + (param.glen + param.llen));
         plist [count] = ((Data *)p)->partition;
         ++count;
      }
   MACRO_TIMER_START (1, "DD_Update timer", 0);
printf("KDDKDD Update1 %d\n", count);
   err = Zoltan_DD_Update (dd, glist, llist, (char*)ulist, plist, count);
   if (err != ZOLTAN_OK)
      ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Update");
   MACRO_TIMER_STOP(1);

   i = 0;
   for (p = store; p < store + param.count * param.slen; p += param.slen)  {
      ZOLTAN_SET_ID (param.glen, glist + i * param.glen, ((Data *)p)->id);
      ++i;
   }
printf("KDDKDD Find1\n");
   MACRO_TIMER_START(2, "DD_Find timer", 0);
   err = Zoltan_DD_Find (dd, glist, llist, (char*)ulist, plist, 
                         param.count, olist);
   if (err != ZOLTAN_OK)
      ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Find");
   MACRO_TIMER_STOP(2);

   errcount = 0;
   for (i = 0; i < param.count; i++)
       if (olist[i] != ((Data *) (store + i * param.slen))->new_owner)
          ++errcount;
   if (errcount)
      sprintf (str, "FIRST TEST FAILED, errcount is %d", errcount);
   else
      sprintf (str, "FIRST TEST SUCCESSFUL");
   ZOLTAN_PRINT_INFO (myproc, yo, str);

   /* repeatedly simulate moving "dots" around the system */
   for (loop = 0; loop < param.nloops; loop++) {
       for (p = store; p < store + param.count * param.slen; p += param.slen)
          ((Data*) p)->old_owner = ((Data*) p)->new_owner;

       /* randomly exchange some dots and randomly reassign others */
       for (i = 0; i < (param.pmove * param.count)/100; i++)  {
          Data *d1, *d2;
          d1 = (Data*) (store + param.slen * (rand() % param.count));
          d2 = (Data*) (store + param.slen * (rand() % param.count));
          new_owner     = d1->new_owner;
          d1->new_owner = d2->new_owner;
          d2->new_owner = new_owner;
       }

       for (i = 0; i < (param.count * param.pscatter)/100; i++)
          ((Data*) (store + param.slen *(rand() % param.count)))->new_owner
           = rand() % nproc;

       /* determine which new GIDs myproc gained, and update directory */
       count = 0;
       for (p = store; p < store + param.count * param.slen; p += param.slen)
           if (((Data*)p)->new_owner == myproc)  {
              ((Data*)p)->id[param.glen] = count;  /* set LID */
              ZOLTAN_SET_ID (param.glen, glist + count * param.glen, 
               ((Data *)p)->id);
              if (param.llen)
                 ZOLTAN_SET_ID (param.llen, llist + count * param.llen, 
                  ((Data*)p)->id + param.glen);
              if (param.ulen)
                  ZOLTAN_SET_ID (param.ulen, ulist + count * param.ulen, 
                   ((Data*)p)->id + (param.glen + param.llen));
              plist [count] = ((Data *)p)->partition;
              ++count;
           }
printf("KDDKDD Update2\n");
       MACRO_TIMER_START(1, "DD_Update", 0);
       err = Zoltan_DD_Update (dd, glist, llist, (char*)ulist, plist, count);
       if (err != ZOLTAN_OK)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Update");
       MACRO_TIMER_STOP(1);

       /* use directory to "find" GIDs  */
       i = 0;
       for (p = store; p < store + param.count * param.slen; p += param.slen) {
          ZOLTAN_SET_ID (param.glen, glist + i * param.glen, ((Data *)p)->id);
          ++i;
       }
 
printf("KDDKDD Find2\n");
       MACRO_TIMER_START(2, "DD_Find timer", 0);
       if (loop % 5 == 0)
          err = Zoltan_DD_Find(dd,glist,NULL,(char*)ulist,
                               plist,param.count,olist);
       else if (loop % 7 == 0)
          err = Zoltan_DD_Find(dd,glist,llist,NULL,plist,param.count,olist);
       else if (loop % 9 == 0)
          err = Zoltan_DD_Find(dd,glist,llist,(char*)ulist,
                               NULL,param.count,olist);
       else if (loop % 2 == 0) 
          err = Zoltan_DD_Find(dd,glist,llist,(char*)ulist,
                               plist,param.count,NULL);
       else
           err = Zoltan_DD_Find(dd,glist,llist,(char*)ulist,
                               plist,param.count,olist);
       if (err != ZOLTAN_OK)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Find");
       MACRO_TIMER_STOP(2);

       if (loop % 2)  {
          errcount = 0 ;
          for (i = 0; i < param.count; i++)
             if (olist[i] != ((Data *)(store + i * param.slen))->new_owner)
                ++errcount;
          if (errcount)
             sprintf (str,"LOOP %d TEST FAILED, errcount is %d",loop,errcount);
          else
             sprintf (str, "LOOP %d TEST SUCCESSFUL", loop);
          ZOLTAN_PRINT_INFO (myproc, yo, str) ;
       }
       else  {
          sprintf (str, "LOOP %d Completed", loop);
          ZOLTAN_PRINT_INFO (myproc, yo, str);
       }

       /* randomly remove a percentage of GIDs from the directory */
       count = 0;
       for (i = 0; i < (param.count * param.pdelete)/100; i++)  {
          data = (Data*) (store + param.slen * (rand() % param.count));
          if (data->new_owner == myproc)  {
             ZOLTAN_SET_ID (param.glen, glist + count * param.glen, data->id);
             ++count;
          }
          data->new_owner = NO_PROC;
       }
printf("KDDKDD Remove1\n");
       MACRO_TIMER_START(3, "DD_Remove timer", 0);
       err = Zoltan_DD_Remove (dd, glist, count);
       MACRO_TIMER_STOP(3);
       if (err != ZOLTAN_OK)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Remove");

       /* update directory (put directory entries back in) */
       for (p = store; p < store + param.count * param.slen; p += param.slen)
          if (((Data*)p)->new_owner == NO_PROC)
              ((Data*)p)->new_owner = loop % nproc; /* place in new location */

       count = 0;
       for (p = store; p < store + param.count * param.slen; p += param.slen)
          if (((Data*)p)->new_owner == myproc)  {
              ZOLTAN_SET_ID(param.glen,glist+count*param.glen,((Data*)p)->id);
              if (param.llen)
                 ZOLTAN_SET_ID (param.llen, llist + count * param.llen, 
                  ((Data*)p)->id + param.glen);
              if (param.ulen)
                  ZOLTAN_SET_ID(param.ulen, ulist + count * param.ulen, 
                   ((Data *)p)->id + (param.glen + param.llen));
              plist [count] = ((Data*)p)->partition;
              ++count;
           }
printf("KDDKDD Update3\n");
       MACRO_TIMER_START(1, "DD_Update timer", 0);
       err = Zoltan_DD_Update (dd, glist, NULL, NULL, NULL, count);
       if (err != ZOLTAN_OK)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Update");
       MACRO_TIMER_STOP(1);
   }

   /* now Find directory info for GIDs myproc now owns and validate */
   count = 0;
   for (i = 0; i < param.count; i++) {
      data = (Data *) (store + i * param.slen);
      if (data->new_owner == myproc)  {
         ZOLTAN_SET_ID (param.glen, glist + count * param.glen, data->id);
         ++count;
      }
   }
printf("KDDKDD Find3\n");
   MACRO_TIMER_START(2, "DD_Find", 0);
   err = Zoltan_DD_Find (dd, glist, NULL, NULL, NULL, count, olist);
   if (err != ZOLTAN_OK)
      ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Find");
   MACRO_TIMER_STOP(2);

   errcount = 0;
   for (i = 0; i < count; i++)
       if (olist[i] != myproc)
          ++errcount;
   if (errcount)  {
      sprintf (str, "errcount is %d", errcount);
      ZOLTAN_PRINT_ERROR (myproc, yo, str);
   }
   else
      ZOLTAN_PRINT_INFO (myproc, yo, "TEST SUCCESSFUL");

printf("KDDKDD Stats\n");
   Zoltan_DD_Stats (dd);

   /* done, now free memory, stop MPI & directory, return */
   ZOLTAN_FREE (&store);
   ZOLTAN_FREE (&glist);
   ZOLTAN_FREE (&plist);
   ZOLTAN_FREE (&olist);

   if (param.llen) ZOLTAN_FREE (&llist);
   if (param.ulen) ZOLTAN_FREE (&ulist);

   ZOLTAN_PRINT_INFO (myproc, yo, "Completing program");
   MACRO_TIMER_START(4, "DD_Destroy", 0);
printf("KDDKDD Destroy\n");
   Zoltan_DD_Destroy (&dd);
   MACRO_TIMER_STOP(4);
   MACRO_TIMER_STOP(5);

   Zoltan_Timer_PrintAll(zt, 0, MPI_COMM_WORLD, stderr);


   MPI_Finalize ();
   return ZOLTAN_OK;
}



static int initialize_data (Param *param, char *store, int nproc)
{
   int i;
   int j;
   Data *data;

   if (param->name_scheme == DD_TEST_NAME_NUMERIC)  {
      for (i = 0; i < param->count; i++)  {
         data = (Data*) (store + i * param->slen);
         for (j = 0; j < param->glen; j++)
            data->id[j] = (ZOLTAN_ID_TYPE) i;

         for (j = 0; j < param->llen; j++)
            data->id[j + param->glen] = (ZOLTAN_ID_TYPE) (i % nproc);

         for (j = 0; j < param->ulen; j++)
            data->id[j+param->glen+param->llen] = (ZOLTAN_ID_TYPE)(i%nproc);

         data->partition = i % 7;
         data->new_owner = i % nproc;
         data->old_owner = i % nproc;
      }
   }

   srand ((unsigned int) param->rseed);

   return ZOLTAN_OK;
}




static int get_params (Param *param)
{
    /* establish default values */
    param->count       = 250000; 
    param->pmove       = 10;
    param->pdelete     = 3;
    param->pscatter    = 4;
    param->nloops      = 20;

    param->glen        = 5;
    param->llen        = 5;
    param->ulen        = 10;
    param->tlen        = 0;

    param->debug_level = 0; 
    param->name_scheme = DD_TEST_NAME_NUMERIC;

    param->rseed = 31415926;

    /* read external file for optional non default parameter values */

    return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
