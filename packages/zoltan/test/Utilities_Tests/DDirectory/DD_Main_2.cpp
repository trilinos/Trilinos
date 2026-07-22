// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*  Program tests the Zoltan Distributed Directory software (stand-alone mode).
*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "zoltan_dd_cpp.h"
#include "zoltan_mem.h"
#include "zoltan_util.h"
#include "zoltan_id.h"
#include "zoltan_align.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#define DD_TEST_NAME_NUMERIC 100  /* arbitrary starting value for defines */
#define NO_PROC -1


typedef struct
{
   int old_owner;
   int new_owner;      /* old != new ==> simultated data move    */
   int partition;
   ZOLTAN_ID_TYPE id[1];   /* malloc'd beyond end of struct              */
                            /* ZOLTAN_ID_TYPE lid[]  at gid + glen        */
                            /* ZOLTAN_ID_TYPE user[] at gid + glen + llen */
} Data;


typedef struct
{
   int count;          /* number of global IDs to simulate          */
   int pdelete;        /* fraction of GIDs removed per cycle x100   */
   int pmove;          /* fraction of GIDs moved per cycle x100     */
   int pscatter;       /* fraction of GIDs scattered per cycle x 100 */
   int nloops;         /* number of times to loop though simulation */
   int rseed;          /* seed for random number generator          */

   int glen;           /* global ID length    */
   int llen;           /* local ID length     */
   int ulen;           /* user data ID length */
   int tlen;           /* hash table length   */
   int slen;           /* length of Data * id stuff */

   int name_scheme;
   int debug_level;
} Param;


static int get_params (Param *p);

static int initialize_data (Param *p, char *store, int nproc);







int main (int argc, char *argv[])
{
   Zoltan_DD *dd;

   int        myproc;         /* MPI rank, my processor's MPI ID      */
   int        nproc;          /* MPI size, number of processors       */
   ZOLTAN_ID_PTR  glist = NULL;   /* list of GIDs                     */
   ZOLTAN_ID_PTR  llist = NULL;   /* list of LIDs                     */
   ZOLTAN_ID_PTR  ulist = NULL;   /* list of user data of type 
                                      ZOLTAN_ID_TYPE                   */
   int       *plist = NULL;   /* list of partitions                   */
   int       *olist = NULL;    /* list of owners                       */
   Param      param;          /* program's adjustable parameters      */
   char      *store = NULL;   /* non directory storage of test data   */
   Data      *data  = NULL;   /* pointer into store                   */

   /* these are all temporary variables:  */
   int   new_owner;
   int   count;
   int   i;
   char *p;
   int   err;
   int   errcount;
   int   loop;


   char str[100];      /* for building output messages */
   char *yo = "DD_Main";


   /* initialize MPI communications, ignore errors */
   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &myproc);
   MPI_Comm_size (MPI_COMM_WORLD, &nproc);

   get_params (&param);     /* read input parameters */

   dd = new Zoltan_DD( MPI_COMM_WORLD, param.glen, param.llen,
    param.ulen*sizeof(ZOLTAN_ID_TYPE), param.tlen, param.debug_level);

   param.slen = sizeof (Data) 
              + (param.glen + param.llen + param.ulen) * sizeof(ZOLTAN_ID_TYPE);
   param.slen = Zoltan_Align(param.slen);
   store = (char *) ZOLTAN_MALLOC (param.count * param.slen);

   /* allocate storage for various lists */
   glist = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC (sizeof(ZOLTAN_ID_TYPE) * param.count
                                                             * param.glen);
   plist = (int *)     ZOLTAN_MALLOC (sizeof (int)        * param.count);
   olist = (int *)     ZOLTAN_MALLOC (sizeof (int)        * param.count);
   if (param.llen != 0)
      llist = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC (sizeof (ZOLTAN_ID_TYPE) * 
                                             param.count * param.llen);
   if (param.ulen != 0)
      ulist = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC (sizeof (ZOLTAN_ID_TYPE) * 
                                             param.count * param.ulen);


   if (store == NULL || glist == NULL
    || (param.llen != 0 && llist == NULL)
    || (param.ulen != 0 && ulist == NULL)
    || plist == NULL || olist == NULL) {
       ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to malloc storage lists");
       return 1;
    }

   initialize_data (&param, store, nproc);

   /* create & update directory with myproc's initial simulated GIDs */
   count = 0;
   for (p = store; p < store + param.count * param.slen; p += param.slen)
      if (((Data *)p)->new_owner == myproc) {
         ZOLTAN_SET_ID (param.glen, glist + count * param.glen, 
                        ((Data *)p)->id);
         if (param.llen)
            ZOLTAN_SET_ID (param.llen, llist + count * param.llen, 
                           ((Data *)p)->id + param.glen);
         if (param.ulen)
            ZOLTAN_SET_ID (param.ulen, ulist + count * param.ulen, 
                           ((Data *)p)->id + (param.glen + param.llen));
         plist [count] = ((Data *)p)->partition;
         count++;
      }
   err = dd->Update (glist, llist, (char *)ulist, plist, count);
   if (err != ZOLTAN_DD_NORMAL_RETURN)
      ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Update");

   i = 0;
   for (p = store; p < store + param.count * param.slen; p += param.slen) {
      ZOLTAN_SET_ID (param.glen, glist + i * param.glen, ((Data *)p)->id);
      i++;
   }
   err = dd->Find (glist, llist, (char *)ulist, plist, param.count, olist);
   if (err != ZOLTAN_DD_NORMAL_RETURN)
      ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Find");

   errcount = 0;
   for (i = 0; i < param.count; i++)
       if (olist[i] != ((Data *) (store + i * param.slen))->new_owner)
          errcount++;
   if (errcount > 0)
      sprintf (str, "FIRST TEST FAILED, errcount is %d", errcount);
   else
      sprintf (str, "FIRST TEST SUCCESSFUL");
   ZOLTAN_PRINT_INFO (myproc, yo, str);

   /* repeatedly simulate moving "dots" around the system */
   for (loop = 0; loop < param.nloops; loop++) {
       for (p = store; p < store + param.count * param.slen; p += param.slen)
          ((Data *) p)->old_owner = ((Data *) p)->new_owner;

       /* randomly exchange some dots and randomly reassign others */
       for (i = 0; i < (param.pmove * param.count)/100; i++) {
          Data *d1, *d2;
          d1 = (Data *) (store + param.slen * (rand() % param.count));
          d2 = (Data *) (store + param.slen * (rand() % param.count));
          new_owner     = d1->new_owner;
          d1->new_owner = d2->new_owner;
          d2->new_owner = new_owner;
       }

       for (i = 0; i < (param.count * param.pscatter)/100; i++)
          ((Data *) (store + param.slen *(rand() % param.count)))->new_owner
            = rand() % nproc;

       /* determine which new GIDs myproc gained, and update directory */
       count = 0;
       for (p = store; p < store + param.count * param.slen; p += param.slen)
           if (((Data *)p)->new_owner == myproc) {
              ((Data *)p)->id[param.glen] = count;  /* set LID */
              ZOLTAN_SET_ID (param.glen, glist + count * param.glen, 
                             ((Data *)p)->id);
              if (param.llen)
                 ZOLTAN_SET_ID (param.llen, llist + count * param.llen, 
                                ((Data *)p)->id + param.glen);
              if (param.ulen)
                  ZOLTAN_SET_ID (param.ulen, ulist + count * param.ulen, 
                                 ((Data *)p)->id + (param.glen + param.llen));
              plist [count] = ((Data *)p)->partition;
              count++;
           }
       err = dd->Update (glist, llist, (char *)ulist, plist, count);
       if (err != ZOLTAN_DD_NORMAL_RETURN)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Update");

       /* use directory to "find" GIDs  */
       i = 0;
       for (p = store; p < store + param.count * param.slen; p += param.slen) {
          ZOLTAN_SET_ID (param.glen, glist + i * param.glen, ((Data *)p)->id);
          i++;
       }
       err = dd->Find (glist, llist, (char *)ulist, plist, param.count, olist);
       if (err != ZOLTAN_DD_NORMAL_RETURN)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Find");

       errcount = 0;
       for (i = 0; i < param.count; i++)
          if (olist[i] != ((Data *)(store + i * param.slen))->new_owner)
             errcount++;
       if (errcount > 0)
          sprintf (str, "LOOP %d TEST FAILED, errcount is %d", loop, errcount);
       else
          sprintf (str, "LOOP %d TEST SUCCESSFUL", loop);
       ZOLTAN_PRINT_INFO (myproc, yo, str);


       /* randomly remove a percentage of GIDs from the directory */
       count = 0;
       for (i = 0; i < (param.count * param.pdelete)/100; i++) {
          data = (Data *) (store + param.slen * (rand() % param.count));

          if (data->new_owner == myproc) {
             ZOLTAN_SET_ID (param.glen, glist + count * param.glen, data->id);
             count++;
          }
          data->new_owner = NO_PROC;
       }
       err = dd->Remove (glist, count);
       if (err != ZOLTAN_DD_NORMAL_RETURN)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Remove");


       /* update directory (put directory entries back in) */
       for (p = store; p < store + param.count * param.slen; p += param.slen)
          if (((Data *)p)->new_owner == NO_PROC)
              ((Data *)p)->new_owner = loop % nproc; /* place in new location */

       count = 0;
       for (p = store; p < store + param.count * param.slen; p += param.slen)
          if (((Data *)p)->new_owner == myproc) {
              ZOLTAN_SET_ID (param.glen, glist + count * param.glen, 
                             ((Data *)p)->id);
              if (param.llen)
                 ZOLTAN_SET_ID (param.llen, llist + count * param.llen, 
                                ((Data *)p)->id + param.glen);
              if (param.ulen)
                  ZOLTAN_SET_ID (param.ulen, ulist + count * param.ulen, 
                                 ((Data *)p)->id + (param.glen + param.llen));
              plist [count] = ((Data *)p)->partition;
              count++;
           }
       err = dd->Update (glist, NULL, NULL, NULL, count);
       if (err != ZOLTAN_DD_NORMAL_RETURN)
          ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Update");
   }

   /* now Find directory info for GIDs myproc now owns and validate */
   count = 0;
   for (i = 0; i < param.count; i++) {
      data = (Data *) (store + i * param.slen);
      if (data->new_owner == myproc) {
         ZOLTAN_SET_ID (param.glen, glist + count * param.glen, data->id);
         count++;
      }
   }
   err = dd->Find (glist, NULL, NULL, NULL, count, olist);
   if (err != ZOLTAN_DD_NORMAL_RETURN)
      ZOLTAN_PRINT_ERROR (myproc, yo, "Failed return from DD Find");

   errcount = 0;
   for (i = 0; i < count; i++)
       if (olist[i] != myproc)
          errcount++;
   if (errcount > 0) {
      sprintf (str, "errcount is %d", errcount);
      ZOLTAN_PRINT_ERROR (myproc, yo, str);
   }
   else
      ZOLTAN_PRINT_INFO (myproc, yo, "TEST SUCCESSFUL");

   dd->Stats ();

   /* done, now free memory, stop MPI & directory, return */
   ZOLTAN_FREE (&store);
   ZOLTAN_FREE (&glist);
   ZOLTAN_FREE (&plist);
   ZOLTAN_FREE (&olist);

   if (param.llen != 0)   ZOLTAN_FREE (&llist);
   if (param.ulen != 0)   ZOLTAN_FREE (&ulist);

   ZOLTAN_PRINT_INFO (myproc, yo, "Completing program");

   delete dd;

   MPI_Finalize ();
   return 0;
}







static int initialize_data (Param *param, char *store, int nproc)
{
   int i;
   int j;
   Data *data;

   if (param->name_scheme == DD_TEST_NAME_NUMERIC)
      {
      for (i = 0; i < param->count; i++)
        {
        data = (Data *) (store + i * param->slen);
        for (j = 0; j < param->glen; j++)
           data->id[j] = (ZOLTAN_ID_TYPE) i;

        for (j = 0; j < param->llen; j++)
           data->id[j + param->glen] = (ZOLTAN_ID_TYPE) (i % nproc);

        for (j = 0; j < param->ulen; j++)
           data->id[j + param->glen + param->llen] = 
              (ZOLTAN_ID_TYPE) (i % nproc);

        data->partition = i % 7;
        data->new_owner = i % nproc;
        data->old_owner = i % nproc;
        }
      }

   srand ((unsigned int) param->rseed);

   return 0;
}




static int get_params (Param *param)
{
    /* establish default values */
/*  param->count       = 250000;  KDD */ /* 1500000; */
    param->count       = 2000; /* 1500000; */
    param->pmove       = 10;
    param->pdelete     = 3;
    param->pscatter    = 4;
/*  param->nloops      = 200; KDD */
    param->nloops      = 20;

    param->glen        = 5;
    param->llen        = 5;
    param->ulen        = 10;
    param->tlen        = 5000;

/*  param->debug_level = 1;  KDD */
    param->debug_level = 3;
    param->name_scheme = DD_TEST_NAME_NUMERIC;

    param->rseed = 31415926;

    /* read external file for optional non default parameter values */

    return 0;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
