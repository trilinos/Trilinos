/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
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
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "irb_const.h"
#include "all_allo_const.h"
#include "par_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "create_proc_list_const.h"
#include "comm_const.h"
#include "ha_const.h"

/* Inertial recursive bisection (IRB) load balancing routine operates on
   "dots" as defined in irb.h */

/* Notes:
   dots are balanced across procs by weight (if used)
   on return, proc owns dotnum "dots" in dense array of max-length dotmax
   input weights (if used) are real numbers > 0.0
   can extend "irb_dot" data structure in calling program, see irb_const.h
   returned tree only contains one cut on each proc,
   need to do MPI_Allgather if wish to collect it on all procs */

#define TINY   1.0e-6

#define IRB_DEFAULT_OVERALLOC 1.0

static void IRB_check(LB *, struct irb_dot *, int, int);
static void IRB_stats(LB *, double, struct irb_dot *,int, double *, int *, int);

static int irb_fn(LB *, int *, LB_GID **, LB_LID **, int **, double,
               int, int, int, int);

/*  IRB_CHECK = 0  No consistency check on input or results */
/*  IRB_CHECK = 1  Check input weights and final results for consistency */
static int IRB_CHECK = 1;

/*  IRB_OUTPUT_LEVEL = 0  No statistics logging */
/*  IRB_OUTPUT_LEVEL = 1  Log times and counts, print summary */
/*  IRB_OUTPUT_LEVEL = 2  Log times and counts, print for each proc */
static int IRB_OUTPUT_LEVEL = 1;

/*  Parameters structure for IRB method.  Used in  */
/*  LB_Set_IRB_Param and LB_IRB.                      */
static PARAM_VARS IRB_params[] = {
                  { "IRB_OVERALLOC", NULL, "DOUBLE" },
                  { "IRB_CHECK", NULL, "INT" },
                  { "IRB_OUTPUT_LEVEL", NULL, "INT" },
                  { "IRB_KEEP_CUTS", NULL, "INT" },
                  { NULL, NULL, NULL } };

/*---------------------------------------------------------------------------*/

int LB_Set_IRB_Param(
     char *name,                 /* name of variable */
     char *val                   /* value of variable */
)
{
     int status;
     PARAM_UTYPE result;         /* value returned from Check_Param */
     int index;                  /* index returned from Check_Param */

     status = LB_Check_Param(name, val, IRB_params, &result, &index);

     return(status);
}

/*---------------------------------------------------------------------------*/

int LB_irb(
     LB *lb,                     /* The load-balancing structure with info for
                                    the IRB balancer.                   */
     int *num_import,            /* Number of non-local objects assigned to
                                    this processor in the new decomposition. */
     LB_GID **import_global_ids, /* Returned value:  array of global IDs for
                                    non-local objects in this processor's new
                                    decomposition.                           */
     LB_LID **import_local_ids,  /* Returned value:  array of local IDs for
                                    non-local objects in this processor's new
                                    decomposition.                           */
     int **import_procs,         /* Returned value:  array of processor IDs for
                                    processors owning the non-local objects in
                                    this processor's new decomposition.      */
     int *num_export,            /* Not computed, set to -1 */
     LB_GID **export_global_ids, /* Not computed. */
     LB_LID **export_local_ids,  /* Not computed. */
     int **export_procs          /* Not computed. */
)
{
     /* Wrapper routine to set parameter values and call the real irb. */
     double overalloc;           /* amount to overallocate by when realloc
                                    of dot array must be done.
                                    1.0 = no extra; 1.5 = 50% extra; etc. */
     int wgtflag;                /* (0) do not (1) do use weights.
                                    Multidimensional weights not supported */
     int check;                  /* Check input & output for consistency? */
     int stats;                  /* Print timing & count summary? */
     int gen_tree;               /* (0) don't (1) generate whole treept to use
                                    later for point and box drop. */

     LB_Bind_Param(IRB_params, "IRB_OVERALLOC", (void *) &overalloc);
     LB_Bind_Param(IRB_params, "IRB_CHECK", (void *) &check);
     LB_Bind_Param(IRB_params, "IRB_OUTPUT_LEVEL", (void *) &stats);
     LB_Bind_Param(IRB_params, "IRB_KEEP_CUTS", (void *) &gen_tree);

     overalloc = IRB_DEFAULT_OVERALLOC;
     check = IRB_CHECK;
     stats = IRB_OUTPUT_LEVEL;
     gen_tree = 0;
     wgtflag = (lb->Obj_Weight_Dim > 0); /* Multidim. weights not accepted */

     LB_Assign_Param_Vals(lb->Params, IRB_params, lb->Debug_Level, lb->Proc,
                          lb->Debug_Proc);

     *num_export = -1;  /* We don't compute the export map. */

     return(irb_fn(lb, num_import, import_global_ids, import_local_ids,
                import_procs, overalloc, wgtflag, check, stats, gen_tree));
}

/*---------------------------------------------------------------------------*/

static int irb_fn(
     LB *lb,                     /* The load-balancing structure with info for
                                    the IRB balancer. */
     int *num_import,            /* Number of non-local objects assigned to
                                    this processor in the new decomposition. */
     LB_GID **import_global_ids, /* Returned value:  array of global IDs for
                                    non-local objects in this processor's new
                                    decomposition. */
     LB_LID **import_local_ids,  /* Returned value:  array of local IDs for
                                    non-local objects in this processor's new
                                    decomposition. */
     int **import_procs,         /* Returned value:  array of processor IDs for
                                    processors owning the non-local objects in
                                    this processor's new decomposition. */
     double overalloc,           /* amount to overallocate by when realloc
                                    of dot array must be done.
                                    1.0 = no extra; 1.5 = 50% extra; etc. */
     int wgtflag,                /* (0) do not (1) do use weights.
                                    Multidimensional weights not supported */
     int check,                  /* Check input & output for consistency? */
     int stats,                  /* Print timing & count summary? */
     int gen_tree                /* (0) do not (1) do generate full treept */
)
{
     char    yo[] = "irb_fn";
     int     proc,nprocs;        /* my proc id, total # of procs */
     struct irb_dot *dotbuf;     /* local dot arrays */
     struct irb_dot *dotpt;      /* local dot arrays */
     int     keep, outgoing;     /* message exchange counters */
     int     incoming;           /* message exchange counters */
     int     pdotnum;            /* # of dots - decomposition changes it */
     int     pdottop;            /* dots >= this index are new */
     int    *dotmark = NULL;     /* which side of median for each dot */
     int     dotnum;             /* number of dots */
     int     dotmax = 0;         /* max # of dots arrays can hold */
     int     dottop;             /* dots >= this index are new */
     int     dotnew;             /* # of new dots after send/recv */
     int     proclower;          /* lower proc in partition */
     int     procmid;            /* 1st proc in upper half of part */
     int     set;                /* which part processor is in = 0/1 */
     int     old_set;            /* part processor was in last cut = 0/1 */
     int     root;               /* processor that stores last cut */
     int     num_procs;          /* number of procs in current part */
     int     ierr;               /* error flag. */
     double *value = NULL;       /* temp array for median_find */
     double *wgts = NULL;        /* temp array for median_find */
     double  valuehalf;          /* median cut position */
     double  fractionlo;         /* desired wt in lower half */
     double  cm[3];              /* Center of mass of objects */
     double  evec[3];            /* Eigenvector defining direction */
     int     first_guess = 0;    /* flag if first guess for median search */
     int     allocflag;          /* have to re-allocate space */
     double  time1,time2;        /* timers */
     double  time3,time4;        /* timers */
     double  timestart,timestop; /* timers */
     double  timers[4];          /* diagnostic timers
                                    0 = start-up time before recursion
                                    1 = time before median iterations
                                    2 = time in median iterations
                                    3 = communication time */
     int     counters[7];        /* diagnostic counts
                                    0 = # of median iterations
                                    1 = # of dots sent
                                    2 = # of dots received
                                    3 = most dots this proc ever owns
                                    4 = most dot memory this proc ever allocs
                                    5 = # of times a previous cut is re-used
                                    6 = # of reallocs of dot array */
     int     i,j,k;              /* local variables */
     struct Comm_Obj *cobj=NULL; /* pointer for communication object */
     int     message_tag;        /* message tag */
     int    *proc_list = NULL;   /* list of processors to send dots to */

     IRB_STRUCT *irb = NULL;     /* Pointer to data structures for IRB */
     struct irb_tree *treept = NULL; /* tree of cuts - single cut on exit*/

     double start_time, end_time;
     double lb_time[2];

     /* MPI data types and user functions */

     MPI_Comm local_comm, tmp_comm;

     LB_TRACE_ENTER(lb, yo);
     if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
        MPI_Barrier(lb->Communicator);
        timestart = time1 = LB_Time();
     }

     /* setup for parallel */

     proc = lb->Proc;
     nprocs = lb->Num_Proc;

     /*
      *  Build the IRB Data structure and
      *  set pointers to information in it.
      */

     start_time = LB_Time();
     ierr = LB_IRB_Build_Structure(lb, &pdotnum, &dotmax, wgtflag);
     if (ierr == LB_FATAL || ierr == LB_MEMERR) {
        fprintf(stderr, "[%d] IRB Error in %s:  Error returned from "
                        "LB_IRB_Build_Structure\n", proc, yo);
        LB_TRACE_EXIT(lb, yo);
        return(ierr);
     }

     irb = (IRB_STRUCT *) (lb->Data_Structure);

     dotpt  = irb->Dots;
     treept = irb->Tree_Ptr;
     end_time = LB_Time();
     lb_time[0] = end_time - start_time;
     start_time = end_time;

     /* local copies of calling parameters */

     dottop = dotnum = pdotnum;

     /* initialize timers and counters */

     if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
        counters[0] = 0;
        counters[1] = 0;
        counters[2] = 0;
        counters[3] = dotnum;
        counters[4] = dotmax;
        counters[5] = 0;
        counters[6] = 0;
        timers[1] = 0.0;
        timers[2] = 0.0;
        timers[3] = 0.0;
     }

     /* create mark and list arrays for dots */

     allocflag = 0;
     if (dotmax > 0) {
        dotmark = (int *) LB_MALLOC(dotmax*sizeof(int));
        if (dotmark == NULL) {
           LB_TRACE_EXIT(lb, yo);
           return LB_MEMERR;
        }
        value = (double *) LB_MALLOC(dotmax*sizeof(double));
        if (value == NULL) {
           LB_FREE(&dotmark);
           LB_TRACE_EXIT(lb, yo);
           return LB_MEMERR;
        }
        wgts = (double *) LB_MALLOC(dotmax*sizeof(double));
        if (wgts == NULL) {
           LB_FREE(&dotmark);
           LB_FREE(&value);
           LB_TRACE_EXIT(lb, yo);
           return LB_MEMERR;
        }
     }
     else {
        dotmark = NULL;
        value = NULL;
        wgts = NULL;
     }

     /* set dot weights = 1.0 if user didn't */

     if (!wgtflag)
        for (i = 0; i < dotnum; i++) dotpt[i].Weight = 1.0;

     /* check that all weights > 0 */

     if (check) {
        for (j = i = 0; i < dotnum; i++) if (dotpt[i].Weight == 0.0) j++;
        MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, lb->Communicator);
        if (k > 0 && proc == 0)
           fprintf(stderr, "IRB WARNING: %d dot weights are equal to 0\n", k);

        for (j = i = 0; i < dotnum; i++) if (dotpt[i].Weight < 0.0) j++;
        MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, lb->Communicator);
        if (k > 0) {
           if (proc == 0)
              fprintf(stderr, "IRB ERROR: %d dot weights are < 0\n", k);
           LB_TRACE_EXIT(lb, yo);
           return LB_FATAL;
        }
     }

     /* create local communicator for use in recursion */

     MPI_Comm_dup(lb->Communicator,&local_comm);

     if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
        time2 = LB_Time();
        timers[0] = time2 - time1;
     }

     /* recursively halve until just one proc in partition */

     num_procs = nprocs;
     root = 0;
     old_set = 1;
     treept[proc].parent = 0;
     treept[proc].left_leaf = 0;

     while (num_procs > 1) {

        if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) time1 = LB_Time();

        ierr = LB_divide_machine(lb, proc, local_comm, &set, &proclower,
                                 &procmid, &num_procs, &fractionlo);
        if (ierr != LB_OK && ierr != LB_WARN) {
           LB_FREE(&dotmark);
           LB_FREE(&value);
           LB_FREE(&wgts);
           LB_TRACE_EXIT(lb, yo);
           return (ierr);
        }

        /* create mark array and active list for dots */

        if (allocflag) {
           allocflag = 0;
           LB_FREE(&dotmark);
           LB_FREE(&value);
           LB_FREE(&wgts);
           dotmark = (int *) LB_MALLOC(dotmax*sizeof(int));
           if (dotmark == NULL) {
              LB_TRACE_EXIT(lb, yo);
              return LB_MEMERR;
           }
           value = (double *) LB_MALLOC(dotmax*sizeof(double));
           if (value == NULL) {
              LB_FREE(&dotmark);
              LB_TRACE_EXIT(lb, yo);
              return LB_MEMERR;
           }
           wgts = (double *) LB_MALLOC(dotmax*sizeof(double));
           if (wgts == NULL) {
              LB_FREE(&dotmark);
              LB_FREE(&value);
              LB_TRACE_EXIT(lb, yo);
              return LB_MEMERR;
           }
        }

        for (i = 0; i < dotnum; i++) {
          wgts[i] = dotpt[i].Weight;
        }
        switch (irb->Num_Geom) {
           case 3:
              ierr = LB_inertial3d(dotpt, dotnum, wgtflag, cm, evec, value,
                                   local_comm);
              break;
           case 2:
              ierr = LB_inertial2d(dotpt, dotnum, wgtflag, cm, evec, value,
                                   local_comm);
              break;
           case 1:
              ierr = LB_inertial1d(dotpt, dotnum, wgtflag, cm, evec, value);
              break;
        }

        if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) time2 = LB_Time();

        if (!LB_find_median(value, wgts, dotmark, dotnum, proc, fractionlo,
                            local_comm, &valuehalf, first_guess,
                            &(counters[0]))) {
           fprintf(stderr, "[%d] %s:IRB Error returned from find_median\n",
                   proc, yo);
           LB_FREE(&dotmark);
           LB_FREE(&value);
           LB_FREE(&wgts);
           LB_TRACE_EXIT(lb, yo);
           return LB_FATAL;
        }

        if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) time3 = LB_Time();

        /* store cut info in tree only if I am procmid, the lowest numbered
           processor in right set.  The left set will have proclower which
           if the lowest numbered processor in either set. */

        if (proc == procmid) {
           treept[proc].cm[0] = cm[0];
           treept[proc].cm[1] = cm[1];
           treept[proc].cm[2] = cm[2];
           treept[proc].ev[0] = evec[0];
           treept[proc].ev[1] = evec[1];
           treept[proc].ev[2] = evec[2];
           treept[proc].cut = valuehalf;
           treept[proc].parent = old_set ? -(root+1) : root+1;
           /* The following two will get overwritten when the information
              is assembled if this is not a terminal cut */
           treept[proc].left_leaf = -proclower;
           treept[proc].right_leaf = -procmid;
        }
        old_set = set;
        root = procmid;

        /* outgoing = number of dots to ship to partner */
        /* dottop = number of dots that have never migrated */

        for (i = 0, keep = 0, outgoing = 0; i < dotnum; i++)
           if (dotmark[i] != set)
              outgoing++;
           else if (i < dottop)
              keep++;
        dottop = keep;

        if (outgoing)
           if ((proc_list = (int *) LB_MALLOC(outgoing*sizeof(int))) == NULL) {
              LB_FREE(&dotmark);
              LB_FREE(&value);
              LB_FREE(&wgts);
              LB_TRACE_EXIT(lb, yo);
              return LB_MEMERR;
           }

        ierr = LB_Create_Proc_List(set, dotnum, outgoing, proc_list,
                                   local_comm);
        if (ierr != LB_OK && ierr != LB_WARN) {
           LB_FREE(&proc_list);
           LB_FREE(&dotmark);
           LB_FREE(&value);
           LB_FREE(&wgts);
           LB_TRACE_EXIT(lb, yo);
           return (ierr);
        }

        incoming = 0;
        message_tag = 1;
        ierr = LB_Comm_Create(&cobj, outgoing, proc_list, local_comm,
                              message_tag, lb->Deterministic, &incoming);
        if (ierr != LB_OK && ierr != LB_WARN) {
           LB_FREE(&proc_list);
           LB_FREE(&dotmark);
           LB_FREE(&value);
           LB_FREE(&wgts);
           LB_TRACE_EXIT(lb, yo);
           return (ierr);
        }

        if (outgoing) LB_FREE(&proc_list);

        /* check if need to malloc more space */

        dotnew = dotnum - outgoing + incoming;

        if (dotnew > dotmax) {
           allocflag = 1;
           dotmax = (int) (overalloc * dotnew);
           if (dotmax < dotnew) dotmax = dotnew;
           dotpt = (struct irb_dot *)
              LB_REALLOC(dotpt,(unsigned) dotmax * sizeof(struct irb_dot));
           if (dotpt == NULL) {
              LB_FREE(&dotmark);
              LB_FREE(&value);
              LB_FREE(&wgts);
              LB_TRACE_EXIT(lb, yo);
              return LB_MEMERR;
           }
           irb->Dots = dotpt;
           if (stats) counters[6]++;
        }

        if (stats) {
           counters[1] += outgoing;
           counters[2] += incoming;
           if (dotnew > counters[3]) counters[3] = dotnew;
           if (dotmax > counters[4]) counters[4] = dotmax;
        }

        /* malloc comm send buffer */

        if (outgoing > 0) {
           dotbuf = (struct irb_dot *) LB_MALLOC(outgoing*sizeof(struct
                                                 irb_dot));
           if (dotbuf == NULL) {
              LB_FREE(&dotmark);
              LB_FREE(&value);
              LB_FREE(&wgts);
              LB_TRACE_EXIT(lb, yo);
              return LB_MEMERR;
           }
        }
        else
          dotbuf = NULL;

        /* fill buffer with dots that are marked for sending */
        /* pack down the unmarked ones */

        keep = outgoing = 0;
        for (i = 0; i < dotnum; i++)
           if (dotmark[i] != set)
              memcpy((char *) &dotbuf[outgoing++], (char *) &dotpt[i],
                     sizeof(struct irb_dot));
           else
              memcpy((char *) &dotpt[keep++], (char *) &dotpt[i],
                     sizeof(struct irb_dot));

        ierr = LB_Comm_Do(cobj, message_tag, (char *) dotbuf,
                          sizeof(struct irb_dot), (char *) (&dotpt[keep]));
        if (ierr != LB_OK && ierr != LB_WARN) {
           LB_FREE(&dotmark);
           LB_FREE(&value);
           LB_FREE(&wgts);
           LB_TRACE_EXIT(lb, yo);
           return (ierr);
        }

        ierr = LB_Comm_Destroy(&cobj);

        LB_FREE(&dotbuf);

        dotnum = dotnew;

        /* create new communicators */

        MPI_Comm_split(local_comm,set,proc,&tmp_comm);
        MPI_Comm_free(&local_comm);
        local_comm = tmp_comm;

        if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
           time4 = LB_Time();
           timers[1] += time2 - time1;
           timers[2] += time3 - time2;
           timers[3] += time4 - time3;
        }

     }

     /* have recursed all the way to final single sub-domain */

     /* free all memory used by IRB and MPI */

     MPI_Comm_free(&local_comm);

     LB_FREE(&value);
     LB_FREE(&wgts);
     LB_FREE(&dotmark);

     end_time = LB_Time();
     lb_time[1] = end_time - start_time;

     if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
        MPI_Barrier(lb->Communicator);
        timestop = time1 = LB_Time();
     }

     /* error checking and statistics */

     if (check) IRB_check(lb, dotpt, dotnum, pdotnum);
     if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME))
        IRB_stats(lb, timestop-timestart, dotpt, dotnum, timers, counters,
                  stats);

     /* update calling routine parameters */

     start_time = LB_Time();

     pdotnum = dotnum;
     pdottop = dottop;

     /*  build return arguments */

     *num_import = dotnum - dottop;
     if (*num_import > 0) {
        if (!LB_Special_Malloc(lb,(void **)import_global_ids,*num_import,
                               LB_SPECIAL_MALLOC_GID)) {
           LB_TRACE_EXIT(lb, yo);
           return LB_MEMERR;
        }
        if (!LB_Special_Malloc(lb,(void **)import_local_ids,*num_import,
                               LB_SPECIAL_MALLOC_LID)) {
           LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID);
           LB_TRACE_EXIT(lb, yo);
           return LB_MEMERR;
        }
        if (!LB_Special_Malloc(lb,(void **)import_procs,*num_import,
                               LB_SPECIAL_MALLOC_INT)) {
           LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID);
           LB_Special_Free(lb,(void **)import_local_ids,LB_SPECIAL_MALLOC_LID);
           LB_TRACE_EXIT(lb, yo);
           return LB_MEMERR;
        }

        for (i = 0; i < *num_import; i++) {
           j = i + dottop;
           LB_SET_GID((*import_global_ids)[i], dotpt[j].Tag.Global_ID);
           LB_SET_LID((*import_local_ids)[i], dotpt[j].Tag.Local_ID);
           (*import_procs)[i] = dotpt[j].Tag.Proc;
        }
     }

     if (gen_tree) {
        MPI_Allgather(&treept[proc], sizeof(struct irb_tree), MPI_BYTE,
                      treept, sizeof(struct irb_tree), MPI_BYTE,
                      lb->Communicator);
        for (i = 1; i < nprocs; i++)
           if (treept[i].parent > 0)
              treept[treept[i].parent - 1].left_leaf = i;
           else if (treept[i].parent < 0)
              treept[-treept[i].parent - 1].right_leaf = i;
     }
     else {
        treept[0].right_leaf = -1;
     }

     end_time = LB_Time();
     lb_time[0] += (end_time - start_time);

     if (lb->Debug_Level >= LB_DEBUG_ATIME) {
        if (lb->Proc == lb->Debug_Proc)
           printf("ZOLTAN IRB Times:  \n");
        LB_Print_Stats(lb, lb_time[0], "ZOLTAN       Build:       ");
        LB_Print_Stats(lb, lb_time[1], "ZOLTAN         IRB:         ");
     }

     if (lb->Debug_Level >= LB_DEBUG_ALL) {
        int kk;
        LB_Print_Sync_Start(lb, TRUE);
        printf("ZOLTAN IRB Proc %d Num_Obj=%d Num_Keep=%d Num_Non_Local=%d\n",
               lb->Proc, pdotnum, pdottop, *num_import);
        printf("  Assigned objects:\n");
        for (kk = 0; kk < pdotnum; kk++) {
           printf("    Obj:  %10d      Orig: %4d\n",
                  irb->Dots[kk].Tag.Global_ID, irb->Dots[kk].Tag.Proc);
        }
        printf("  Non_locals:\n");
        for (kk = 0; kk < *num_import; kk++) {
           printf("    Obj:  %10d      Orig: %4d\n", (*import_global_ids)[kk],
                  (*import_procs)[kk]);
        }
        LB_Print_Sync_End(lb, TRUE);
     }

     LB_TRACE_EXIT(lb, yo);
     /* Temporary return value until error codes are fully implemented */
     return(LB_OK);
}


/* ----------------------------------------------------------------------- */

/* consistency checks on IRB results */

static void IRB_check(LB *lb, struct irb_dot *dotpt, int dotnum, int dotorig)
{
     int i, proc, nprocs, total1, total2;
     double weight, wtmax, wtmin, wtone, tolerance;

     MPI_Comm_rank(lb->Communicator,&proc);
     MPI_Comm_size(lb->Communicator,&nprocs);

     /* check that total # of dots remained the same */

     MPI_Allreduce(&dotorig,&total1,1,MPI_INT,MPI_SUM,lb->Communicator);
     MPI_Allreduce(&dotnum,&total2,1,MPI_INT,MPI_SUM,lb->Communicator);
     if (total1 != total2) {
        if (proc == 0)
           fprintf(stderr, "IRB ERROR: Points before IRB = %d, "
                           "Points after IRB = %d\n", total1,total2);
     }

     /* check that result is load-balanced within log2(P)*max-wt */

     weight = wtone = 0.0;
     for (i = 0; i < dotnum; i++) {
        weight += dotpt[i].Weight;
        if (dotpt[i].Weight > wtone) wtone = dotpt[i].Weight;
     }

     MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
     MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
     MPI_Allreduce(&wtone,&tolerance,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);

     /* i = smallest power-of-2 >= nprocs */
     /* tolerance = largest-single-weight*log2(nprocs) */

     for (i = 0; (nprocs >> i) != 0; i++);
     tolerance = tolerance * i * (1.0 + TINY);

     if (wtmax - wtmin > tolerance) {
        if (proc == 0)
           fprintf(stderr,"IRB ERROR: Load-imbalance > tolerance of %g\n",
                          tolerance);
        MPI_Barrier(lb->Communicator);
        if (weight == wtmin)
           fprintf(stderr, "  Proc %d has weight = %g\n",proc,weight);
        if (weight == wtmax)
           fprintf(stderr, "  Proc %d has weight = %g\n",proc,weight);
     }

     MPI_Barrier(lb->Communicator);
}


/* IRB statistics */

static void IRB_stats(LB *lb, double timetotal, struct irb_dot *dotpt,
                      int dotnum, double *timers, int *counters, int stats)
{
     int i, proc, nprocs, sum, min, max, print_proc;
     double ave, rsum, rmin, rmax;
     double weight, wttot, wtmin, wtmax;

     MPI_Comm_rank(lb->Communicator,&proc);
     MPI_Comm_size(lb->Communicator,&nprocs);
     print_proc = lb->Debug_Proc;

     if (proc == print_proc) printf("IRB total time: %g (secs)\n",timetotal);

     if (stats) {
        if (proc == print_proc) printf("IRB Statistics:\n");

        MPI_Barrier(lb->Communicator);

        /* distribution info */

        for (i = 0, weight = 0.0; i < dotnum; i++) weight += dotpt[i].Weight;
        MPI_Allreduce(&weight,&wttot,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
        MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
        MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);

        if (proc == print_proc) {
           printf(" Total weight of dots = %g\n",wttot);
           printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
                  wttot/nprocs,wtmax,wtmin);
        }

        MPI_Barrier(lb->Communicator);
        if (stats > 1) printf("    Proc %d has weight = %g\n",proc,weight);

        for (i = 0, weight = 0.0; i < dotnum; i++)
           if (dotpt[i].Weight > weight) weight = dotpt[i].Weight;
        MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);

        if (proc == print_proc)
           printf(" Maximum weight of single dot = %g\n",wtmax);

        MPI_Barrier(lb->Communicator);
        if (stats > 1) printf("    Proc %d max weight = %g\n",proc,weight);

        /* counter info */

        MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
        MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
        MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
        ave = ((double) sum)/nprocs;
        if (proc == print_proc)
           printf(" Median iter: ave = %g, min = %d, max = %d\n",ave,min,max);
        MPI_Barrier(lb->Communicator);
        if (stats > 1)
           printf("    Proc %d median count = %d\n",proc,counters[0]);

        MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
        MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
        MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
        ave = ((double) sum)/nprocs;
        if (proc == print_proc)
           printf(" Send count: ave = %g, min = %d, max = %d\n",ave,min,max);
        MPI_Barrier(lb->Communicator);
        if (stats > 1)
           printf("    Proc %d send count = %d\n",proc,counters[1]);

        MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
        MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
        MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
        ave = ((double) sum)/nprocs;
        if (proc == print_proc)
           printf(" Recv count: ave = %g, min = %d, max = %d\n",ave,min,max);
        MPI_Barrier(lb->Communicator);
        if (stats > 1)
           printf("    Proc %d recv count = %d\n",proc,counters[2]);

        MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
        MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
        MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
        ave = ((double) sum)/nprocs;
        if (proc == print_proc)
           printf(" Max dots: ave = %g, min = %d, max = %d\n",ave,min,max);
        MPI_Barrier(lb->Communicator);
        if (stats > 1)
           printf("    Proc %d max dots = %d\n",proc,counters[3]);

        MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
        MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
        MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
        ave = ((double) sum)/nprocs;
        if (proc == print_proc)
           printf(" Max memory: ave = %g, min = %d, max = %d\n",ave,min,max);
        MPI_Barrier(lb->Communicator);
        if (stats > 1)
           printf("    Proc %d max memory = %d\n",proc,counters[4]);

        MPI_Allreduce(&counters[6],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
        MPI_Allreduce(&counters[6],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
        MPI_Allreduce(&counters[6],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
        ave = ((double) sum)/nprocs;
        if (proc == print_proc)
           printf(" # of OverAlloc: ave = %g, min = %d, max = %d\n",
                  ave,min,max);
        MPI_Barrier(lb->Communicator);
        if (stats > 1)
           printf("    Proc %d # of OverAlloc = %d\n",proc,counters[6]);
     }

     /* timer info */

     MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
     MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
     MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
     ave = rsum/nprocs;
     if (proc == print_proc)
        printf(" Start-up time %%: ave = %g, min = %g, max = %g\n",
               ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
     MPI_Barrier(lb->Communicator);
     if (stats > 1)
        printf("    Proc %d start-up time = %g\n",proc,timers[0]);

     MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
     MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
     MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
     ave = rsum/nprocs;
     if (proc == print_proc)
        printf(" Pre-median time %%: ave = %g, min = %g, max = %g\n",
               ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
     MPI_Barrier(lb->Communicator);
     if (stats > 1)
        printf("    Proc %d pre-median time = %g\n",proc,timers[1]);

     MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
     MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
     MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
     ave = rsum/nprocs;
     if (proc == print_proc)
        printf(" Median time %%: ave = %g, min = %g, max = %g\n",
               ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
     MPI_Barrier(lb->Communicator);
     if (stats > 1)
        printf("    Proc %d median time = %g\n",proc,timers[2]);

     MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
     MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
     MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
     ave = rsum/nprocs;
     if (proc == print_proc)
        printf(" Comm time %%: ave = %g, min = %g, max = %g\n",
               ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
     MPI_Barrier(lb->Communicator);
     if (stats > 1)
        printf("    Proc %d comm time = %g\n",proc,timers[3]);

     MPI_Barrier(lb->Communicator);
}
