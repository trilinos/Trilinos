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
#include "rcb_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "ha_const.h"
#include "par_median_const.h"

/* Recursive coordinate bisectioning (RCB) routine
   operates on "dots" as defined in shared_const.h
*/

/* Steve Plimpton, Sandia National Labs, ABQ, NM  87185
   Dept 9221, MS 1111
   (505) 845-7873
   sjplimp@cs.sandia.gov
*/

/* Notes:
   dots are balanced across procs by weight (if used)
   on return, proc owns dotnum "dots" in dense array of max-length dotmax
   all dots will be inside (or on surface of) 3-d box defined by rcbbox
   input weights (if used) are real numbers > 0.0
   can extend "Dot_Struct" data structure in calling program, see shared_const.h
   returned RCB tree only contains one cut on each proc,
     need to do MPI_Allgather if wish to collect it on all procs
*/
#define MYHUGE 1.0e30

/*  RCB_DEFAULT_OUTPUT_LEVEL = 0  No statistics logging */
/*  RCB_DEFAULT_OUTPUT_LEVEL = 1  Log times and counts, print summary */
/*  RCB_DEFAULT_OUTPUT_LEVEL = 2  Log times and counts, print for each proc */
#define RCB_DEFAULT_OUTPUT_LEVEL 0
#define RCB_DEFAULT_OVERALLOC 1.0
#define RCB_DEFAULT_REUSE FALSE

/* function prototypes */

static int rcb_fn(LB *, int *, LB_ID_PTR *, LB_ID_PTR *, int **, double,
                 int, int, int, int, int);


/*  Parameters structure for RCB method.  Used in  */
/*  LB_Set_RCB_Param and LB_rcb.                   */
static PARAM_VARS RCB_params[] = {
                  { "RCB_OVERALLOC", NULL, "DOUBLE" },
                  { "RCB_REUSE", NULL, "INT" },
                  { "CHECK_GEOM", NULL, "INT" },
                  { "RCB_OUTPUT_LEVEL", NULL, "INT" },
                  { "KEEP_CUTS", NULL, "INT" },
                  { NULL, NULL, NULL } };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int LB_Set_RCB_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */

    status = LB_Check_Param(name, val, RCB_params, &result, &index);

    return(status);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int LB_rcb(
  LB *lb,                       /* The load-balancing structure with info for
                                   the RCB balancer.                         */
  int *num_import,              /* Number of non-local objects assigned to this
                                   processor in the new decomposition.       */
  LB_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  LB_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  int **import_procs,           /* Returned value:  array of processor IDs for
                                   processors owning the non-local objects in
                                   this processor's new decomposition.       */
  int *num_export,              /* Not computed, set to -1 */
  LB_ID_PTR *export_global_ids, /* Not computed. */
  LB_ID_PTR *export_local_ids,  /* Not computed. */
  int **export_procs            /* Not computed. */
)
{
    /* Wrapper routine to set parameter values and call the real rcb. */
    double overalloc;         /* amount to overallocate by when realloc
                                 of dot array must be done.     
                                 1.0 = no extra; 1.5 = 50% extra; etc. */
    int reuse;                /* (0) don't use (1) use previous cuts
                                 stored in treept at initial guesses.  */
    int wgtflag;              /* (0) do not (1) do use weights.
                                 Multidimensional weights not supported */
    int check_geom;           /* Check input & output for consistency? */
    int stats;                /* Print timing & count summary? */
    int gen_tree;             /* (0) don't (1) generate whole treept to use
                                 later for point and box drop. */


    LB_Bind_Param(RCB_params, "RCB_OVERALLOC", (void *) &overalloc);
    LB_Bind_Param(RCB_params, "RCB_REUSE", (void *) &reuse);
    LB_Bind_Param(RCB_params, "CHECK_GEOM", (void *) &check_geom);
    LB_Bind_Param(RCB_params, "RCB_OUTPUT_LEVEL", (void *) &stats);
    LB_Bind_Param(RCB_params, "KEEP_CUTS", (void *) &gen_tree);

    overalloc = RCB_DEFAULT_OVERALLOC;
    reuse = RCB_DEFAULT_REUSE;
    check_geom = DEFAULT_CHECK_GEOM;
    stats = RCB_DEFAULT_OUTPUT_LEVEL;
    gen_tree = 0;
    wgtflag = (lb->Obj_Weight_Dim > 0); /* Multidim. weights not accepted */

    LB_Assign_Param_Vals(lb->Params, RCB_params, lb->Debug_Level, lb->Proc,
                         lb->Debug_Proc);

    /* Initializations in case of early exit. */
    *num_import = -1;
    *num_export = -1;  /* We don't compute the export map. */

    return(rcb_fn(lb, num_import, import_global_ids, import_local_ids,
		 import_procs, overalloc, reuse, wgtflag,
                 check_geom, stats, gen_tree));
}

/*---------------------------------------------------------------------------*/

static int rcb_fn(
  LB *lb,                       /* The load-balancing structure with info for
                                   the RCB balancer.                         */
  int *num_import,              /* Number of non-local objects assigned to this
                                   processor in the new decomposition.       */
  LB_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  LB_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  int **import_procs,           /* Returned value:  array of processor IDs for
                                   processors owning the non-local objects in
                                   this processor's new decomposition.       */
  double overalloc,             /* amount to overallocate by when realloc
                                   of dot array must be done.     
                                   1.0 = no extra; 1.5 = 50% extra; etc.     */
  int reuse,                    /* (0) don't use (1) use previous cuts
                                   stored in treept at initial guesses.      */
  int wgtflag,                  /* (0) do not (1) do use weights.
                                   Multidimensional weights not supported    */
  int check_geom,               /* Check input & output for consistency?     */
  int stats,                    /* Print timing & count summary?             */
  int gen_tree                  /* (0) do not (1) do generate full treept    */
)
{
  char    yo[] = "rcb_fn";
  int     proc,nprocs;              /* my proc id, total # of procs */
  LB_ID_PTR gidpt;                  /* pointer to rcb->Global_IDs. */
  LB_ID_PTR lidpt;                  /* pointer to rcb->Local_IDs. */
  struct Dot_Struct *dotpt;         /* pointer to rcb->Dots. */
  struct rcb_box boxtmp;            /* tmp rcb box */
  int     pdotnum;                  /* # of dots - decomposition changes it */
  int     pdottop;                  /* dots >= this index are new */
  int    *dotmark = NULL;           /* which side of median for each dot */
  int     dotnum;                   /* number of dots */
  int     dotmax = 0;               /* max # of dots arrays can hold */
  int     dottop;                   /* dots >= this index are new */
  int     proclower;                /* lower proc in partition */
  int     procmid;                  /* 1st proc in upper half of part */
  int     set;                      /* which part processor is in = 0/1 */
  int     old_set;                  /* part processor was in last cut = 0/1 */
  int     root;                     /* processor that stores last cut */
  int     num_procs;                /* number of procs in current part */
  int     dim;                      /* which of 3 axes median cut is on */
  int     ierr;                     /* error flag. */
  int    *proc_list = NULL;         /* temp array for reusing old cuts */
  int     outgoing;                 /* number of outgoing dots for reuse */
  double *coord = NULL;             /* temp array for median_find */
  double *wgts = NULL;              /* temp array for median_find */
  double  valuehalf;                /* median cut position */
  double  fractionlo;               /* desired wt in lower half */
  int     first_guess;              /* flag if first guess for median search */
  int     allocflag;                /* have to re-allocate space */
  double  time1,time2,time3,time4;  /* timers */
  double  timestart,timestop;       /* timers */
  double  timers[4]={0.,0.,0.,0.};  /* diagnostic timers 
			              0 = start-up time before recursion
				      1 = time before median iterations
				      2 = time in median iterations
				      3 = communication time */
  int     counters[7];              /* diagnostic counts
			              0 = # of median iterations
				      1 = # of dots sent
				      2 = # of dots received
				      3 = most dots this proc ever owns
				      4 = most dot memory this proc ever allocs
				      5 = # of times a previous cut is re-used
				      6 = # of reallocs of dot array */
  int     reuse_count[7];           /* counter (as above) for reuse to record
                                       the number of dots premoved */
  int     i,j;                      /* local variables */

  RCB_STRUCT *rcb = NULL;           /* Pointer to data structures for RCB.  */
  struct rcb_box *rcbbox = NULL;    /* bounding box of final RCB sub-domain */
  struct rcb_tree *treept = NULL;   /* tree of RCB cuts - only single cut on 
                                      exit */

  double start_time, end_time;
  double lb_time[2];
  int tmp_nprocs;                   /* added for Tflops_Special */
  int old_nprocs;                   /* added for Tflops_Special */

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;
  MPI_Op box_op;
  MPI_Datatype box_type;
  MPI_User_function LB_rcb_box_merge;

  LB_TRACE_ENTER(lb, yo);
  if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
    MPI_Barrier(lb->Communicator);
    timestart = time1 = LB_Time(lb->Timer);
  }

  /* setup for parallel */

  proc = lb->Proc;
  nprocs = lb->Num_Proc;

  /*
   *  Build the RCB Data structure and 
   *  set pointers to information in it.
   */

  start_time = LB_Time(lb->Timer);
  ierr = LB_RCB_Build_Structure(lb, &pdotnum, &dotmax, wgtflag);
  if (ierr == LB_FATAL || ierr == LB_MEMERR) {
    LB_PRINT_ERROR(proc, yo, 
                   "Error returned from LB_RCB_Build_Structure.");
    LB_TRACE_EXIT(lb, yo);
    return(ierr);
  }

  rcb = (RCB_STRUCT *) (lb->Data_Structure);

  gidpt  = rcb->Global_IDs;
  lidpt  = rcb->Local_IDs;
  dotpt  = rcb->Dots; 
  rcbbox = rcb->Box;
  treept = rcb->Tree_Ptr;
  end_time = LB_Time(lb->Timer);
  lb_time[0] = end_time - start_time;
  start_time = end_time;

  /* local copies of calling parameters */

  dottop = dotnum = pdotnum;

  /* initialize counters */

  counters[0] = 0;
  counters[1] = 0;
  counters[2] = 0;
  counters[3] = dotnum;
  counters[4] = dotmax;
  counters[5] = 0;
  counters[6] = 0;
  for (i = 0; i < 7; i++) reuse_count[i] = 0;


  /* create mark and list arrays for dots */

  allocflag = 0;
  if (dotmax > 0) {
    dotmark = (int *) LB_MALLOC(dotmax*sizeof(int));
    if (dotmark == NULL) {
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    coord = (double *) LB_MALLOC(dotmax*sizeof(double));
    if (coord == NULL) {
      LB_FREE(&dotmark);
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    wgts = (double *) LB_MALLOC(dotmax*sizeof(double));
    if (wgts == NULL) {
      LB_FREE(&dotmark);
      LB_FREE(&coord);
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
  }
  else {
    dotmark = NULL;
    coord = NULL;
    wgts = NULL;
  }

  /* if reuse is turned on, turn on gen_tree since it is needed. */
  /* Also, if this is not first time through, send dots to previous proc. */
  if (reuse) {
    gen_tree = 1;

    if (treept[0].dim != -1) {
      /* find previous location of dots */
      for (outgoing = i = 0; i < dotnum; i++) {
        ierr = LB_Point_Assign(lb, dotpt[i].X, &dotmark[i]);
        if (dotmark[i] != proc) outgoing++;
      }

      if (outgoing)
        if ((proc_list = (int *) LB_MALLOC(outgoing*sizeof(int))) == NULL) {
          LB_FREE(&dotmark);
          LB_FREE(&coord);
          LB_FREE(&wgts);
          LB_TRACE_EXIT(lb, yo);
          return LB_MEMERR;
        }

      for (dottop = j = i = 0; i < dotnum; i++)
        if (dotmark[i] != proc)
          proc_list[j++] = dotmark[i];
        else
          dottop++;

      /* move dots */
      allocflag = 0;
      ierr = LB_RB_Send_Dots(lb, &gidpt, &lidpt, &dotpt, dotmark, proc_list,
                             outgoing, &dotnum, &dotmax, proc, &allocflag,
                             overalloc, stats, reuse_count, lb->Communicator);
      if (ierr) {
        LB_PRINT_ERROR(proc, yo, "Error returned from LB_RB_Send_Dots.");
        LB_FREE(&proc_list);
        LB_FREE(&dotmark);
        LB_FREE(&coord);
        LB_FREE(&wgts);
        LB_TRACE_EXIT(lb, yo);
        return (ierr);
      }

      if (allocflag) {
        /*
         * gidpt, lidpt and dotpt were reallocated in LB_RB_Send_Dots;
         * store their values in rcb.
         */
        rcb->Global_IDs = gidpt;
        rcb->Local_IDs = lidpt;
        rcb->Dots = dotpt;
      }

      /* update counters */
      if (dotnum > counters[3]) counters[3] = dotnum;
      if (dotmax > counters[4]) counters[4] = dotmax;
      counters[6] += reuse_count[6];

      if (outgoing) LB_FREE(&proc_list);
    }
  }

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);

  MPI_Op_create(&LB_rcb_box_merge,1,&box_op);

  /* set dot weights = 1.0 if user didn't */

  if (!wgtflag)
    for (i = 0; i < dotnum; i++) dotpt[i].Weight = 1.0;

  if (check_geom) {
    ierr = LB_RB_check_geom_input(lb, dotpt, dotnum);
    if (ierr == LB_FATAL) {
      LB_PRINT_ERROR(proc, yo, "Error returned from LB_RB_check_geom_input");
      LB_TRACE_EXIT(lb, yo);
      return(ierr);
    }
  }

  /* initialize sub-domain bounding box to entire domain */

  boxtmp.lo[0] = boxtmp.lo[1] = boxtmp.lo[2] = MYHUGE;
  boxtmp.hi[0] = boxtmp.hi[1] = boxtmp.hi[2] = -MYHUGE;
  
  for (i = 0; i < dotnum; i++) {
    for (j = 0; j < 3; j++) {
      if (dotpt[i].X[j] < boxtmp.lo[j])
	boxtmp.lo[j] = dotpt[i].X[j];
      if (dotpt[i].X[j] > boxtmp.hi[j])
	boxtmp.hi[j] = dotpt[i].X[j];
    }
  }

  MPI_Allreduce(&boxtmp,rcbbox,1,box_type,box_op,lb->Communicator);

  /* create local communicator for use in recursion */

  if (lb->Tflops_Special)
     local_comm = lb->Communicator;
  else
     MPI_Comm_dup(lb->Communicator,&local_comm);

  if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
    time2 = LB_Time(lb->Timer);
    timers[0] = time2 - time1;
  }

  /* recursively halve until just one proc in partition */
  
  old_nprocs = num_procs = nprocs;
  root = 0;
  old_set = 1;
  treept[proc].parent = 0;
  treept[proc].left_leaf = 0;
  if (lb->Tflops_Special) {
     proclower = 0;
     tmp_nprocs = nprocs;
  }

  while (num_procs > 1 || (lb->Tflops_Special && tmp_nprocs > 1)) {
    /* tmp_nprocs is size of largest partition - force all processors to go
       through all levels of rcb */
    if (lb->Tflops_Special) tmp_nprocs = (tmp_nprocs + 1)/2;

    if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) time1 = LB_Time(lb->Timer);

    ierr = LB_divide_machine(lb, proc, local_comm, &set, &proclower,
                             &procmid, &num_procs, &fractionlo);
    if (ierr != LB_OK && ierr != LB_WARN) {
      LB_FREE(&dotmark);
      LB_FREE(&coord);
      LB_FREE(&wgts);
      LB_TRACE_EXIT(lb, yo);
      return (ierr);
    }
    
    /* dim = dimension (xyz = 012) to bisect on */

    dim = 0;
    if (rcbbox->hi[1] - rcbbox->lo[1] >
	rcbbox->hi[0] - rcbbox->lo[0])
      dim = 1;
    if (dim == 0 && rcbbox->hi[2] - rcbbox->lo[2] >
	rcbbox->hi[0] - rcbbox->lo[0])
      dim = 2;
    if (dim == 1 && rcbbox->hi[2] - rcbbox->lo[2] >
	rcbbox->hi[1] - rcbbox->lo[1])
      dim = 2;
    
    /* create mark array and active list for dots */

    if (allocflag) {
      allocflag = 0;
      LB_FREE(&dotmark);
      LB_FREE(&coord);
      LB_FREE(&wgts);
      dotmark = (int *) LB_MALLOC(dotmax*sizeof(int));
      if (dotmark == NULL) {
        LB_TRACE_EXIT(lb, yo);
        return LB_MEMERR;
      }
      coord = (double *) LB_MALLOC(dotmax*sizeof(double));
      if (coord == NULL) {
        LB_FREE(&dotmark);
        LB_TRACE_EXIT(lb, yo);
        return LB_MEMERR;
      }
      wgts = (double *) LB_MALLOC(dotmax*sizeof(double));
      if (wgts == NULL) {
        LB_FREE(&dotmark);
        LB_FREE(&coord);
        LB_TRACE_EXIT(lb, yo);
        return LB_MEMERR;
      }
    }

    /* copy correct coordinate value into the temporary array */
    for (i = 0; i < dotnum; i++) {
      coord[i] = dotpt[i].X[dim];
      wgts[i] = dotpt[i].Weight;
    }

    /* determine if there is a first guess to use */
    if (reuse && dim == treept[procmid].dim) {
      if (stats) counters[5]++;
      valuehalf = treept[procmid].cut;
      first_guess = 1;
    }
    else first_guess = 0;

    if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) time2 = LB_Time(lb->Timer);

    if (!LB_find_median(lb, coord, wgts, dotmark, dotnum, proc, fractionlo,
                        local_comm, &valuehalf, first_guess, &(counters[0]),
                        nprocs, old_nprocs, proclower)) {
      LB_PRINT_ERROR(proc, yo, "Error returned from LB_find_median.");
      LB_FREE(&dotmark);
      LB_FREE(&coord);
      LB_FREE(&wgts);
      LB_TRACE_EXIT(lb, yo);
      return LB_FATAL;
    }

    if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) time3 = LB_Time(lb->Timer);

    /* store cut info in tree only if I am procmid */

    if (proc == procmid) {
      treept[proc].dim = dim;
      treept[proc].cut = valuehalf;
      treept[proc].parent = old_set ? -(root+1) : root+1;
      /* The following two will get overwritten when the information is
         assembled if this is not a terminal cut */
      treept[proc].left_leaf = -proclower;
      treept[proc].right_leaf = -procmid;
    }
    old_set = set;
    root = procmid;
    
    /* use cut to shrink RCB domain bounding box */

    if (old_nprocs > 1) {
      if (!set)
        rcbbox->hi[dim] = valuehalf;
      else
        rcbbox->lo[dim] = valuehalf;
    }

    allocflag = 0;
    ierr = LB_RB_Send_Outgoing(lb, &gidpt, &lidpt, &dotpt, dotmark, &dottop,
                               &dotnum, &dotmax, set, &allocflag, overalloc,
                               stats, counters, local_comm, proclower,
                               old_nprocs);
    if (ierr) {
      LB_PRINT_ERROR(proc, yo, "Error returned from LB_RB_Send_Outgoing.");
      LB_FREE(&dotmark);
      LB_FREE(&coord);
      LB_FREE(&wgts);
      LB_TRACE_EXIT(lb, yo);
      return (ierr);
    }

    if (allocflag) {
      /* 
       * gidpt, lidpt and dotpt were reallocated in LB_RB_Send_Outgoing;
       * store their values in rcb.
       */
      rcb->Global_IDs = gidpt;
      rcb->Local_IDs = lidpt;
      rcb->Dots = dotpt;
    }

    /* create new communicators */

    if (lb->Tflops_Special) {
       if (set) proclower = procmid;
       old_nprocs = num_procs;
    }
    else {
       MPI_Comm_split(local_comm,set,proc,&tmp_comm);
       MPI_Comm_free(&local_comm);
       local_comm = tmp_comm;
    }

    if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
      time4 = LB_Time(lb->Timer);
      timers[1] += time2 - time1;
      timers[2] += time3 - time2;
      timers[3] += time4 - time3;
    }
  }

  /* have recursed all the way to final single sub-domain */

  /* free all memory used by RCB and MPI */

  if (!lb->Tflops_Special) MPI_Comm_free(&local_comm);
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);

  LB_FREE(&coord);
  LB_FREE(&wgts);
  LB_FREE(&dotmark);

  end_time = LB_Time(lb->Timer);
  lb_time[1] = end_time - start_time;

  if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
    MPI_Barrier(lb->Communicator);
    timestop = time1 = LB_Time(lb->Timer);
  }

  /* error checking and statistics */

  if (check_geom) {
    ierr = LB_RB_check_geom_output(lb, dotpt,dotnum,pdotnum,rcbbox);
    if (ierr == LB_FATAL) {
      LB_PRINT_ERROR(proc, yo, "Error returned from LB_RB_check_geom_output");
      LB_TRACE_EXIT(lb, yo);
      return(ierr);
    }
  }

  if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) 
    LB_RB_stats(lb, timestop-timestart,dotpt,dotnum,
                timers,counters,stats,reuse_count,rcbbox,reuse);

  /* update calling routine parameters */
  
  start_time = LB_Time(lb->Timer);

  pdotnum = dotnum;
  pdottop = dottop;

  /*  build return arguments */

  if (lb->Return_Lists) {
    *num_import = dotnum - dottop;
    if (*num_import > 0) {
      ierr = LB_RB_Return_Arguments(lb, gidpt, lidpt, dotpt, *num_import,
                                    import_global_ids, import_local_ids,
                                    import_procs, dottop);
      if (ierr) {
         LB_PRINT_ERROR(proc,yo,"Error returned from LB_RB_Return_Arguments.");
         LB_TRACE_EXIT(lb, yo);
         return ierr;
      }
    }
  }

  if (gen_tree) {
    MPI_Allgather(&treept[proc], sizeof(struct rcb_tree), MPI_BYTE,
       treept, sizeof(struct rcb_tree), MPI_BYTE, lb->Communicator);
    treept[0].dim = 0;
    for (i = 1; i < nprocs; i++)
       if (treept[i].parent > 0)
          treept[treept[i].parent - 1].left_leaf = i;
       else if (treept[i].parent < 0)
          treept[-treept[i].parent - 1].right_leaf = i;
  }
  else {
    treept[0].dim = -1;
  }

  end_time = LB_Time(lb->Timer);
  lb_time[0] += (end_time - start_time);

  if (lb->Debug_Level >= LB_DEBUG_ATIME) {
    if (lb->Proc == lb->Debug_Proc) {
      printf("ZOLTAN RCB Times:  \n");
    }
    LB_Print_Stats(lb->Communicator, lb->Debug_Proc, lb_time[0], 
                   "ZOLTAN     Build:       ");
    LB_Print_Stats(lb->Communicator, lb->Debug_Proc, lb_time[1], 
                   "ZOLTAN     RCB:         ");
  }

  if (lb->Debug_Level >= LB_DEBUG_ALL) {
    LB_RB_Print_All(lb, rcb->Global_IDs, rcb->Dots, 
                    pdotnum, pdottop, *num_import,
                    *import_global_ids, *import_procs);
  }

  LB_TRACE_EXIT(lb, yo);
  /* Temporary return value until error codes are fully implemented */
  return(LB_OK);  
}



/* ----------------------------------------------------------------------- */

/* MPI user-defined reduce operations */

/* min/max merge of each component of a rcb_box */

void LB_rcb_box_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

{
  int i;
  struct rcb_box *box1,*box2;

  box1 = (struct rcb_box *) in;
  box2 = (struct rcb_box *) inout;

  for (i = 0; i < 3; i++) {
    if (box1->lo[i] < box2->lo[i])
      box2->lo[i] = box1->lo[i];
    if (box1->hi[i] > box2->hi[i])
      box2->hi[i] = box1->hi[i];
  }
}
