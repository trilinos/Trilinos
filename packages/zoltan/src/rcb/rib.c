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
 *    Revision: 1.6.2.2 $
 ****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "rib_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "ha_const.h"
#include "par_median_const.h"

/* Inertial recursive bisection (RIB) load balancing routine operates on
   "dots" as defined in shared_const.h */

/* Notes:
   dots are balanced across procs by weight (if used)
   on return, proc owns dotnum "dots" in dense array of max-length dotmax
   input weights (if used) are real numbers > 0.0
   can extend "Dot_Struct" data structure in calling program, see shared_const.h
   returned tree only contains one cut on each proc,
   need to do MPI_Allgather if wish to collect it on all procs */
#define MYHUGE 1.0e30

/*  RIB_OUTPUT_LEVEL = 0  No statistics logging */
/*  RIB_OUTPUT_LEVEL = 1  Log times and counts, print summary */
/*  RIB_OUTPUT_LEVEL = 2  Log times and counts, print for each proc */
#define RIB_DEFAULT_OUTPUT_LEVEL 0
#define RIB_DEFAULT_OVERALLOC 1.0


static int rib_fn(LB *, int *, LB_ID_PTR *, LB_ID_PTR *, int **, double,
            int, int, int, int);
static void print_rib_tree(LB *, struct rib_tree *);
/* for Tflops_Special */
static void LB_min_max(double *, double *, int, int, int, MPI_Comm);

/*  Parameters structure for RIB method.  Used in  */
/*  LB_Set_RIB_Param and LB_RIB.                      */
static PARAM_VARS RIB_params[] = {
               { "RIB_OVERALLOC", NULL, "DOUBLE" },
               { "CHECK_GEOM", NULL, "INT" },
               { "RIB_OUTPUT_LEVEL", NULL, "INT" },
               { "KEEP_CUTS", NULL, "INT" },
               { NULL, NULL, NULL } };

/*---------------------------------------------------------------------------*/

int LB_Set_RIB_Param(
  char *name,                 /* name of variable */
  char *val                   /* value of variable */
)
{
  int status;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */

  status = LB_Check_Param(name, val, RIB_params, &result, &index);

  return(status);
}

/*---------------------------------------------------------------------------*/

int LB_rib(
  LB *lb,                       /* The load-balancing structure with info for
                                the RIB balancer.                       */
  int *num_import,              /* Number of non-local objects assigned to
                                this processor in the new decomposition.*/
  LB_ID_PTR *import_global_ids, /* Returned value: array of global IDs for
                                non-local objects in this processor's new
                                decomposition.                          */
  LB_ID_PTR *import_local_ids,  /* Returned value: array of local IDs for
                                non-local objects in this processor's new
                                decomposition.                          */
  int **import_procs,           /* Returned value: array of processor IDs for
                                processors owning the non-local objects in
                                this processor's new decomposition.     */
  int *num_export,              /* Not computed, set to -1 */
  LB_ID_PTR *export_global_ids, /* Not computed. */
  LB_ID_PTR *export_local_ids,  /* Not computed. */
  int **export_procs            /* Not computed. */
)
{
  /* Wrapper routine to set parameter values and call the real rib. */
  double overalloc;           /* amount to overallocate by when realloc
                              of dot array must be done.
                              1.0 = no extra; 1.5 = 50% extra; etc. */
  int wgtflag;                /* (0) do not (1) do use weights.
                              Multidimensional weights not supported */
  int check_geom;             /* Check input & output for consistency? */
  int stats;                  /* Print timing & count summary? */
  int gen_tree;               /* (0) don't (1) generate whole treept to use
                              later for point and box drop. */

  LB_Bind_Param(RIB_params, "RIB_OVERALLOC", (void *) &overalloc);
  LB_Bind_Param(RIB_params, "CHECK_GEOM", (void *) &check_geom);
  LB_Bind_Param(RIB_params, "RIB_OUTPUT_LEVEL", (void *) &stats);
  LB_Bind_Param(RIB_params, "KEEP_CUTS", (void *) &gen_tree);

  overalloc = RIB_DEFAULT_OVERALLOC;
  check_geom = DEFAULT_CHECK_GEOM;
  stats = RIB_DEFAULT_OUTPUT_LEVEL;
  gen_tree = 0;
  wgtflag = (lb->Obj_Weight_Dim > 0); /* Multidim. weights not accepted */

  LB_Assign_Param_Vals(lb->Params, RIB_params, lb->Debug_Level, lb->Proc,
                    lb->Debug_Proc);

  /* Initializations in case of early exit. */
  *num_import = -1;
  *num_export = -1;  /* We don't compute the export map. */

  return(rib_fn(lb, num_import, import_global_ids, import_local_ids,
                import_procs, overalloc, wgtflag, check_geom, stats, gen_tree));
}

/*---------------------------------------------------------------------------*/

static int rib_fn(
  LB *lb,                       /* The load-balancing structure with info for
                                the RIB balancer. */
  int *num_import,              /* Number of non-local objects assigned to
                                this processor in the new decomposition.*/
  LB_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                non-local objects in this processor's new
                                decomposition. */
  LB_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                non-local objects in this processor's new
                                decomposition. */
  int **import_procs,           /* Returned value: array of processor IDs for
                                processors owning the non-local objects in
                                this processor's new decomposition. */
  double overalloc,             /* amount to overallocate by when realloc
                                of dot array must be done.
                                  1.0 = no extra; 1.5 = 50% extra; etc. */
  int wgtflag,                  /* (0) do not (1) do use weights.
                                Multidimensional weights not supported */
  int check_geom,               /* Check input & output for consistency? */
  int stats,                    /* Print timing & count summary? */
  int gen_tree                  /* (0) do not (1) do generate full treept */
)
{
  char    yo[] = "rib_fn";
  int     proc,nprocs;        /* my proc id, total # of procs */
  LB_ID_PTR gidpt = NULL;     /* local global IDs array. */
  LB_ID_PTR lidpt = NULL;     /* local local IDs array. */
  struct Dot_Struct *dotpt;   /* local dot arrays */
  int     pdotnum;            /* # of dots - decomposition changes it */
  int     pdottop;            /* dots >= this index are new */
  int    *dotmark = NULL;     /* which side of median for each dot */
  int     dotnum;             /* number of dots */
  int     dotmax = 0;         /* max # of dots arrays can hold */
  int     dottop;             /* dots >= this index are new */
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
  double  timers[4]={0.,0.,0.,0.}; 
                              /* diagnostic timers
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
  int     i;                  /* local variables */
  int     use_ids;            /* When true, global and local IDs will be
                                 stored along with dots in the RCB_STRUCT.
                                 When false, storage, manipulation, and
                                 communication of IDs is avoided.     
                                 Set by call to LB_Use_IDs().         */


  RIB_STRUCT *rib = NULL;     /* Pointer to data structures for RIB */
  struct rib_tree *treept = NULL; /* tree of cuts - single cut on exit*/

  double start_time, end_time;
  double lb_time[2];
  int tmp_nprocs;             /* added for Tflops_Special */
  int old_nprocs;             /* added for Tflops_Special */
  double valuelo;             /* smallest value of value[i] */
  double valuehi;             /* largest value of value[i] */
  double weight;              /* weight for current partition */
  double weightlo;            /* weight of lower side of cut */
  double weighthi;            /* weight of upper side of cut */
  int *dotlist;               /* list of dots for find_median */

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;

  LB_TRACE_ENTER(lb, yo);
  if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
    MPI_Barrier(lb->Communicator);
    timestart = time1 = LB_Time(lb->Timer);
  }

  /* setup for parallel */

  proc = lb->Proc;
  nprocs = lb->Num_Proc;

  /*
   * Determine whether to store, manipulate, and communicate global and
   * local IDs.
   */
  use_ids = LB_Use_IDs(lb);

  /*
   *  Build the RIB Data structure and
   *  set pointers to information in it.
   */

  start_time = LB_Time(lb->Timer);
  ierr = LB_RIB_Build_Structure(lb, &pdotnum, &dotmax, wgtflag, use_ids);
  if (ierr == LB_FATAL || ierr == LB_MEMERR) {
    LB_PRINT_ERROR(proc, yo, "Error returned from LB_RIB_Build_Structure.");
    LB_TRACE_EXIT(lb, yo);
    return(ierr);
  }

  rib = (RIB_STRUCT *) (lb->Data_Structure);

  gidpt = rib->Global_IDs;
  lidpt = rib->Local_IDs;
  dotpt  = rib->Dots;
  treept = rib->Tree_Ptr;
  end_time = LB_Time(lb->Timer);
  lb_time[0] = end_time - start_time;
  start_time = end_time;

  /* local copies of calling parameters */

  dottop = dotnum = pdotnum;

  /* initialize timers and counters */

  counters[0] = 0;
  counters[1] = 0;
  counters[2] = 0;
  counters[3] = dotnum;
  counters[4] = dotmax;
  counters[5] = 0;
  counters[6] = 0;

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
    dotlist = (int *) LB_MALLOC(dotmax*sizeof(int));
    if (dotlist == NULL) {
      LB_FREE(&wgts);
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
    dotlist = NULL;
  }

  /* set dot weights = 1.0 if user didn't and determine total weight */

  if (!wgtflag) {
    for (i = 0; i < dotnum; i++) dotpt[i].Weight = 1.0;
    weightlo = (double) dotnum;
  }
  else
    for (weightlo = 0.0, i=0; i < dotnum; i++)
       weightlo += dotpt[i].Weight;
  MPI_Allreduce(&weightlo, &weight, 1, MPI_DOUBLE, MPI_SUM, lb->Communicator);

  if (check_geom) {
    ierr = LB_RB_check_geom_input(lb, dotpt, dotnum);
    if (ierr == LB_FATAL) {
      LB_PRINT_ERROR(proc, yo, "Error returned from LB_RB_check_geom_input");
      LB_TRACE_EXIT(lb, yo);
      return(ierr);
    }
  }

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
       through all levels of rib */
    if (lb->Tflops_Special) tmp_nprocs = (tmp_nprocs + 1)/2;

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
      LB_FREE(&dotlist);
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
      dotlist = (int *) LB_MALLOC(dotmax*sizeof(int));
      if (dotlist == NULL) {
        LB_FREE(&wgts);
        LB_FREE(&dotmark);
        LB_FREE(&value);
        LB_TRACE_EXIT(lb, yo);
        return LB_MEMERR;
      }
    }

    for (i = 0; i < dotnum; i++) {
      wgts[i] = dotpt[i].Weight;
    }
    if (old_nprocs > 1) {    /* for Tflops_Special */
       switch (rib->Num_Geom) {
          case 3:
             ierr = LB_inertial3d(lb, dotpt, dotnum, wgtflag, cm, evec, value,
                                  local_comm, proc, old_nprocs, proclower);
             break;
          case 2:
             ierr = LB_inertial2d(lb, dotpt, dotnum, wgtflag, cm, evec, value,
                                  local_comm, proc, old_nprocs, proclower);
             break;
          case 1:
             ierr = LB_inertial1d(dotpt, dotnum, wgtflag, cm, evec, value);
             break;
       }
       valuelo = MYHUGE;
       valuehi = -MYHUGE;
       for (i = 0; i < dotnum; i++) {
          if (value[i] < valuelo) valuelo = value[i];
          if (value[i] > valuehi) valuehi = value[i];
       }
       if (lb->Tflops_Special)
          LB_min_max(&valuelo, &valuehi, proclower, proc, old_nprocs,
                     local_comm);
       else {
          valuehalf = valuehi;
          MPI_Allreduce(&valuehalf,&valuehi,1,MPI_DOUBLE,MPI_MAX,local_comm);
          valuehalf = valuelo;
          MPI_Allreduce(&valuehalf,&valuelo,1,MPI_DOUBLE,MPI_MIN,local_comm);
       }
    }
    else {                    /* For Tflops_Special initialize value */
       for (i = 0; i < dotmax; i++)
          value[i] = 0.0;
       valuelo = valuehi = 0.0;
    }

    if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) 
      time2 = LB_Time(lb->Timer);

    if (!LB_find_median(lb, value, wgts, dotmark, dotnum, proc, fractionlo,
                        local_comm, &valuehalf, first_guess,
                        &(counters[0]), nprocs, old_nprocs, proclower,
                        wgtflag, valuelo, valuehi, weight, &weightlo,
                        &weighthi, dotlist)) {
      LB_PRINT_ERROR(proc, yo, "Error returned from find_median.");
      LB_FREE(&dotmark);
      LB_FREE(&value);
      LB_FREE(&wgts);
      LB_TRACE_EXIT(lb, yo);
      return LB_FATAL;
    }

    if (set)    /* set weight for current partition */
       weight = weighthi;
    else
       weight = weightlo;

    if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) 
      time3 = LB_Time(lb->Timer);

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

    ierr = LB_RB_Send_Outgoing(lb, &gidpt, &lidpt, &dotpt, dotmark, &dottop,
                               &dotnum, &dotmax, set, &allocflag, overalloc,
                               stats, counters, use_ids, local_comm, proclower,
                               old_nprocs);
    if (ierr) {
      LB_PRINT_ERROR(proc, yo, "Error returned from LB_RB_Send_Outgoing.");
      LB_FREE(&dotmark);
      LB_FREE(&value);
      LB_FREE(&wgts);
      LB_TRACE_EXIT(lb, yo);
      return ierr;
    }
    
    if (allocflag) {
      /* 
       * gidpt, lidpt and dotpt were reallocated in LB_RB_Send_Outgoing;
       * store their values in rib.
       */
      rib->Global_IDs = gidpt;
      rib->Local_IDs = lidpt;
      rib->Dots = dotpt;
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

  /* free all memory used by RIB and MPI */

  if (!lb->Tflops_Special) MPI_Comm_free(&local_comm);

  LB_FREE(&value);
  LB_FREE(&wgts);
  LB_FREE(&dotmark);
  LB_FREE(&dotlist);

  end_time = LB_Time(lb->Timer);
  lb_time[1] = end_time - start_time;

  if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME)) {
    MPI_Barrier(lb->Communicator);
    timestop = time1 = LB_Time(lb->Timer);
  }

  /* error checking and statistics */

  if (check_geom) {
    ierr = LB_RB_check_geom_output(lb, dotpt,dotnum,pdotnum,NULL);
    if (ierr == LB_FATAL) {
      LB_PRINT_ERROR(proc, yo, "Error returned from LB_RB_check_geom_output");
      LB_TRACE_EXIT(lb, yo);
      return(ierr);
    }
  }

  if (stats || (lb->Debug_Level >= LB_DEBUG_ATIME))
    LB_RB_stats(lb, timestop-timestart, dotpt, dotnum, timers, counters,
                stats, NULL, NULL, FALSE);

  /* update calling routine parameters */

  start_time = LB_Time(lb->Timer);

  pdotnum = dotnum;
  pdottop = dottop;

  /*  build return arguments */

  if (lb->Return_Lists) {
    /* lb->Return_Lists is true ==> use_ids is true */
    *num_import = dotnum - dottop;
    if (*num_import > 0) {
      ierr = LB_RB_Return_Arguments(lb, gidpt, lidpt, dotpt, *num_import,
                                    import_global_ids, import_local_ids,
                                    import_procs, dottop);
      if (ierr) {
        LB_PRINT_ERROR(proc, yo,
                       "Error returned from LB_RB_Return_Arguments.");
        LB_TRACE_EXIT(lb, yo);
        return ierr;
      }
    }
  }

  if (lb->Debug_Level >= LB_DEBUG_ALL)
    print_rib_tree(lb, &(treept[proc]));

  if (gen_tree) {
    MPI_Allgather(&treept[proc], sizeof(struct rib_tree), MPI_BYTE,
                  treept, sizeof(struct rib_tree), MPI_BYTE,
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

  end_time = LB_Time(lb->Timer);
  lb_time[0] += (end_time - start_time);

  if (lb->Debug_Level >= LB_DEBUG_ATIME) {
    if (lb->Proc == lb->Debug_Proc)
      printf("ZOLTAN RIB Times:  \n");
    LB_Print_Stats(lb->Communicator, lb->Debug_Proc, lb_time[0], 
                   "ZOLTAN       Build:       ");
    LB_Print_Stats(lb->Communicator, lb->Debug_Proc, lb_time[1], 
                   "ZOLTAN         RIB:         ");
  }

  if (lb->Debug_Level >= LB_DEBUG_ALL) {
    /* lb->Debug_Level >= LB_DEBUG_ALL ==> use_ids is true */
    LB_RB_Print_All(lb, rib->Global_IDs, rib->Dots, 
                    pdotnum, pdottop, *num_import, 
                    *import_global_ids, *import_procs);
  }

  /* Free memory allocated by the algorithm. */
  if (!gen_tree) {
    /* Free all memory used. */
    LB_RIB_Free_Structure(lb);
  }
  else {
    /* Free only Dots and IDs; keep other structures. */
    LB_FREE(&(rib->Global_IDs));
    LB_FREE(&(rib->Local_IDs));
    LB_FREE(&(rib->Dots));
  }

  LB_TRACE_EXIT(lb, yo);
  /* Temporary return value until error codes are fully implemented */
  return(LB_OK);
}

static void print_rib_tree(LB *lb, struct rib_tree *treept)
{
  LB_Print_Sync_Start(lb->Communicator, TRUE);
  printf("Proc %d:  Tree Struct:\n", lb->Proc);
  printf("          cm         = (%e,%e,%e)\n", 
         treept->cm[0], treept->cm[1], treept->cm[2]);
  printf("          ev         = (%e,%e,%e)\n", 
         treept->ev[0], treept->ev[1], treept->ev[2]);
  printf("          cut        = %e\n", treept->cut);
  printf("          parent     = %d\n", treept->parent);
  printf("          left_leaf  = %d\n", treept->left_leaf);
  printf("          right_leaf = %d\n", treept->right_leaf);
  LB_Print_Sync_End(lb->Communicator, TRUE);
}

static void LB_min_max(
   double   *min,             /* minimum value */
   double   *max,             /* maximum value */
   int      proclower,        /* smallest processor in partition */
   int      proc,             /* rank of processor in global partition */
   int      nprocs,           /* number of processors in partition */
   MPI_Comm comm
)
{
   double   tmp[2], tmp1[2];  /* temporaries for min/max */
   int      tag = 32100;      /* message tag */
   int      rank;             /* rank of processor in partition */
   int      partner;          /* message partner in binary exchange */
   int      to;               /* message partner not in binary exchange */
   int      mask;             /* mask to determine communication partner */
   int      nprocs_small;     /* largest power of 2 contained in nprocs */
   int      hbit;             /* 2^hbit = nproc_small */
   MPI_Status status;

   /* This routine finds the global min of min and the global max of max */

   rank = proc - proclower;
   tmp[0] = *min;
   tmp[1] = *max;

   /* Find next lower power of 2. */
   for (hbit = 0; (nprocs >> hbit) != 1; hbit++);
 
   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nprocs) {
      nprocs_small *= 2;
      hbit++;
   }
 
   to = proclower + (rank ^ nprocs_small);
   if (rank & nprocs_small) {  /* processors greater than largest power of 2 */
      MPI_Send(tmp, 2, MPI_DOUBLE, to, tag, comm);
      tag += hbit + 1;
      MPI_Recv(tmp, 2, MPI_DOUBLE, to, tag, comm, &status);
   }
   else {   /* processors within greatest power of 2 */
      if (rank + nprocs_small < nprocs) {
         MPI_Recv(tmp1, 2, MPI_DOUBLE, to, tag, comm, &status);
         if (tmp1[0] < tmp[0]) tmp[0] = tmp1[0];
         if (tmp1[1] > tmp[1]) tmp[1] = tmp1[1];
      }  
      for (mask = nprocs_small >> 1; mask; mask >>= 1) { /* binary exchange */
         tag++;
         partner = proclower + (rank ^ mask);
         MPI_Send(tmp, 2, MPI_DOUBLE, partner, tag, comm);
         MPI_Recv(tmp1, 2, MPI_DOUBLE, partner, tag, comm, &status);
         if (tmp1[0] < tmp[0]) tmp[0] = tmp1[0];
         if (tmp1[1] > tmp[1]) tmp[1] = tmp1[1];
      }  
      tag++;
      if (rank + nprocs_small < nprocs)
         MPI_Send(tmp, 2, MPI_DOUBLE, to, tag, comm);
   }

   *min = tmp[0];
   *max = tmp[1];
}
