/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    Revision: 1.6.2.2 $
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "zz_const.h"
#include "rib.h"
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


static int rib_fn(ZZ *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **,
                  double, int, int, int, int, float *);
static void print_rib_tree(ZZ *, struct rib_tree *);
/* for Tflops_Special */
static void Zoltan_RIB_min_max(double *, double *, int, int, int, MPI_Comm);

/*  Parameters structure for RIB method.  Used in  */
/*  Zoltan_RIB_Set_Param and Zoltan_RIB.                      */
static PARAM_VARS RIB_params[] = {
               { "RIB_OVERALLOC", NULL, "DOUBLE", 0 },
               { "CHECK_GEOM", NULL, "INT", 0 },
               { "RIB_OUTPUT_LEVEL", NULL, "INT", 0 },
               { "KEEP_CUTS", NULL, "INT", 0 },
               { NULL, NULL, NULL, 0 } };

/*---------------------------------------------------------------------------*/

int Zoltan_RIB_Set_Param(
  char *name,                 /* name of variable */
  char *val                   /* value of variable */
)
{
  int status;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */

  status = Zoltan_Check_Param(name, val, RIB_params, &result, &index);

  return(status);
}

/*---------------------------------------------------------------------------*/

int Zoltan_RIB(
  ZZ *zz,                       /* The Zoltan structure with info for
                                the RIB balancer.                       */
  float *part_sizes,            /* Input:  Array of size zz->Num_Global_Parts
                                containing the percentage of work to be
                                assigned to each partition.               */
  int *num_import,              /* Returned value: Number of non-local 
                                objects assigned to
                                this processor in the new decomposition.*/
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value: array of global IDs for
                                non-local objects in this processor's new
                                decomposition.                          */
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value: array of local IDs for
                                non-local objects in this processor's new
                                decomposition.                          */
  int **import_procs,           /* Returned value: array of processor IDs for
                                processors owning the non-local objects in
                                this processor's new decomposition.     */
  int **import_to_part,         /* Returned value:  array of partitions to
                                 which imported objects should be assigned.
                                 KDDKDD  Currently assumes #parts == #proc. */
  int *num_export,              /* Not computed, set to -1 */
  ZOLTAN_ID_PTR *export_global_ids, /* Not computed. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Not computed. */
  int **export_procs,           /* Not computed. */
  int **export_to_part          /* Not computed. */
)
{
  /* Wrapper routine to set parameter values and call the real rib. */
  char *yo = "Zoltan_RIB";
  double overalloc;           /* amount to overallocate by when realloc
                              of dot array must be done.
                              1.0 = no extra; 1.5 = 50% extra; etc. */
  int wgtflag;                /* (0) do not (1) do use weights.
                              Multidimensional weights not supported */
  int check_geom;             /* Check input & output for consistency? */
  int stats;                  /* Print timing & count summary? */
  int gen_tree;               /* (0) don't (1) generate whole treept to use
                              later for point and box drop. */
  float *work_fraction=NULL;/* Partition sizes -- one weight per partition */
  int i, ierr;

  Zoltan_Bind_Param(RIB_params, "RIB_OVERALLOC", (void *) &overalloc);
  Zoltan_Bind_Param(RIB_params, "CHECK_GEOM", (void *) &check_geom);
  Zoltan_Bind_Param(RIB_params, "RIB_OUTPUT_LEVEL", (void *) &stats);
  Zoltan_Bind_Param(RIB_params, "KEEP_CUTS", (void *) &gen_tree);

  overalloc = RIB_DEFAULT_OVERALLOC;
  check_geom = DEFAULT_CHECK_GEOM;
  stats = RIB_DEFAULT_OUTPUT_LEVEL;
  gen_tree = 0;
  wgtflag = (zz->Obj_Weight_Dim > 0); /* Multidim. weights not accepted */

  Zoltan_Assign_Param_Vals(zz->Params, RIB_params, zz->Debug_Level, zz->Proc,
                    zz->Debug_Proc);

  /* Initializations in case of early exit. */
  *num_import = -1;
  *num_export = -1;  /* We don't compute the export map. */

  if (zz->Obj_Weight_Dim <= 1)
    work_fraction = part_sizes;
  else {
    /* Multi-dimensional partition sizes not permitted yet;
       reduce to one weight per partition.  
       This reduction will not be needed once multi-constraint partitioning
       is implemented. */
    work_fraction = (float *) ZOLTAN_MALLOC(sizeof(float) *
                                            zz->LB.Num_Global_Parts);
    if (work_fraction == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      return(ZOLTAN_MEMERR);
    }

    /* RIB supports only the first weight; pick out the appropriate
       part_sizes entries for the first weight. */
    for (i = 0 ; i < zz->LB.Num_Global_Parts ; i++)
       work_fraction[i] = part_sizes[i*zz->Obj_Weight_Dim];
  }

  ierr = rib_fn(zz, num_import, import_global_ids, import_local_ids,
                import_procs, import_to_part,
                overalloc, wgtflag, check_geom, stats, gen_tree, work_fraction);

  if (zz->Obj_Weight_Dim > 1) ZOLTAN_FREE(&work_fraction);
  return(ierr);

}

/*---------------------------------------------------------------------------*/

static int rib_fn(
  ZZ *zz,                       /* The Zoltan structure with info for
                                the RIB balancer. */
  int *num_import,              /* Number of non-local objects assigned to
                                this processor in the new decomposition.*/
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                non-local objects in this processor's new
                                decomposition. */
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                non-local objects in this processor's new
                                decomposition. */
  int **import_procs,           /* Returned value: array of processor IDs for
                                processors owning the non-local objects in
                                this processor's new decomposition. */
  int **import_to_part,         /* Returned value: array of partitions to
                                which objects are imported.
                                KDDKDD Currently assumes #parts == #procs */
  double overalloc,             /* amount to overallocate by when realloc
                                of dot array must be done.
                                  1.0 = no extra; 1.5 = 50% extra; etc. */
  int wgtflag,                  /* (0) do not (1) do use weights.
                                Multidimensional weights not supported */
  int check_geom,               /* Check input & output for consistency? */
  int stats,                    /* Print timing & count summary? */
  int gen_tree,                 /* (0) do not (1) do generate full treept */
  float *part_sizes             /* Input:  Array of size zz->Num_Global_Parts
                                containing the percentage of work to be
                                assigned to each partition.               */
)
{
  char    yo[] = "rib_fn";
  int     proc,nprocs;        /* my proc id, total # of procs */
  ZOLTAN_ID_PTR gidpt = NULL;     /* local global IDs array. */
  ZOLTAN_ID_PTR lidpt = NULL;     /* local local IDs array. */
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
                                 Set by call to Zoltan_RB_Use_IDs().         */


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
  int rectilinear_blocks = 0; /* parameter for find_median (not used by rib) */

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    MPI_Barrier(zz->Communicator);
    timestart = time1 = Zoltan_Time(zz->Timer);
  }

  /* setup for parallel */

  proc = zz->Proc;
  nprocs = zz->Num_Proc;

  /*
   * Determine whether to store, manipulate, and communicate global and
   * local IDs.
   */
  use_ids = Zoltan_RB_Use_IDs(zz);

  /*
   *  Build the RIB Data structure and
   *  set pointers to information in it.
   */

  start_time = Zoltan_Time(zz->Timer);
  ierr = Zoltan_RIB_Build_Structure(zz, &pdotnum, &dotmax, wgtflag, use_ids);
  if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RIB_Build_Structure.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }

  rib = (RIB_STRUCT *) (zz->LB.Data_Structure);

  gidpt = rib->Global_IDs;
  lidpt = rib->Local_IDs;
  dotpt  = rib->Dots;
  treept = rib->Tree_Ptr;
  end_time = Zoltan_Time(zz->Timer);
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
    dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
    if (dotmark == NULL) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    value = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
    if (value == NULL) {
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    wgts = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
    if (wgts == NULL) {
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
    if (dotlist == NULL) {
      ZOLTAN_FREE(&wgts);
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
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
  MPI_Allreduce(&weightlo, &weight, 1, MPI_DOUBLE, MPI_SUM, zz->Communicator);

  if (check_geom) {
    ierr = Zoltan_RB_check_geom_input(zz, dotpt, dotnum);
    if (ierr == ZOLTAN_FATAL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RB_check_geom_input");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ierr);
    }
  }

  /* create local communicator for use in recursion */

  if (zz->Tflops_Special)
     local_comm = zz->Communicator;
  else
     MPI_Comm_dup(zz->Communicator,&local_comm);

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    time2 = Zoltan_Time(zz->Timer);
    timers[0] = time2 - time1;
  }

  /* recursively halve until just one proc in partition */

  old_nprocs = num_procs = nprocs;
  root = 0;
  old_set = 1;
  treept[proc].parent = 0;
  treept[proc].left_leaf = 0;
  if (zz->Tflops_Special) {
     proclower = 0;
     tmp_nprocs = nprocs;
  }

  while (num_procs > 1 || (zz->Tflops_Special && tmp_nprocs > 1)) {
    /* tmp_nprocs is size of largest partition - force all processors to go
       through all levels of rib */
    if (zz->Tflops_Special) tmp_nprocs = (tmp_nprocs + 1)/2;

    ierr = Zoltan_Divide_Machine(zz, part_sizes, proc, local_comm, 
                                 &set, &proclower,
                                 &procmid, &num_procs, &fractionlo);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_FREE(&wgts);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }

    /* create mark array and active list for dots */

    if (allocflag) {
      allocflag = 0;
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_FREE(&wgts);
      ZOLTAN_FREE(&dotlist);
      dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
      if (dotmark == NULL) {
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      value = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
      if (value == NULL) {
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      wgts = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
      if (wgts == NULL) {
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_FREE(&value);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
      if (dotlist == NULL) {
        ZOLTAN_FREE(&wgts);
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_FREE(&value);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
    }

    for (i = 0; i < dotnum; i++) {
      wgts[i] = dotpt[i].Weight;
    }
    if (old_nprocs > 1) {    /* for Tflops_Special */
       switch (rib->Num_Geom) {
          case 3:
             ierr = Zoltan_RIB_inertial3d(zz, dotpt, dotnum, wgtflag, cm, evec, value,
                                  local_comm, proc, old_nprocs, proclower);
             break;
          case 2:
             ierr = Zoltan_RIB_inertial2d(zz, dotpt, dotnum, wgtflag, cm, evec, value,
                                  local_comm, proc, old_nprocs, proclower);
             break;
          case 1:
             ierr = Zoltan_RIB_inertial1d(dotpt, dotnum, wgtflag, cm, evec, value);
             break;
       }
       valuelo = MYHUGE;
       valuehi = -MYHUGE;
       for (i = 0; i < dotnum; i++) {
          if (value[i] < valuelo) valuelo = value[i];
          if (value[i] > valuehi) valuehi = value[i];
       }
       if (zz->Tflops_Special)
          Zoltan_RIB_min_max(&valuelo, &valuehi, proclower, proc, old_nprocs,
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

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time2 = Zoltan_Time(zz->Timer);

    if (!Zoltan_RB_find_median(
                   zz->Tflops_Special, value, wgts, dotmark, dotnum, proc, 
                   fractionlo, local_comm, &valuehalf, first_guess,
                   &(counters[0]), nprocs, old_nprocs, proclower,
                   wgtflag, valuelo, valuehi, weight, &weightlo,
                   &weighthi, dotlist, rectilinear_blocks)) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from find_median.");
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_FREE(&wgts);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_FATAL;
    }

    if (set)    /* set weight for current partition */
      weight = weighthi;
    else
      weight = weightlo;

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time3 = Zoltan_Time(zz->Timer);

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

    ierr = Zoltan_RB_Send_Outgoing(zz, &gidpt, &lidpt, &dotpt, dotmark, &dottop,
                               &dotnum, &dotmax, set, &allocflag, overalloc,
                               stats, counters, use_ids, local_comm, proclower,
                               old_nprocs);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RB_Send_Outgoing.");
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_FREE(&wgts);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ierr;
    }
    
    if (allocflag) {
      /* 
       * gidpt, lidpt and dotpt were reallocated in Zoltan_RB_Send_Outgoing;
       * store their values in rib.
       */
      rib->Global_IDs = gidpt;
      rib->Local_IDs = lidpt;
      rib->Dots = dotpt;
    }

    /* create new communicators */

    if (zz->Tflops_Special) {
       if (set) proclower = procmid;
       old_nprocs = num_procs;
    }
    else {
       MPI_Comm_split(local_comm,set,proc,&tmp_comm);
       MPI_Comm_free(&local_comm);
       local_comm = tmp_comm;
    }

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
      time4 = Zoltan_Time(zz->Timer);
      timers[1] += time2 - time1;
      timers[2] += time3 - time2;
      timers[3] += time4 - time3;
    }
  }

  /* have recursed all the way to final single sub-domain */

  /* free all memory used by RIB and MPI */

  if (!zz->Tflops_Special) MPI_Comm_free(&local_comm);

  ZOLTAN_FREE(&value);
  ZOLTAN_FREE(&wgts);
  ZOLTAN_FREE(&dotmark);
  ZOLTAN_FREE(&dotlist);

  end_time = Zoltan_Time(zz->Timer);
  lb_time[1] = end_time - start_time;

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    MPI_Barrier(zz->Communicator);
    timestop = time1 = Zoltan_Time(zz->Timer);
  }

  /* error checking and statistics */

  if (check_geom) {
    ierr = Zoltan_RB_check_geom_output(zz, dotpt, part_sizes, 
                                       dotnum, pdotnum, NULL);
    if (ierr == ZOLTAN_FATAL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RB_check_geom_output");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ierr);
    }
  }

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME))
    Zoltan_RB_stats(zz, timestop-timestart, dotpt, dotnum, timers, counters,
                stats, NULL, NULL, FALSE);

  /* update calling routine parameters */

  start_time = Zoltan_Time(zz->Timer);

  pdotnum = dotnum;
  pdottop = dottop;

  /*  build return arguments */

  if (zz->LB.Return_Lists) {
    /* zz->LB.Return_Lists is true ==> use_ids is true */
    *num_import = dotnum - dottop;
    if (*num_import > 0) {
      ierr = Zoltan_RB_Return_Arguments(zz, gidpt, lidpt, dotpt, *num_import,
                                    import_global_ids, import_local_ids,
                                    import_procs, import_to_part, dottop);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(proc, yo,
                       "Error returned from Zoltan_RB_Return_Arguments.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ierr;
      }
    }
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    print_rib_tree(zz, &(treept[proc]));

  if (gen_tree) {
    MPI_Allgather(&treept[proc], sizeof(struct rib_tree), MPI_BYTE,
                  treept, sizeof(struct rib_tree), MPI_BYTE,
                  zz->Communicator);
    for (i = 1; i < nprocs; i++)
      if (treept[i].parent > 0)
        treept[treept[i].parent - 1].left_leaf = i;
      else if (treept[i].parent < 0)
        treept[-treept[i].parent - 1].right_leaf = i;
  }
  else {
    treept[0].right_leaf = -1;
  }

  end_time = Zoltan_Time(zz->Timer);
  lb_time[0] += (end_time - start_time);

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
    if (zz->Proc == zz->Debug_Proc)
      printf("ZOLTAN RIB Times:  \n");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, lb_time[0], 
                   "ZOLTAN       Build:       ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, lb_time[1], 
                   "ZOLTAN         RIB:         ");
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    /* zz->Debug_Level >= ZOLTAN_DEBUG_ALL ==> use_ids is true */
    Zoltan_RB_Print_All(zz, rib->Global_IDs, rib->Dots, 
                    pdotnum, pdottop, *num_import, 
                    *import_global_ids, *import_procs);
  }

  /* Free memory allocated by the algorithm. */
  if (!gen_tree) {
    /* Free all memory used. */
    Zoltan_RIB_Free_Structure(zz);
  }
  else {
    /* Free only Dots and IDs; keep other structures. */
    ZOLTAN_FREE(&(rib->Global_IDs));
    ZOLTAN_FREE(&(rib->Local_IDs));
    ZOLTAN_FREE(&(rib->Dots));
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  /* Temporary return value until error codes are fully implemented */
  return(ZOLTAN_OK);
}

static void print_rib_tree(ZZ *zz, struct rib_tree *treept)
{
  Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
  printf("Proc %d:  Tree Struct:\n", zz->Proc);
  printf("          cm         = (%e,%e,%e)\n", 
         treept->cm[0], treept->cm[1], treept->cm[2]);
  printf("          ev         = (%e,%e,%e)\n", 
         treept->ev[0], treept->ev[1], treept->ev[2]);
  printf("          cut        = %e\n", treept->cut);
  printf("          parent     = %d\n", treept->parent);
  printf("          left_leaf  = %d\n", treept->left_leaf);
  printf("          right_leaf = %d\n", treept->right_leaf);
  Zoltan_Print_Sync_End(zz->Communicator, TRUE);
}

static void Zoltan_RIB_min_max(
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

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
