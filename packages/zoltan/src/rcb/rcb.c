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
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "zz_const.h"
#include "rcb.h"
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

static int rcb_fn(ZZ *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **,
                  double, int, int, int, int, int, int, int, int);
static void print_rcb_tree(ZZ *, struct rcb_tree *);


/*  Parameters structure for RCB method.  Used in  */
/*  Zoltan_RCB_Set_Param and Zoltan_RCB.                   */
static PARAM_VARS RCB_params[] = {
                  { "RCB_OVERALLOC", NULL, "DOUBLE" },
                  { "RCB_REUSE", NULL, "INT" },
                  { "CHECK_GEOM", NULL, "INT" },
                  { "RCB_OUTPUT_LEVEL", NULL, "INT" },
                  { "KEEP_CUTS", NULL, "INT" },
                  { "RCB_LOCK_DIRECTIONS", NULL, "INT" },
                  { "RCB_SET_DIRECTIONS", NULL, "INT" },
                  { "RCB_RECTILINEAR_BLOCKS", NULL, "INT" },
                  { NULL, NULL, NULL } };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int Zoltan_RCB_Set_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, RCB_params, &result, &index);

    return(status);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int Zoltan_RCB(
  ZZ *zz,                       /* The Zoltan structure with info for
                                   the RCB balancer.                         */
  int *num_import,              /* Number of non-local objects assigned to this
                                   processor in the new decomposition.       */
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  int **import_procs,           /* Returned value:  array of processor IDs for
                                   processors owning the non-local objects in
                                   this processor's new decomposition.       */
  int **import_to_part,         /* Returned value:  array of partitions to
                                   which imported objects should be assigned. 
                                   KDDKDD  Currently assumes #parts == #procs.*/
  int *num_export,              /* Not computed, set to -1 */
  ZOLTAN_ID_PTR *export_global_ids, /* Not computed. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Not computed. */
  int **export_procs,           /* Not computed. */
  int **export_to_part          /* Not computed. */
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
    int reuse_dir;            /* (0) don't (1) reuse directions determined in
                                 the first iteration for future iterations. */
    int preset_dir;           /* Set order of directions: 0: don't set
                                 1: xyz,        2: xzy,      3: yzx,
                                 4: yxz,        5: zxy,      6: zyx  */
    int rectilinear_blocks;   /* (0) do (1) don't break ties in find_median */


    Zoltan_Bind_Param(RCB_params, "RCB_OVERALLOC", (void *) &overalloc);
    Zoltan_Bind_Param(RCB_params, "RCB_REUSE", (void *) &reuse);
    Zoltan_Bind_Param(RCB_params, "CHECK_GEOM", (void *) &check_geom);
    Zoltan_Bind_Param(RCB_params, "RCB_OUTPUT_LEVEL", (void *) &stats);
    Zoltan_Bind_Param(RCB_params, "KEEP_CUTS", (void *) &gen_tree);
    Zoltan_Bind_Param(RCB_params, "RCB_LOCK_DIRECTIONS", (void *) &reuse_dir);
    Zoltan_Bind_Param(RCB_params, "RCB_SET_DIRECTIONS", (void *) &preset_dir);
    Zoltan_Bind_Param(RCB_params, "RCB_RECTILINEAR_BLOCKS",
                              (void *) &rectilinear_blocks);

    overalloc = RCB_DEFAULT_OVERALLOC;
    reuse = RCB_DEFAULT_REUSE;
    check_geom = DEFAULT_CHECK_GEOM;
    stats = RCB_DEFAULT_OUTPUT_LEVEL;
    gen_tree = 0;
    wgtflag = (zz->Obj_Weight_Dim > 0); /* Multidim. weights not accepted */
    reuse_dir = 0;
    preset_dir = 0;
    rectilinear_blocks = 0;

    Zoltan_Assign_Param_Vals(zz->Params, RCB_params, zz->Debug_Level, zz->Proc,
                         zz->Debug_Proc);

    /* Initializations in case of early exit. */
    *num_import = -1;
    *num_export = -1;  /* We don't compute the export map. */

    return(rcb_fn(zz, num_import, import_global_ids, import_local_ids,
		 import_procs, import_to_part, overalloc, reuse, wgtflag,
                 check_geom, stats, gen_tree, reuse_dir, preset_dir,
                 rectilinear_blocks));
}

/*---------------------------------------------------------------------------*/

static int rcb_fn(
  ZZ *zz,                       /* The Zoltan structure with info for
                                   the RCB balancer.                         */
  int *num_import,              /* Number of non-local objects assigned to this
                                   processor in the new decomposition.       */
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                   non-local objects in this processor's new
                                   decomposition.                            */
  int **import_procs,           /* Returned value:  array of processor IDs for
                                   processors owning the non-local objects in
                                   this processor's new decomposition.       */
  int **import_to_part,         /* Returned value:  array of partitions
                                   to which objects are imported.
                                   KDDKDD Currently Assume #parts == #procs. */
  double overalloc,             /* amount to overallocate by when realloc
                                   of dot array must be done.     
                                   1.0 = no extra; 1.5 = 50% extra; etc.     */
  int reuse,                    /* (0) don't use (1) use previous cuts
                                   stored in treept at initial guesses.      */
  int wgtflag,                  /* (0) do not (1) do use weights.
                                   Multidimensional weights not supported    */
  int check_geom,               /* Check input & output for consistency?     */
  int stats,                    /* Print timing & count summary?             */
  int gen_tree,                 /* (0) do not (1) do generate full treept    */
  int reuse_dir,                /* (0) don't (1) reuse directions determined in
                                   the first iteration for future iterations. */
  int preset_dir,               /* Set order of directions:     0: don't set,
                                    1: xyz,        2: xzy,      3: yzx,
                                    4: yxz,        5: zxy,      6: zyx  */
  int rectilinear_blocks        /* (0) do (1) don't break ties in find_median*/
)
{
  char    yo[] = "rcb_fn";
  int     proc,nprocs;              /* my proc id, total # of procs */
  ZOLTAN_ID_PTR gidpt = NULL;           /* pointer to rcb->Global_IDs. */
  ZOLTAN_ID_PTR lidpt = NULL;           /* pointer to rcb->Local_IDs. */
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
  int     use_ids;                  /* When true, global and local IDs will be
                                       stored along with dots in the RCB_STRUCT.
                                       When false, storage, manipulation, and
                                       communication of IDs is avoided.     
                                       Set by call to Zoltan_RB_Use_IDs().         */

  RCB_STRUCT *rcb = NULL;           /* Pointer to data structures for RCB.  */
  struct rcb_box *rcbbox = NULL;    /* bounding box of final RCB sub-domain */
  struct rcb_tree *treept = NULL;   /* tree of RCB cuts - only single cut on 
                                       exit */

  double start_time, end_time;
  double lb_time[2];
  int    tmp_nprocs;                /* added for Tflops_Special */
  int    old_nprocs;                /* added for Tflops_Special */
  double weight;                    /* weight for current partition */
  double weightlo;                  /* weight of lower side of cut */
  double weighthi;                  /* weight of upper side of cut */
  int *dotlist;                     /* list of dots for find_median */
  int lock_direction = 0;           /* flag to determine direction after first
                                       iteration */
  int level;                        /* recursion level of RCB for preset_dir */
  int *dim_spec = NULL;             /* specified direction for preset_dir */
  int ix[3];                        /* temporaries for preset_dir */
  double wx, wy, wz;                /* width for preset_dir */

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;
  MPI_Op box_op;
  MPI_Datatype box_type;
  MPI_User_function Zoltan_RCB_box_merge;

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    MPI_Barrier(zz->Communicator);
    timestart = time1 = Zoltan_Time(zz->Timer);
  }

  /* setup for parallel */

  proc = zz->Proc;
  nprocs = zz->Num_Proc;

  /* initializations */

  /* 
   * Determine whether to store, manipulate, and communicate global and 
   * local IDs.
   */
  use_ids = Zoltan_RB_Use_IDs(zz);

  /*
   *  Build the RCB Data structure and 
   *  set pointers to information in it.
   */

  start_time = Zoltan_Time(zz->Timer);
  ierr = Zoltan_RCB_Build_Structure(zz, &pdotnum, &dotmax, wgtflag, use_ids);
  if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(proc, yo, 
                   "Error returned from Zoltan_RCB_Build_Structure.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }

  rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);

  gidpt  = rcb->Global_IDs;
  lidpt  = rcb->Local_IDs;
  dotpt  = rcb->Dots; 
  rcbbox = rcb->Box;
  treept = rcb->Tree_Ptr;
  end_time = Zoltan_Time(zz->Timer);
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
    dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
    if (dotmark == NULL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    coord = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
    if (coord == NULL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    wgts = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
    if (wgts == NULL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&coord);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
    if (dotlist == NULL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&wgts);
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&coord);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  }
  else {
    dotmark = NULL;
    coord = NULL;
    wgts = NULL;
    dotlist = NULL;
  }

  /* if reuse is turned on, turn on gen_tree since it is needed. */
  /* Also, if this is not first time through, send dots to previous proc. */
  if (reuse) {
    gen_tree = 1;

    if (treept[0].dim != -1) {
      /* find previous location of dots */
      for (outgoing = i = 0; i < dotnum; i++) {
        ierr = Zoltan_LB_Point_Assign(zz, dotpt[i].X, &dotmark[i]);
        if (dotmark[i] != proc) outgoing++;
      }

      if (outgoing)
        if ((proc_list = (int *) ZOLTAN_MALLOC(outgoing*sizeof(int))) == NULL) {
          ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
          ZOLTAN_FREE(&dotmark);
          ZOLTAN_FREE(&coord);
          ZOLTAN_FREE(&wgts);
          ZOLTAN_TRACE_EXIT(zz, yo);
          return ZOLTAN_MEMERR;
        }

      for (dottop = j = i = 0; i < dotnum; i++)
        if (dotmark[i] != proc)
          proc_list[j++] = dotmark[i];
        else
          dottop++;

      /* move dots */
      allocflag = 0;
      ierr = Zoltan_RB_Send_Dots(zz, &gidpt, &lidpt, &dotpt, dotmark, proc_list,
                             outgoing, &dotnum, &dotmax, proc, &allocflag,
                             overalloc, stats, reuse_count, use_ids,
                             zz->Communicator);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RB_Send_Dots.");
        ZOLTAN_FREE(&proc_list);
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_FREE(&coord);
        ZOLTAN_FREE(&wgts);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ierr);
      }

      if (allocflag) {
        /*
         * gidpt, lidpt and dotpt were reallocated in Zoltan_RB_Send_Dots;
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

      if (outgoing) ZOLTAN_FREE(&proc_list);
    }
  }

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);

  MPI_Op_create(&Zoltan_RCB_box_merge,1,&box_op);

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

  MPI_Allreduce(&boxtmp,rcbbox,1,box_type,box_op,zz->Communicator);

  /* if preset_dir is turned on, count number of levels of recursion,
     determine number of cuts in each direction, and then assign those
     cuts to an order according to order of directions */
  if (preset_dir) {
     if (preset_dir < 0 || preset_dir > 6) {
        ZOLTAN_PRINT_ERROR(proc, yo, 
                       "Error: parameter RCB_SET_DIRECTIONS out of bounds");
        preset_dir = 1;
     }
     wx = rcbbox->hi[0] - rcbbox->lo[0];
     wy = rcbbox->hi[1] - rcbbox->lo[1];
     wz = rcbbox->hi[2] - rcbbox->lo[2];
     for (i = 0; i < 3; i++) ix[i] = 0;
     for (tmp_nprocs = nprocs, level = 0; tmp_nprocs > 1; level++) {
        tmp_nprocs = (tmp_nprocs + 1)/2;
        if (wz > wx && wz > wy) {
           ix[2]++;
           wz /= 2.0;
        }
        else
           if (wy > wx && wy > wz) {
              ix[1]++;
              wy /= 2.0;
           }
           else {
              ix[0]++;
              wx /= 2.0;
           }
     }
     dim_spec = (int *) ZOLTAN_MALLOC(level*sizeof(int));
     if (dim_spec == NULL) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
        ZOLTAN_FREE(&wgts);
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_FREE(&coord);
        ZOLTAN_FREE(&dotlist);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
     }
     for (i = j = 0; i < level; i++) {
        if (j == 0) {
           if (preset_dir < 3)
              dim = 0;
           else if (preset_dir < 5)
              dim = 1;
           else
              dim = 2;
           if (ix[dim] == 0)
              j++;
        }
        if (j == 1) {
           if (preset_dir == 1 || preset_dir == 6)
              dim = 1;
           else if (preset_dir < 4)
              dim = 2;
           else
              dim = 0;
           if (ix[dim] == 0)
              j++;
        }
        if (j == 2) {
           if (preset_dir == 3 || preset_dir == 6)
              dim = 0;
           else if (preset_dir == 2 || preset_dir == 5)
              dim = 1;
           else
              dim = 0;
        }
        dim_spec[i] = dim;
        ix[dim]--;
     }
  }

  /* if reuse_dir is turned on, turn on gen_tree since it is needed. */
  /* Also, if this is not first time through, set lock_direction to reuse
     the directions determined previously */
  if (reuse_dir) {
     gen_tree = 1;
     if (treept[0].dim != -1)
        lock_direction = 1;
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
  level = 0;

  while (num_procs > 1 || (zz->Tflops_Special && tmp_nprocs > 1)) {
    /* tmp_nprocs is size of largest partition - force all processors to go
       through all levels of rcb */
    if (zz->Tflops_Special) tmp_nprocs = (tmp_nprocs + 1)/2;

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) time1 = Zoltan_Time(zz->Timer);

    ierr = Zoltan_Divide_Machine(zz, proc, local_comm, &set, &proclower,
                             &procmid, &num_procs, &fractionlo);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error in Zoltan_Divide_Machine.");
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&coord);
      ZOLTAN_FREE(&wgts);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
    
    /* dim = dimension (xyz = 012) to bisect on */

    if (lock_direction)
       dim = treept[procmid].dim;
    else if (preset_dir)
       dim = dim_spec[level++];
    else {
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
    }
    
    /* create mark array and active list for dots */

    if (allocflag) {
      allocflag = 0;
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&coord);
      ZOLTAN_FREE(&wgts);
      ZOLTAN_FREE(&dotlist);
      dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
      if (dotmark == NULL) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      coord = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
      if (coord == NULL) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      wgts = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double));
      if (wgts == NULL) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_FREE(&coord);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int));
      if (dotlist == NULL) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
        ZOLTAN_FREE(&wgts);
        ZOLTAN_FREE(&dotmark);
        ZOLTAN_FREE(&coord);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
    }

    /* copy correct coordinate value into the temporary array */
    for (i = 0; i < dotnum; i++) {
      coord[i] = dotpt[i].X[dim];
      wgts[i] = dotpt[i].Weight;
    }

    /* determine if there is a first guess to use */
    /* The test on old_nprocs is for the TFLOPS_SPECIAL flag */
    if (old_nprocs > 1 && reuse && dim == treept[procmid].dim) {
      if (stats) counters[5]++;
      valuehalf = treept[procmid].cut;
      first_guess = 1;
    }
    else first_guess = 0;

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time2 = Zoltan_Time(zz->Timer);

    if (!Zoltan_RB_find_median(
                zz->Tflops_Special, coord, wgts, dotmark, dotnum, proc, 
                fractionlo, local_comm, &valuehalf, first_guess, &(counters[0]),
                nprocs, old_nprocs, proclower, wgtflag, rcbbox->lo[dim],
                rcbbox->hi[dim], weight, &weightlo, &weighthi,
                dotlist, rectilinear_blocks)) {
      ZOLTAN_PRINT_ERROR(proc, yo,"Error returned from Zoltan_RB_find_median.");
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&coord);
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
    ierr = Zoltan_RB_Send_Outgoing(zz, &gidpt, &lidpt, &dotpt, dotmark, &dottop,
                               &dotnum, &dotmax, set, &allocflag, overalloc,
                               stats, counters, use_ids, local_comm, proclower,
                               old_nprocs);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RB_Send_Outgoing.");
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&coord);
      ZOLTAN_FREE(&wgts);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }

    if (allocflag) {
      /* 
       * gidpt, lidpt and dotpt were reallocated in Zoltan_RB_Send_Outgoing;
       * store their values in rcb.
       */
      rcb->Global_IDs = gidpt;
      rcb->Local_IDs = lidpt;
      rcb->Dots = dotpt;
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

  /* free all memory used by RCB and MPI */

  if (!zz->Tflops_Special) MPI_Comm_free(&local_comm);
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);

  ZOLTAN_FREE(&coord);
  ZOLTAN_FREE(&wgts);
  ZOLTAN_FREE(&dotmark);
  ZOLTAN_FREE(&dotlist);
  if (preset_dir) ZOLTAN_FREE(&dim_spec);

  end_time = Zoltan_Time(zz->Timer);
  lb_time[1] = end_time - start_time;

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    MPI_Barrier(zz->Communicator);
    timestop = time1 = Zoltan_Time(zz->Timer);
  }

  /* error checking and statistics */

  if (check_geom) {
    ierr = Zoltan_RB_check_geom_output(zz, dotpt,dotnum,pdotnum,rcbbox);
    if (ierr == ZOLTAN_FATAL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RB_check_geom_output");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ierr);
    }
  }

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
    Zoltan_RB_stats(zz, timestop-timestart,dotpt,dotnum,
                timers,counters,stats,reuse_count,rcbbox,reuse);

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
         ZOLTAN_PRINT_ERROR(proc,yo,"Error returned from Zoltan_RB_Return_Arguments.");
         ZOLTAN_TRACE_EXIT(zz, yo);
         return ierr;
      }
    }
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    print_rcb_tree(zz, &(treept[proc]));

  if (gen_tree) {
    MPI_Allgather(&treept[proc], sizeof(struct rcb_tree), MPI_BYTE,
      treept, sizeof(struct rcb_tree), MPI_BYTE, zz->Communicator);
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

  end_time = Zoltan_Time(zz->Timer);
  lb_time[0] += (end_time - start_time);

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
    if (zz->Proc == zz->Debug_Proc) {
      printf("ZOLTAN RCB Times:  \n");
    }
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, lb_time[0], 
                   "ZOLTAN     Build:       ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, lb_time[1], 
                   "ZOLTAN     RCB:         ");
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    /* zz->Debug_Level >= ZOLTAN_DEBUG_ALL ==> use_ids is true */
    Zoltan_RB_Print_All(zz, rcb->Global_IDs, rcb->Dots, 
                    pdotnum, pdottop, *num_import,
                    *import_global_ids, *import_procs);
  }

  /* Free memory allocated by the algorithm. */
  if (!reuse && !gen_tree) {
    /* Free all memory used. */
    Zoltan_RCB_Free_Structure(zz);
  }
  else {
    /* Free only Dots and IDs; keep other structures. */
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    ZOLTAN_FREE(&(rcb->Dots));
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  /* Temporary return value until error codes are fully implemented */
  return(ZOLTAN_OK);  
}



/* ----------------------------------------------------------------------- */

/* MPI user-defined reduce operations */

/* min/max merge of each component of a rcb_box */

void Zoltan_RCB_box_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

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

static void print_rcb_tree(ZZ *zz, struct rcb_tree *treept)
{
  Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
  printf("Proc %d:  Tree Struct:\n", zz->Proc);
  printf("          cut        = %e\n", treept->cut);
  printf("          dim        = %d\n", treept->dim);
  printf("          parent     = %d\n", treept->parent);
  printf("          left_leaf  = %d\n", treept->left_leaf);
  printf("          right_leaf = %d\n", treept->right_leaf);
  Zoltan_Print_Sync_End(zz->Communicator, TRUE);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
