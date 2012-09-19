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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <float.h>
#include "zz_const.h"
#include "rcb.h"
#include "params_const.h"
#include "ha_const.h"
#include "par_median_const.h"
#include "par_bisect_const.h"
#include "rcb_params.h"

/* Recursive coordinate bisectioning (RCB) routine
   operates on "dots" as defined in shared_const.h
*/

/* Original version was written by
   Steve Plimpton, Sandia National Labs, ABQ, NM  87185
   Dept 9221, MS 1111
   (505) 845-7873
   sjplimp@cs.sandia.gov

   Heavily modified and enhanced for Zoltan.
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
/*****************************************************************************/

/*  RCB_DEFAULT_OUTPUT_LEVEL = 0  No statistics logging */
/*  RCB_DEFAULT_OUTPUT_LEVEL = 1  Log times and counts, print summary */
/*  RCB_DEFAULT_OUTPUT_LEVEL = 2  Log times and counts, print for each proc */
#define RCB_DEFAULT_OUTPUT_LEVEL 0
#define RCB_DEFAULT_OVERALLOC 1.2
#define RCB_DEFAULT_REUSE FALSE

/* The median of an array of floating point values is found by iterating 
 * through candidates.  We have two methods for choosing candidates:
 *
 *   BISECTION - Candidate within region [L, R] is value closest to (L+R)/2
 *
 *   RANDOM - First candidate is the median of the collection of local
 *     medians on each process.  Then we do binary search through the
 *     remaining local medians.  If median is still not found, each process
 *     supplies a candidate chosen "randomly" from its values in the
 *     region of interest, and we do a binary search through these.  Repeat
 *     this last step with narrowed region until median is found.
 *
 * BISECTION works very well if the geometry is symmetric, because the initial
 *  guess is very close to the true median.  
 *
 * RANDOM attempts to do less communication at the expense of more local
 *  computation, and is a good choice when running on large numbers of
 *  processors, unless the geometry is very symmetric.
 *
 * For a given size of array, the number of iterations required by BISECTION
 * to find the median can vary greatly, from remarkably few to remarkably many,
 * just due to the layout of the values.  In contrast, the number of iterations
 * required by RANDOM on same-sized arrays varies much less.
 */
#define PIVOT_CHOICE_BISECTION 1
#define PIVOT_CHOICE_RANDOM    2

/*****************************************************************************/
/* function prototypes */

static int rcb_fn(ZZ *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **,
  int *, ZOLTAN_ID_PTR *,
  double, int, int, int, int, int, int, int, int, int, int, double, int, int,
  int, float *);
static void print_rcb_tree(ZZ *, int, int, struct rcb_tree *);
static int cut_dimension(int, struct rcb_tree *, int, int, int *, int *, 
  struct rcb_box *);
static int set_preset_dir(int, int, int, struct rcb_box *, int **);
static int serial_rcb(ZZ *, struct Dot_Struct *, int *, int *, int, int,
  struct rcb_box *, double *, int, int, int *, int *, int, int, int, int,
  int, int, int, int, int, int, int, int, MPI_Op, MPI_Datatype,
  int, ZOLTAN_GNO_TYPE *, struct rcb_tree *, int *, int, double *, double *,
  float *, double *, int, int, double, double *);
static void compute_RCB_box(struct rcb_box *, int, struct Dot_Struct *, int *,
  MPI_Op, MPI_Datatype, MPI_Comm, int, int, int, int);

/*****************************************************************************/

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
  float *part_sizes,            /* Input:  Array of size zz->LB.Num_Global_Parts
                                   * zz->Obj_Weight_Dim
                                   containing the percentage of work to be
                                   assigned to each part.               */
  int *num_import,              /* Returned value:  Number of non-local 
                                   objects assigned to this
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
  int **import_to_part,         /* Returned value:  array of parts to
                                   which imported objects should be assigned. */
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
    int wgtflag;              /* no. of weights per dot. */
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
    int obj_wgt_comp;         /* 1 if all (multi-)weights for an object 
                                 have same units (comparable), 0 otherwise. */
    int mcnorm;               /* norm (1,2,3) to use in multicriteria alg. */
    double max_aspect_ratio;  /* maximum allowed ratio of box dimensions */
    int recompute_box;        /* (0) do (1) don't recompute bounding box for
                                  partition sets at each level of recursion */
    int average_cuts;         /* Flag forcing median line to be drawn halfway
                                 between two closest objects. */
    int pivot_choice;
    int idummy;
    int final_output;
    int ierr=ZOLTAN_OK;
    double ddummy;


    Zoltan_Bind_Param(RCB_params, "RCB_OVERALLOC", (void *) &overalloc);
    Zoltan_Bind_Param(RCB_params, "RCB_REUSE", (void *) &reuse);
    Zoltan_Bind_Param(RCB_params, "CHECK_GEOM", (void *) &check_geom);
    Zoltan_Bind_Param(RCB_params, "RCB_OUTPUT_LEVEL", (void *) &stats);
    Zoltan_Bind_Param(RCB_params, "FINAL_OUTPUT", (void *) &final_output);
    Zoltan_Bind_Param(RCB_params, "KEEP_CUTS", (void *) &gen_tree);
    Zoltan_Bind_Param(RCB_params, "RCB_LOCK_DIRECTIONS", (void *) &reuse_dir);
    Zoltan_Bind_Param(RCB_params, "RCB_SET_DIRECTIONS", (void *) &preset_dir);
    Zoltan_Bind_Param(RCB_params, "RCB_RECTILINEAR_BLOCKS",
                              (void *) &rectilinear_blocks);
    Zoltan_Bind_Param(RCB_params, "OBJ_WEIGHTS_COMPARABLE",
                              (void *) &obj_wgt_comp);
    Zoltan_Bind_Param(RCB_params, "RCB_MULTICRITERIA_NORM",
                              (void *) &mcnorm);
    Zoltan_Bind_Param(RCB_params, "RCB_MAX_ASPECT_RATIO",
                              (void *) &max_aspect_ratio);
    Zoltan_Bind_Param(RCB_params, "RCB_RECOMPUTE_BOX",
                              (void *) &recompute_box);
    Zoltan_Bind_Param(RCB_params, "REDUCE_DIMENSIONS",
                              (void *) &idummy);
    Zoltan_Bind_Param(RCB_params, "DEGENERATE_RATIO",
                              (void *) &ddummy);
    Zoltan_Bind_Param(RCB_params, "AVERAGE_CUTS",
                              (void *) &average_cuts);
    Zoltan_Bind_Param(RCB_params, "RANDOM_PIVOTS",
                              (void *) &pivot_choice);

    /* Set default values. */
    overalloc = RCB_DEFAULT_OVERALLOC;
    reuse = RCB_DEFAULT_REUSE;
    check_geom = DEFAULT_CHECK_GEOM;
    stats = RCB_DEFAULT_OUTPUT_LEVEL;
    gen_tree = 0;
    wgtflag = zz->Obj_Weight_Dim;
    reuse_dir = 0;
    preset_dir = 0;
    rectilinear_blocks = 0;
    obj_wgt_comp = 0;     
    mcnorm = 1; 
    max_aspect_ratio = 10.;
    recompute_box = 0;
    idummy = 0;
    final_output = 0;
    ddummy = 0.0;
    average_cuts = 0;
    pivot_choice = 0;

    Zoltan_Assign_Param_Vals(zz->Params, RCB_params, zz->Debug_Level, zz->Proc,
                         zz->Debug_Proc);

    /* Need to keep the tree if reuse or reuse_dir, so turn gen_tree on. */
    if (reuse || reuse_dir) gen_tree = 1;

    /* Initializations in case of early exit. */
    *num_import = -1;
    *num_export = -1;  /* We don't compute the export map. */

    if (pivot_choice == 0){
      pivot_choice = PIVOT_CHOICE_BISECTION;
    }
    else{
      pivot_choice = PIVOT_CHOICE_RANDOM;
    }

    if (final_output && (stats < 1)){
      /* FINAL_OUTPUT is a graph/phg param, corresponds to our OUTPUT_LEVEL 1 */
      stats = 1;
    }

    ierr = rcb_fn(zz, num_import, import_global_ids, import_local_ids,
		 import_procs, import_to_part, num_export, export_global_ids,
                 overalloc, reuse, wgtflag,
                 check_geom, stats, gen_tree, reuse_dir, preset_dir,
                 rectilinear_blocks, obj_wgt_comp, mcnorm, 
                 max_aspect_ratio, recompute_box, average_cuts, pivot_choice,
                 part_sizes);

    return(ierr);
}

/*---------------------------------------------------------------------------*/

static int rcb_fn(
  ZZ *zz,                       /* The Zoltan structure with info for
                                   the RCB balancer.                         */
  int *num_import,              /* Number of non-local objects assigned to this
                                   processor in the new decomposition.      
                                   When LB.Return_Lists==CANDIDATE_LISTS,
                                   num_import returns the number of input
                                   objects as given by ZOLTAN_NUM_OBJ_FN. */
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                   non-local objects in this processor's new
                                   decomposition. 
                                   When LB.Return_Lists==CANDIDATE_LISTS,
                                   this array contains GIDs for all input 
                                   objs as given by ZOLTAN_OBJ_LIST_FN. */
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                   non-local objects in this processor's new
                                   decomposition.                           
                                   When LB.Return_Lists==CANDIDATE_LISTS,
                                   this array contains LIDs for all input 
                                   objs as given by ZOLTAN_OBJ_LIST_FN. */
  int **import_procs,           /* Returned value:  array of processor IDs for
                                   processors owning the non-local objects in
                                   this processor's new decomposition.      
                                   When LB.Return_Lists==CANDIDATE_LISTS,
                                   the returned array is NULL. */
  int **import_to_part,         /* Returned value:  array of parts
                                   to which objects are imported.           
                                   When LB.Return_Lists==CANDIDATE_LISTS,
                                   the returned array is NULL. */
  int *num_export,              /* Returned value only when 
                                   LB.Return_Lists==CANDIDATE_LISTS; number of
                                   input objs as given by ZOLTAN_NUM_OBJ_FN */
  ZOLTAN_ID_PTR *export_global_ids, /* Returned value only when
                                   LB.Return_Lists==CANDIDATE_LISTS; for each
                                   input obj (from ZOLTAN_OBJ_LIST_FN), 
                                   return a candidate obj from the part to which
                                   the obj is assigned; used in PHG matching */
  double overalloc,             /* amount to overallocate by when realloc
                                   of dot array must be done.     
                                   1.0 = no extra; 1.5 = 50% extra; etc.     */
  int reuse,                    /* (0) don't use (1) use previous cuts
                                   stored in treept at initial guesses.      */
  int wgtflag,                  /* No. of weights per dot provided by user    */
  int check_geom,               /* Check input & output for consistency?     */
  int stats,                    /* Print timing & count summary?             */
  int gen_tree,                 /* (0) do not (1) do generate full treept    */
  int reuse_dir,                /* (0) don't (1) reuse directions determined in
                                   the first iteration for future iterations. */
  int preset_dir,               /* Set order of directions:     0: don't set,
                                    1: xyz,        2: xzy,      3: yzx,
                                    4: yxz,        5: zxy,      6: zyx  */
  int rectilinear_blocks,       /* (0) do (1) don't break ties in find_median*/
  int obj_wgt_comp,             /* (0) obj wgts not comparable, use autoscaling
                                   (1) obj wgts have same units, no scaling */
  int mcnorm,                   /* norm (1-3) for multicriteria bisection */
  double max_aspect_ratio,      /* max allowed ratio of box dimensions */
  int recompute_box,            /* (0) do (1) don't recompute bounding box for
                                   partition sets at each level of recursion */
  int average_cuts,             /* Flag forcing median line to be drawn halfway
                                   between two closest objects. */
  int pivot_choice, 
  float *part_sizes             /* Input: Array of size 
                                   zz->LB.Num_Global_Parts * wgtflag 
                                   containing the percentage of work 
                                   to be assigned to each part.    */
)
{
  char    yo[] = "rcb_fn";
  int     proc,nprocs;              /* my proc id, total # of procs */
  struct Dot_Struct *dotpt;         /* temporary pointer to rcb->Dots. */
  int     pdotnum;                  /* # of dots - decomposition changes it */
  int    *dotmark = NULL;           /* which side of median for each dot */
  int     dotnum;                   /* number of dots */
  int     dotmax = 0;               /* max # of dots arrays can hold */
  int     dottop;                   /* dots >= this index are new */
  int     proclower;                /* 1st proc in lower set */
  int     procmid;                  /* 1st proc in upper set */
  int     partlower;                /* 1st part in lower set */
  int     partmid;                  /* 1st part in upper set */
  int     set;                      /* which set processor is in = 0/1 */
  int     old_set;                  /* set processor was in last cut = 0/1 */
  int     root;                     /* part that stores last cut */
  int     num_procs;                /* number of procs in current set */
  int     num_parts;                /* number of parts in current set */
  int     dim;                      /* which of 3 axes median cut is on */
  int     ierr = ZOLTAN_OK;         /* error flag. */
  int    *proc_list = NULL;         /* temp array for reusing old cuts */
  int     outgoing;                 /* number of outgoing dots for reuse */
  double *coord = NULL;             /* temp array for serial_rcb */
  double *wgts = NULL;              /* temp array for serial_rcb */
  double *pts ;
  double  valuehalf;                /* median cut position */
  int     first_guess;              /* flag if first guess for median search */
  int     allocflag;                /* have to re-allocate space */
  double  time1=0,time2=0,time3=0,time4=0;  /* timers */
  double  timestart=0,timestop=0;   /* timers */
  double  timers[4]={0.,0.,0.,0.};  /* diagnostic timers 
			              0 = start-up time before recursion
				      1 = time before median iterations
				      2 = time in median iterations
				      3 = communication time */
  ZOLTAN_GNO_TYPE counters[7];        /* diagnostic counts
			              0 = unused
				      1 = # of dots sent
				      2 = # of dots received
				      3 = most dots this proc ever owns
				      4 = most dot memory this proc ever allocs
				      5 = # of times a previous cut is re-used
				      6 = # of reallocs of dot array  */
  ZOLTAN_GNO_TYPE reuse_count[7];   /* counter (as above) for reuse to record
                                       the number of dots premoved */
  int     i,j;                      /* local variables */
  int     use_ids;                  /* When true, global and local IDs will be
                                       stored along with dots in the RCB_STRUCT.
                                       When false, storage, manipulation, and
                                       communication of IDs is avoided.     
                                       Set by call to Zoltan_RB_Use_IDs(). */
  RCB_STRUCT *rcb = NULL;           /* Pointer to data structures for RCB.  */
  struct rcb_box *rcbbox = NULL;    /* bounding box of final RCB sub-domain */
  struct rcb_tree *treept = NULL;   /* tree of RCB cuts - only single cut on 
                                       exit */

  double start_time=0, end_time=0;
  double lb_time[2]={0,0};
  int    tfs_disregard_results = 0; /* added for Tflops_Special; all procs
                                       must enter Zoltan_RB_find_bisector
                                       when Tflops_Special==1.  This flag
                                       indicates that the aspect ratio test
                                       failed on this processor, so its 
                                       results from Zoltan_RB_find_bisector
                                       should be ignored. */
  int    tfs[2], tmp_tfs[2];        /* added for Tflops_Special; max number
                                       of procs and parts over all processors
                                       in each iteration (while loop) of 
                                       parallel partitioning.  */
  int    old_nprocs;                /* added for Tflops_Special */
  int    old_nparts;                /* added for Tflops_Special */
  double newwgt;
  double weight[RB_MAX_WGTS];       /* weight for current set */
  double weightlo[RB_MAX_WGTS];     /* weight of lower side of cut */
  double weighthi[RB_MAX_WGTS];     /* weight of upper side of cut */
  double weightlo_best[RB_MAX_WGTS];     /* temp weightlo */
  double weighthi_best[RB_MAX_WGTS];     /* temp weighthi */
  double wgtscale[RB_MAX_WGTS];     /* weight scaling factors */
  double fraclo[RB_MAX_WGTS];       /* desired weight in lower half */
  int *dotlist = NULL;              /* list of dots used only in find_median;
                                       allocated above find_median for 
                                       better efficiency (don't necessarily
                                       have to realloc for each find_median).*/
  int lock_direction = 0;           /* flag to determine direction after first
                                       iteration */
  int one_cut_dir= 0;               /* try only one cut direction */
  int level;                        /* recursion level of RCB for preset_dir */
  int *dim_spec = NULL;             /* specified direction for preset_dir */
  int fp=0;                         /* first part assigned to this proc */
  int np=0;                         /* number of parts assigned to this proc */
  int wgtdim;                       /* max(wgtflag,1) */
  int breakflag;                    /* flag for exiting loop */
  int *dotmark0 = NULL;             /* temp dotmark array */
  int *dotmark_best = NULL;         /* temp dotmark array */
  double valuehalf_best= 0.0L;      /* temp valuehalf */
  int dim_best;                     /* best cut dimension  */
  double norm_max, norm_best;       /* norm of largest half after bisection */
  double max_box;                   /* largest length of bbox */
  char msg[128];                    /* buffer for error messages */
  double pt[3];
  int *dindx = NULL, *tmpdindx = NULL;

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;
  int free_comm = FALSE;            /* Flag indicating whether MPI_Comm_free
                                       should be called on local_comm at end. */
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
  num_parts = zz->LB.Num_Global_Parts;

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);
  MPI_Op_create(&Zoltan_RCB_box_merge,1,&box_op);

  /* initializations */

  if ((wgtflag > 1) && (pivot_choice == PIVOT_CHOICE_RANDOM)){
    /* If RANDOM turns out to be wanted for wgtflag>1, we can implement it */
    ZOLTAN_PRINT_WARN(proc, yo, 
      "random_pivots turned off because it is not implemented for multiple weights");
    pivot_choice = PIVOT_CHOICE_BISECTION;
  }

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
  ierr = Zoltan_RCB_Build_Structure(zz, &pdotnum, &dotmax, wgtflag, overalloc,
                                    use_ids, gen_tree);

  if (ierr < 0) {
    ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from Zoltan_RCB_Build_Structure.");
    goto End;
  }

  rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);

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

  /* Ensure there are dots */

  MPI_Allreduce(&dotnum, &i, 1, MPI_INT, MPI_MAX, zz->Communicator);

  if (i == 0){
    if (proc == 0){
      ZOLTAN_PRINT_WARN(proc, yo, "RCB partitioning called with no objects");
    }
    rcbbox->lo[0] = rcbbox->lo[1] = rcbbox->lo[2] = 0;
    rcbbox->hi[0] = rcbbox->hi[1] = rcbbox->hi[2] = 0;
    timestart = timestop = 0;
    goto EndReporting;
  }

  /* If using RCB for matching, need to generate candidate lists.
   * Candidate lists include input GIDs, LIDs as provided by the application.
   * We need to capture that input here before we move any dots!
   * We return it in the import lists.
   * Candidates will be computed after partitioning and returned in the
   * export lists.
   */
  if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
    ierr = Zoltan_RB_Candidates_Copy_Input(zz, dotnum,
                                           rcb->Global_IDs, rcb->Local_IDs,
                                           &rcb->Dots, 
                                           num_import,
                                           import_global_ids, import_local_ids,
                                           import_procs, import_to_part);
    if (ierr < 0) {
       ZOLTAN_PRINT_ERROR(proc,yo,
                          "Error returned from Zoltan_RB_Return_Arguments.");
       goto End;
    }
  }

  /* create mark and list arrays for dots */

  allocflag = 0;
  wgtdim = (wgtflag>0 ? wgtflag : 1);
  if (dotmax > 0) {
    if (!(dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))
     || !(dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }
  else {
    dotmark = NULL;
    dotlist = NULL;
  }

  /* Also, if this is not first time through, send dots to previous proc. */
  if (reuse) {
    dotpt  = &rcb->Dots; 
    pt[1] = pt[2] = 0;

    if (treept[0].dim != -1) {
      /* find previous location of dots */
      for (outgoing = i = 0; i < dotnum; i++) {
        pt[0] = dotpt->X[i];
        if (dotpt->Y){
          pt[1] = dotpt->Y[i];
          if (dotpt->Z){
            pt[2] = dotpt->Z[i];
          }
        }
        ierr = Zoltan_RB_Point_Assign(zz, pt, &dotmark[i], NULL);
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(proc, yo, 
                             "Error returned from Zoltan_RB_Point_Assign");
          goto End;
        }
        if (dotmark[i] != proc) outgoing++;
      }

      if (outgoing)
        if ((proc_list = (int *) ZOLTAN_MALLOC(outgoing*sizeof(int))) == NULL) {
          ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
          ierr = ZOLTAN_MEMERR;
          goto End;
        }

      for (dottop = j = i = 0; i < dotnum; i++)
        if (dotmark[i] != proc)
          proc_list[j++] = dotmark[i];
        else
          dottop++;

      /* move dots */
      allocflag = 0;
#if 0
      ierr = Zoltan_RB_Send_Dots(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                                 &(rcb->Dots), &dotmark, 
                                 proc_list, outgoing, 
                                 &dotnum, &dotmax, proc, &allocflag,
                                 overalloc, stats, reuse_count, use_ids,
                                 zz->Communicator);
#else
      ierr = Zoltan_RB_Send_Dots_less_memory(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                                 &(rcb->Dots), &dotmark, 
                                 proc_list, outgoing, 
                                 &dotnum, &dotmax, proc, &allocflag,
                                 overalloc, stats, reuse_count, use_ids,
                                 zz->Communicator);
#endif
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(proc, yo, 
                           "Error returned from Zoltan_RB_Send_Dots.");
        goto End;
      }

      /* update counters */
      if (dotnum > counters[3]) counters[3] = dotnum;
      if (dotmax > counters[4]) counters[4] = dotmax;
      counters[6] += reuse_count[6];

      if (outgoing) ZOLTAN_FREE(&proc_list);
    }
  }

  dotpt = &rcb->Dots;

  if (wgtflag == 0){   /* no user supplied weights */
    weightlo[0] = (double) dotnum;
  }
  else {
    /* put sum of weights in weightlo */
    double *wgt;
    for (j=0; j<dotpt->nWeights; j++){
      weightlo[j] = 0.0;
      wgt = dotpt->Weight + j;
      for (i=0; i < dotnum; i++){
        weightlo[j] += *wgt;
        wgt += wgtflag;
      }
    }
  }
  /* Let weight be the global sum of all weights. */
  MPI_Allreduce(weightlo, weight, wgtdim, MPI_DOUBLE, MPI_SUM, zz->Communicator);

  /* Set weight scaling factors. */

  wgtscale[0] = 1.0;
  if (wgtdim > 1){
    for (j=0; j<wgtdim; j++){
      if (obj_wgt_comp || (weight[j]==0.0))
        wgtscale[j] = 1.0;
      else{
        wgtscale[j] = 1.0/weight[j]; /* normalize to make sum 1.0 */
        weight[j] = 1.0;
      }
    }
  }

  if (wgtflag == 0){   /* no user supplied weights */
    dotpt->uniformWeight = wgtscale[0];
  }
  else{
    /* scale the weights */
    double *wgt;
    for (j=0; j<dotpt->nWeights; j++){
      wgt = dotpt->Weight + j;
      for (i=0; i < dotnum; i++){
        newwgt = *wgt * wgtscale[j];
        *wgt = newwgt;
        wgt += wgtflag;
      }
    }
  }

  if (check_geom) {
    ierr = Zoltan_RB_check_geom_input(zz, dotpt, dotnum);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo, 
                         "Error returned from Zoltan_RB_check_geom_input");
      goto End;
    }
  }

  /* initialize sub-domain bounding box to entire domain */
  compute_RCB_box(rcbbox, dotnum, dotpt, NULL, box_op, box_type, 
                  zz->Communicator, zz->Num_Proc, zz->Proc, zz->Proc, 
                  zz->Tflops_Special);

  /* if preset_dir is turned on, assign cut order according to order of 
     directions */
  if (preset_dir) {
     ierr = set_preset_dir(proc, num_parts, preset_dir, rcbbox, &dim_spec);
     if (ierr < 0) 
       goto End;
  }


  /* Also, if this is not first time through, set lock_direction to reuse
     the directions determined previously */
  if (reuse_dir) {
     if (treept[0].dim != -1)
        lock_direction = 1;
  }

  /* create local communicator for use in recursion */

  if (zz->Tflops_Special)
     local_comm = zz->Communicator;
  else {
     MPI_Comm_dup(zz->Communicator,&local_comm);
     free_comm = TRUE;
  }

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    time2 = Zoltan_Time(zz->Timer);
    timers[0] = time2 - time1;
  }

  /* Main loop: recursively halve until just one part or proc in set */
  
  old_nprocs = num_procs = nprocs;
  old_nparts = num_parts;
  partlower = 0;
  root = 0;
  old_set = 1;
  ierr = Zoltan_LB_Proc_To_Part(zz, proc, &np, &fp);
  if (treept) {
    for (i = fp; i < (fp + np); i++) {
      treept[i].parent = 0;
      treept[i].left_leaf = 0;
    }
  }
  if (zz->Tflops_Special) {
    proclower = 0;
    tfs[0] = nprocs;
    tfs[1] = num_parts;
  }
  level = 0;

  while ((num_parts > 1 && num_procs > 1) || 
         (zz->Tflops_Special && tfs[0] > 1 && tfs[1] > 1)) {

    sprintf(msg, "In main RCB loop: num_parts=%d, num_procs=%d\n", 
            num_parts, num_procs);
    ZOLTAN_TRACE_DETAIL(zz, yo, msg);

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time1 = Zoltan_Time(zz->Timer);

    ierr = Zoltan_Divide_Machine(zz, zz->Obj_Weight_Dim, part_sizes, 
                                 proc, local_comm, &set, 
                                 &proclower, &procmid, &num_procs, 
                                 &partlower, &partmid, &num_parts, 
                                 fraclo);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error in Zoltan_Divide_Machine.");
      goto End;
    }

    /* tfs[0] is max number of processors in all sets over all processors -
     * tfs[1] is max number of parts in all sets over all processors -
     * force all processors to go through all levels of parallel rcb */
    if (zz->Tflops_Special) {
      tmp_tfs[0] = num_procs;
      tmp_tfs[1] = num_parts;
      MPI_Allreduce(tmp_tfs, tfs, 2, MPI_INT, MPI_MAX, local_comm);
    }

    /* create mark array and active list for dots */

    if (allocflag) {
      allocflag = 0;
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&dotlist);
      if (!(dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))
       || !(dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }

    if (recompute_box) {
      /* Compute bounding box for partition set */
      compute_RCB_box(rcbbox, dotnum, dotpt, NULL, box_op, box_type,
                      local_comm, old_nprocs, proc - proclower, proc, 
                      zz->Tflops_Special);
    }

    /* Compute max box length. */
    max_box = 0.0;
    for (dim=0; dim<3; dim++)
      if (rcbbox->hi[dim] - rcbbox->lo[dim] > max_box) 
        max_box = rcbbox->hi[dim] - rcbbox->lo[dim];

    /* try all cut directions and pick best one. */
    breakflag= 0;
    dim_best = -1;
    norm_best = -1.;
    one_cut_dir = (wgtdim<=1) || lock_direction || preset_dir;
    if (!one_cut_dir){
      if (!(dotmark0 = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))
       || !(dotmark_best = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))){
        ZOLTAN_PRINT_ERROR(proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
      for (j=0; j<dotnum; j++)
        dotmark0[j] = dotmark[j];
    }

    for (dim=0; dim<rcb->Num_Dim; dim++){

      /* One cut direction only */
      if (one_cut_dir){
        dim = cut_dimension(lock_direction, treept, partmid, 
                        preset_dir, dim_spec, &level, rcbbox);
        breakflag= 1;
      }
      else {
        /* Do not cut along this dimension if box is too thin. */
        if (zz->Tflops_Special) tfs_disregard_results = 0;

        if (rcbbox->hi[dim] - rcbbox->lo[dim] < max_box/max_aspect_ratio) {
          if (zz->Tflops_Special)/* All procs must participate in find_bisector;
                                    compute cut but don't use the results. */
            tfs_disregard_results = 1;
          else 
            continue;
        }

        /* Restore original dotmark array. */
        for (j=0; j<dotnum; j++)
          dotmark[j] = dotmark0[j];
      }

      dotpt = &rcb->Dots;

      pts = ((dim==0) ? dotpt->X : ((dim==1) ? dotpt->Y : dotpt->Z));

      /* determine if there is a first guess to use */
      /* The test on old_nparts is for the TFLOPS_SPECIAL flag */
      if (old_nparts > 1 && reuse && dim == treept[partmid].dim) {
        if (stats) counters[5]++;
        valuehalf = treept[partmid].cut;
        first_guess = 1;
      }
      else first_guess = 0;
  
      if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
        time2 = Zoltan_Time(zz->Timer);
  
      if (wgtflag <= 1){
        if (pivot_choice == PIVOT_CHOICE_BISECTION){
          if (!Zoltan_RB_find_median(          
               zz->Tflops_Special, pts, dotpt->Weight, dotpt->uniformWeight, dotmark, dotnum, proc, 
               fraclo, local_comm, &valuehalf, first_guess,
               nprocs, old_nprocs, proclower, old_nparts, 
               wgtflag, rcbbox->lo[dim], rcbbox->hi[dim], 
               weight[0], weightlo, weighthi,
               dotlist, rectilinear_blocks, average_cuts)) {
            ZOLTAN_PRINT_ERROR(proc, yo,"Error returned from Zoltan_RB_find_median.");
            ierr = ZOLTAN_FATAL;
            goto End;
          }
        }
        else{
          if (!Zoltan_RB_find_median_randomized(
               zz->Tflops_Special, pts, dotpt->Weight, dotpt->uniformWeight, dotmark, dotnum, proc, 
               fraclo, local_comm, &valuehalf, first_guess,
               nprocs, old_nprocs, proclower, old_nparts, 
               wgtflag, rcbbox->lo[dim], rcbbox->hi[dim], 
               weight[0], weightlo, weighthi,
               dotlist, rectilinear_blocks, average_cuts)) {
            ZOLTAN_PRINT_ERROR(proc, yo,"Error returned from Zoltan_RB_find_median_randomized.");
            ierr = ZOLTAN_FATAL;
            goto End;
          }
        }
      }
      else { 
        if (Zoltan_RB_find_bisector(
               zz, zz->Tflops_Special, pts, dotpt->Weight, dotpt->uniformWeight, dotmark, dotnum, 
               wgtflag, mcnorm, fraclo, local_comm, 
               &valuehalf, first_guess,
               old_nprocs, proclower, old_nparts, 
               rcbbox->lo[dim], rcbbox->hi[dim], 
               weight, weightlo, weighthi, &norm_max,
               dotlist, rectilinear_blocks, average_cuts)
          != ZOLTAN_OK) {
          ZOLTAN_PRINT_ERROR(proc, yo,"Error returned from Zoltan_RB_find_bisector.");
          ierr = ZOLTAN_FATAL;
          goto End;
        }

        /* test for better balance */
        if ((!one_cut_dir) && 
            ((norm_best<0.) ||
              (!tfs_disregard_results && (norm_max < norm_best)))) {
          norm_best = norm_max; 
          dim_best = dim;
          for (j=0; j<wgtdim; j++){
            weightlo_best[j] = weightlo[j];
            weighthi_best[j] = weighthi[j];
          }
          for (j=0; j<dotnum; j++)
            dotmark_best[j] = dotmark[j];
          valuehalf_best = valuehalf;
        }
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] Debug: cut dim=%1d, norm_max=%f, dim_best=%1d, norm_best=%f, cut value=%f\n", 
            proc, dim, norm_max, dim_best, norm_best, valuehalf);
          if (wgtflag>1)
            printf("[%1d] Debug: weightlo=(%f,%f), weighthi=(%f,%f)\n",
              proc, weightlo[0], weightlo[1],  weighthi[0], weighthi[1]);
        }
      }
      if (breakflag) break; /* if one_cut_dir is true */
    }

    if (!one_cut_dir){
      /* We have tried all cut directions. Restore best result. */
      if (dim_best<0){
        ZOLTAN_PRINT_ERROR(proc, yo, 
                           "Zoltan could not find any valid RCB cut.");
        ierr = ZOLTAN_FATAL;
        goto End;
      }
      dim = dim_best;
      for (j=0; j<wgtdim; j++){
        weightlo[j] = weightlo_best[j];
        weighthi[j] = weighthi_best[j];
      }
      for (j=0; j<dotnum; j++)
        dotmark[j] = dotmark_best[j];
      valuehalf = valuehalf_best;
      /* free temp arrays */
      ZOLTAN_FREE(&dotmark0);
      ZOLTAN_FREE(&dotmark_best);

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] Debug: BEST dim=%1d, cut value=%f, norm_max=%f \n", 
          proc, dim, valuehalf, norm_max);
        if (wgtflag>1)
          printf("[%1d] Debug: BEST weightlo=(%f,%f), weighthi=(%f,%f)\n",
            proc, weightlo[0], weightlo[1],  weighthi[0], weighthi[1]);
      }
    }

    if (set)    /* set weight for current part */
      for (j=0; j<wgtdim; j++) weight[j] = weighthi[j];
    else
      for (j=0; j<wgtdim; j++) weight[j] = weightlo[j];

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time3 = Zoltan_Time(zz->Timer);

    /* store cut info in tree only if proc "owns" partmid */
    /* test of partmid > 0 prevents treept[0] being set when this cut is 
       only removing low-numbered processors (proclower to procmid-1) that
       have no parts in them from the processors remaining to 
       be partitioned. */
    if (treept && partmid > 0 && partmid == fp) {
      treept[partmid].dim = dim;
      treept[partmid].cut = valuehalf;
      treept[partmid].parent = old_set ? -(root+1) : root+1;
      /* The following two will get overwritten when the information is
         assembled if this is not a terminal cut */
      treept[partmid].left_leaf = -partlower;
      treept[partmid].right_leaf = -partmid;
    }
    if (old_nprocs > 1 && partmid > 0 && partmid != partlower + old_nparts) {  
      /* old_nprocs > 1 test: Don't reset these values if proc is in loop only 
       * because of other procs for Tflops_Special.
       * partmid > 0 test:  Don't reset these values if low-numbered processors
       * (proclower to procmid-1) have zero parts and this cut is removing
       * them from the processors remaining to be partitioned. 
       * partmid != partlower + old_nparts test:  Don't reset these values if
       * cut is removing high-numbered processors with zero parts from
       * the processors remaining to be partitioned.
       */
      old_set = set;
      root = partmid;
    }
    
    /* use cut to shrink RCB domain bounding box */

    if (old_nprocs > 1) {
      if (!set)
        rcbbox->hi[dim] = valuehalf;
      else
        rcbbox->lo[dim] = valuehalf;
    }

    allocflag = 0;

    ierr = Zoltan_RB_Send_Outgoing(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), &dotmark,
                               &dottop, &dotnum, &dotmax,
                               set, &allocflag, overalloc,
                               stats, counters, use_ids, 
                                local_comm, proclower,
                               old_nprocs, partlower, partmid);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo,
                         "Error returned from Zoltan_RB_Send_Outgoing.");
      goto End;
    }

    /* create new communicators */

    if (zz->Tflops_Special) {
       if (set) {
         proclower = procmid;
         partlower = partmid;
       }
       old_nprocs = num_procs;
       old_nparts = num_parts;
    }
    else {
       if (set) partlower = partmid;
       MPI_Comm_split(local_comm,set,proc,&tmp_comm);
       MPI_Comm_free(&local_comm);
       local_comm = tmp_comm;
       old_nprocs = num_procs;
       old_nparts = num_parts;
    }

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
      time4 = Zoltan_Time(zz->Timer);
      timers[1] += time2 - time1;
      timers[2] += time3 - time2;
      timers[3] += time4 - time3;
    }
  }

  /* have recursed all the way to a single processor sub-domain */

  /* Send dots to correct processors for their parts.  This is needed
     most notably when a processor has zero parts on it, but still has
     some dots after the parallel partitioning. */

  ierr = Zoltan_RB_Send_To_Part(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), &dotmark, &dottop,
                               &dotnum, &dotmax, &allocflag, overalloc,
                               stats, counters, use_ids);
  if (ierr < 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Error returned from Zoltan_RB_Send_To_Part");
    goto End;
  }

  /* All dots are now on the processors they will end up on; now generate
   * more partitions if needed. */

  if (num_parts > 1) {
    dindx = (int *) ZOLTAN_MALLOC(dotnum * 2 * sizeof(int));
    tmpdindx = dindx + dotnum;
    if (allocflag) {
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&dotlist);
      if (!(dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))
       || !(dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }

    coord = (double *) ZOLTAN_MALLOC(dotnum*sizeof(double));
    wgts = (double *) ZOLTAN_MALLOC(wgtflag*dotnum*sizeof(double));

    if (!coord || (wgtflag && !wgts)){
      ZOLTAN_PRINT_ERROR(proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < dotnum; i++)
      dindx[i] = i;

    ierr = serial_rcb(zz, &rcb->Dots, dotmark, dotlist, old_set, root,
               rcbbox, weight, dotnum, num_parts,
               &(dindx[0]), &(tmpdindx[0]), partlower, 
               proc, wgtflag, lock_direction, reuse, stats, gen_tree, 
               preset_dir, 
               rectilinear_blocks, obj_wgt_comp, mcnorm, 
               recompute_box,
               box_op, box_type, average_cuts, 
               counters, treept, dim_spec, level,
               coord, wgts, part_sizes, wgtscale, rcb->Num_Dim, pivot_choice, 
               max_aspect_ratio, timers);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from serial_rcb");
      goto End;
    }
  }

  end_time = Zoltan_Time(zz->Timer);
  lb_time[1] = end_time - start_time;

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    MPI_Barrier(zz->Communicator);
    timestop = time1 = Zoltan_Time(zz->Timer);
  }

  /* error checking and statistics */

  if (check_geom) {
    ierr = Zoltan_RB_check_geom_output(zz, &rcb->Dots, part_sizes, np, fp,
                                       dotnum, pdotnum, rcbbox);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo, 
                         "Error returned from Zoltan_RB_check_geom_output");
      goto End;
    }
  }

EndReporting:

  /* update calling routine parameters */
  
  start_time = Zoltan_Time(zz->Timer);

  pdotnum = dotnum;

  /* Perform remapping (if requested) */

  if (zz->LB.Remap_Flag) {
    ierr = Zoltan_RB_Remap(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                           &(rcb->Dots), &dotnum, &dotmax,
                &allocflag, overalloc, stats, counters, use_ids);
    /* Note:  dottop is no longer valid after remapping.  Remapping might
       destroy the nice local-followed-by-non-local ordering of the 
       dots array.  Do not use dottop after remapping. */
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_RB_Remap.");
      goto End;
    }
  }

  /*  build return arguments */

  if (zz->LB.Return_Lists != ZOLTAN_LB_NO_LISTS && 
      zz->LB.Return_Lists != ZOLTAN_LB_CANDIDATE_LISTS) {
    /* zz->LB.Return_Lists is true ==> use_ids is true */
    ierr = Zoltan_RB_Return_Arguments(zz, rcb->Global_IDs, rcb->Local_IDs, 
                                      &rcb->Dots, num_import,
                                      import_global_ids, import_local_ids,
                                      import_procs, import_to_part, 
                                      dotnum);
    if (ierr < 0) {
       ZOLTAN_PRINT_ERROR(proc,yo,
                          "Error returned from Zoltan_RB_Return_Arguments.");
       goto End;
    }
  }
  else if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
    /* Select a candidate for each part and return it in the export_GIDs. */
    ierr = Zoltan_RB_Candidates_Output(zz, dotnum, dindx,
                                       rcb->Global_IDs, rcb->Local_IDs,
                                       &rcb->Dots, 
                                       *num_import, *import_global_ids,
                                       num_export, export_global_ids);
    if (ierr < 0) {
       ZOLTAN_PRINT_ERROR(proc,yo,
                          "Error returned from Zoltan_RB_Return_Candidates.");
       goto End;
    }
  }
  ZOLTAN_FREE(&dindx);

  if (gen_tree) {
    int *displ, *recvcount;
    int sendcount;
    struct rcb_tree *treetmp = NULL;  /* temporary tree of RCB cuts for use in
                                         MPI_Allgatherv */

    treetmp = (struct rcb_tree *)
	      ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(struct rcb_tree));
    displ = (int *) ZOLTAN_MALLOC(2 * zz->Num_Proc * sizeof(int));
    if (!displ || !treetmp) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    recvcount = displ + zz->Num_Proc;

    ierr = Zoltan_RB_Tree_Gatherv(zz, sizeof(struct rcb_tree), &sendcount,
                                  recvcount, displ);

    /* 
     * Create copy of treept so that MPI_Allgatherv doesn't use same
     * memory for sending and receiving; removes valgrind warning.
     */
    for (i = 0; i < zz->LB.Num_Global_Parts; i++)
      treetmp[i] = treept[i];

    MPI_Allgatherv(&treetmp[fp], sendcount, MPI_BYTE, treept, recvcount, displ,
                   MPI_BYTE, zz->Communicator);
    treept[0].dim = 0;
    for (i = 1; i < zz->LB.Num_Global_Parts; i++)
      if (treept[i].parent > 0)
        treept[treept[i].parent - 1].left_leaf = i;
      else if (treept[i].parent < 0)
        treept[-treept[i].parent - 1].right_leaf = i;
    ZOLTAN_FREE(&displ);
    ZOLTAN_FREE(&treetmp);
  }

  if (treept && zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    print_rcb_tree(zz, np, fp, &(treept[fp]));

  end_time = Zoltan_Time(zz->Timer);
  lb_time[0] += (end_time - start_time);

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
    Zoltan_RB_stats(zz, timestop-timestart,&rcb->Dots,dotnum,part_sizes,
                    timers,counters,stats,reuse_count,rcbbox,reuse);

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
    Zoltan_RB_Print_All(zz, rcb->Global_IDs, &rcb->Dots, 
                    dotnum, *num_import,
                    *import_global_ids, *import_procs);
  }

End:

  /* Free memory allocated by the algorithm. */

  if (free_comm) MPI_Comm_free(&local_comm);
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);

  ZOLTAN_FREE(&dotmark);
  ZOLTAN_FREE(&coord);
  ZOLTAN_FREE(&wgts);
  ZOLTAN_FREE(&dotlist);
  if (preset_dir) ZOLTAN_FREE(&dim_spec);
  if (dotmark0) ZOLTAN_FREE(&dotmark0);
  if (dotmark_best) ZOLTAN_FREE(&dotmark_best);

  if (!reuse && !gen_tree && rcb && (rcb->Tran.Target_Dim < 0)) {
    /* Free all memory used. */
    Zoltan_RCB_Free_Structure(zz);
  }
  else if (rcb != NULL) {
    /* Free only Dots and IDs; keep other structures. */
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    Zoltan_Free_And_Reset_Dot_Structure(&rcb->Dots);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);  
}



/******************************************************************************/

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

/******************************************************************************/

static void print_rcb_tree(ZZ *zz, int np, int fp, struct rcb_tree *treept_arr)
{
int i;
struct rcb_tree *treept = treept_arr;

  Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
  for (i = fp; i < fp+np; i++) {
    printf("Proc %d, Part %d:  Tree Struct:\n", zz->Proc, i);
    printf("                   cut        = %e\n", treept->cut);
    printf("                   dim        = %d\n", treept->dim);
    printf("                   parent     = %d\n", treept->parent);
    printf("                   left_leaf  = %d\n", treept->left_leaf);
    printf("                   right_leaf = %d\n", treept->right_leaf);
    treept++;
  }
  Zoltan_Print_Sync_End(zz->Communicator, TRUE);
}

/******************************************************************************/

static int cut_dimension(
  int lock_direction,        /* flag indicating whether directions are locked */
  struct rcb_tree *treept,   /* tree of RCB cuts */
  int partmid,               /* lowest part in set 1. */
  int preset_dir,            /* Set order of directions:     0: don't set,
                                 1: xyz,        2: xzy,      3: yzx,
                                 4: yxz,        5: zxy,      6: zyx  */
  int *dim_spec,             /* specified direction for preset_dir */
  int *level,                /* recursion level of RCB for preset_dir */
  struct rcb_box *rcbbox     /* bounding box of RCB sub-domain */
)
{
/* compute dim = dimension (xyz = 012) to bisect on */
int dim;

  if (lock_direction)
    dim = treept[partmid].dim;
  else if (preset_dir)
    dim = dim_spec[(*level)++];
  else {
    dim = 0;
    if (rcbbox->hi[1] - rcbbox->lo[1] > rcbbox->hi[0] - rcbbox->lo[0])
      dim = 1;
    if (dim == 0 && 
        rcbbox->hi[2] - rcbbox->lo[2] > rcbbox->hi[0] - rcbbox->lo[0])
      dim = 2;
    if (dim == 1 &&
        rcbbox->hi[2] - rcbbox->lo[2] > rcbbox->hi[1] - rcbbox->lo[1])
      dim = 2;
  }

  return dim;
}

/******************************************************************************/

static int set_preset_dir(
  int proc,                  /* Current processor */
  int nparts,                /* total number of parts */
  int preset_dir,            /* Set order of directions:     0: don't set,
                                 1: xyz,        2: xzy,      3: yzx,
                                 4: yxz,        5: zxy,      6: zyx  */
  struct rcb_box *rcbbox,    /* bounding box of RCB sub-domain */
  int **dim_spec             /* specified direction for preset_dir */
)
{
/* if preset_dir is turned on, count number of levels of recursion,
   determine number of cuts in each direction, and then assign those
   cuts to an order according to order of directions */

char *yo = "set_preset_dir";
int ix[3];                        /* temporaries for preset_dir */
double wx, wy, wz;                /* width for preset_dir */
int tmp_nparts;
int level;
int i, j;
int ierr = ZOLTAN_OK;
int dim = 0;

  if (preset_dir < 0 || preset_dir > 6) {
    ZOLTAN_PRINT_WARN(proc, yo, 
                    "Parameter RCB_SET_DIRECTIONS out of bounds; reset to 1.");
    ierr = ZOLTAN_WARN;
    preset_dir = 1;
  }

  wx = rcbbox->hi[0] - rcbbox->lo[0];
  wy = rcbbox->hi[1] - rcbbox->lo[1];
  wz = rcbbox->hi[2] - rcbbox->lo[2];

  for (i = 0; i < 3; i++) ix[i] = 0;

  for (tmp_nparts = nparts, level = 0; tmp_nparts > 1; level++) {
    tmp_nparts = (tmp_nparts + 1)/2;
    if (wz > wx && wz > wy) {
      ix[2]++;
      wz /= 2.0;
    }
    else if (wy > wx && wy > wz) {
      ix[1]++;
      wy /= 2.0;
    }
    else {
      ix[0]++;
      wx /= 2.0;
    }
  }

  if (level > 0) {
    *dim_spec = (int *) ZOLTAN_MALLOC(level*sizeof(int));
    if (*dim_spec == NULL) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
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
      (*dim_spec)[i] = dim;
      ix[dim]--;
    }
  }

End:
  return ierr;
}

/******************************************************************************/

static int serial_rcb(
  ZZ *zz,
  struct Dot_Struct *dotpt,  /* pointer to rcb->Dots. */
  int *dotmark,              /* which side of median for each dot */
  int *dotlist,              /* list of dots used only in find_median;
                                allocated above find_median for 
                                better efficiency (don't necessarily
                                have to realloc for each find_median).*/
  int old_set,               /* Set the objects to be partitioned were in 
                                for last cut */
  int root,                  /* part for which last cut was stored. */
  struct rcb_box *rcbbox,    /* bounding box of RCB sub-domain */
  double *weight,            /* weight(s) for current set */
  int dotnum,                /* number of input dots */
  int num_parts,             /* number of parts to create. */
  int *dindx,                /* Index into dotpt for dotnum dots to be 
                                partitioned; reordered in serial_rcb so set0
                                dots are followed by set1 dots. */
  int *tmpdindx,             /* Temporary memory used in reordering dindx. */
  int partlower,             /* smallest part number to be created. */
  int proc,                  /* processor number. */
  int wgtflag,               /* No. of weights per dot. */
  int lock_direction,        /* flag indicating whether directions are locked */
  int reuse,                 /* (0) don't use (1) use previous cuts
                                stored in treept at initial guesses.      */
  int stats,                 /* Print timing & count summary?             */
  int gen_tree,              /* (0) do not (1) do generate full treept    */
  int preset_dir,            /* Set order of directions:     0: don't set,
                                 1: xyz,        2: xzy,      3: yzx,
                                 4: yxz,        5: zxy,      6: zyx  */
  int rectilinear_blocks,    /* (0) do (1) don't break ties in find_median*/
  int obj_wgt_comp,          /* (1) obj wgts have same units, no scaling */
  int mcnorm,                /* norm for multicriteria bisection */
  int recompute_box,         /* (0) do (1) don't recompute bounding box for
                                partition sets at each level of recursion */
  MPI_Op box_op,             /* Operator needed if recompute_box */
  MPI_Datatype box_type,     /* Data type needed if recompute_box */
  int average_cuts,          /* Flag forcing median line to be drawn halfway
                                between two closest objects. */
  ZOLTAN_GNO_TYPE counters[],            /* diagnostic counts */
  struct rcb_tree *treept,   /* tree of RCB cuts */
  int *dim_spec,             /* specified direction for preset_dir */
  int level,                 /* recursion level of RCB for preset_dir */
  double *coord,
  double *wgts,
  float *part_sizes,         /* Array of size zz->LB.Num_Global_Parts
                                containing the percentage of work to be
                                assigned to each part.               */
  double *wgtscale,          /* Array of size wgtflag that gives the
                                scaling factors for each weight dimension. */
  int ndim,                  /* number of geometric dimensions */
  int pivot_choice, 
  double max_aspect_ratio,
  double timers[]            /* as in rcb_fn */
)
{
  char *yo = "serial_rcb";
  int ierr = ZOLTAN_OK;
  int i, j;
  int dim;
  int first_guess;
  int new_nparts;
  int partmid;
  int wgtdim;
  double valuehalf;
  double fractionlo[RB_MAX_WGTS];    /* desired weight in lower half */
  double weightlo[RB_MAX_WGTS];      /* weight in lower half */
  double weighthi[RB_MAX_WGTS];      /* weight in upper half */
  double weightlo_best[RB_MAX_WGTS]; /* temp weightlo */
  double weighthi_best[RB_MAX_WGTS]; /* temp weighthi */
  int set0, set1, breakflag, one_cut_dir;
  struct rcb_box tmpbox;
  int *dotmark0 = NULL;             /* temp dotmark array */
  int *dotmark_best = NULL;         /* temp dotmark array */
  double valuehalf_best = 0.0L;     /* temp valuehalf */
  int dim_best;                     /* best cut dimension  */
  double norm_max, norm_best;       /* norm of largest half after bisection */
  double max_box;                   /* largest length of bbox */
  double *c, *w;
  double uniformWeight;
  double start_time=0., end_time;

  wgtdim = (wgtflag>0 ? wgtflag : 1);

  if (num_parts == 1) {
    for (i = 0; i < dotnum; i++)
      dotpt->Part[dindx[i]] = partlower;
  }
  else {
    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      start_time = Zoltan_Time(zz->Timer);

    ierr = Zoltan_Divide_Parts(zz, zz->Obj_Weight_Dim, part_sizes, num_parts,
                               &partlower, &partmid, fractionlo);

    /* Compute max box length. */
    max_box = 0.0;
    for (dim=0; dim<3; dim++)
      if (rcbbox->hi[dim] - rcbbox->lo[dim] > max_box) 
        max_box = rcbbox->hi[dim] - rcbbox->lo[dim];

    /* try all cut directions and pick best one. */
    breakflag= 0;
    dim_best = -1;
    norm_best = -1.;
    one_cut_dir = (wgtflag<=1) || lock_direction || preset_dir;
    if (!one_cut_dir){
      if (!(dotmark0 = (int *) ZOLTAN_MALLOC(dotnum*sizeof(int)))
       || !(dotmark_best = (int *) ZOLTAN_MALLOC(dotnum*sizeof(int)))){
        ZOLTAN_PRINT_ERROR(proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
      for (j=0; j<dotnum; j++)
        dotmark0[j] = dotmark[j];
    }

    if (recompute_box) {
      /* Compute bounding box for partition set */
      compute_RCB_box(rcbbox, dotnum, dotpt, dindx, 
                      box_op, box_type, MPI_COMM_SELF, 1, 0, proc, 0);
    }

    for (dim=0; dim<ndim; dim++){

      /* One cut direction only */
      if (one_cut_dir){
        dim = cut_dimension(lock_direction, treept, partmid, 
                            preset_dir, dim_spec, &level, rcbbox);
        breakflag= 1;
      }
      else {
        /* Do not cut along this dimension if box is too thin. */

        if (rcbbox->hi[dim] - rcbbox->lo[dim] < max_box/max_aspect_ratio) 
          continue;

        /* Restore original dotmark array. */
        for (j=0; j<dotnum; j++)
          dotmark[j] = dotmark0[j];
      }

      c = ((dim==0) ? dotpt->X : ((dim==1) ? dotpt->Y : dotpt->Z));
      uniformWeight = dotpt->uniformWeight;

      /* copy correct coordinate value into the temporary array */
      /* dotpt = rcb->Dots; */

      for (i = 0; i < dotnum; i++) {
        coord[i] = c[dindx[i]];        /* INVALID WRITE */
      }

      if (dotpt->nWeights > 0){
        w = wgts;
        for (i = 0; i < dotnum; i++) {
          for (j=0; j<dotpt->nWeights; j++){
            *w++ = dotpt->Weight[dindx[i] * dotpt->nWeights + j];
          }
        }
      }

      /* determine if there is a first guess to use */
      if (reuse && dim == treept[partmid].dim) {
        if (stats) counters[5]++;
        valuehalf = treept[partmid].cut;
        first_guess = 1;
      }
      else first_guess = 0;

      if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
        end_time = Zoltan_Time(zz->Timer);
        timers[1] += end_time - start_time;
        start_time = end_time;
      }
  
      if (wgtflag <= 1){
        /* Call find_median with Tflops_Special == 0; avoids communication */
        if (pivot_choice == PIVOT_CHOICE_BISECTION){
          if (!Zoltan_RB_find_median(
                 0, coord, wgts, uniformWeight, dotmark, dotnum, proc, 
                 fractionlo, MPI_COMM_SELF, &valuehalf, 
                 first_guess, zz->Num_Proc, 1, zz->Proc, num_parts,
                 wgtflag, rcbbox->lo[dim], rcbbox->hi[dim], 
                 weight[0], weightlo, weighthi,
                 dotlist, rectilinear_blocks, average_cuts)) {
            ZOLTAN_PRINT_ERROR(proc, yo,"Error returned from Zoltan_RB_find_median.");
            ierr = ZOLTAN_FATAL;
            goto End;
          }
        }
        else{  
          if (!Zoltan_RB_find_median_randomized(
                 0, coord, wgts, uniformWeight,  dotmark, dotnum, proc, 
                 fractionlo, MPI_COMM_SELF, &valuehalf, 
                 first_guess,
                 zz->Num_Proc, 1, zz->Proc, num_parts,
                 wgtflag, rcbbox->lo[dim], rcbbox->hi[dim], 
                 weight[0], weightlo, weighthi,
                 dotlist, rectilinear_blocks, average_cuts)) {
            ZOLTAN_PRINT_ERROR(proc, yo,"Error returned from Zoltan_RB_find_median.");
            ierr = ZOLTAN_FATAL;
            goto End;
          }
        }
      }
      else { 
        /* Call find_bisector with Tflops_Special == 0; avoids communication */
        if (Zoltan_RB_find_bisector(
             zz, 0, coord, wgts, uniformWeight, dotmark, dotnum, 
             wgtflag, mcnorm, fractionlo, MPI_COMM_SELF, 
             &valuehalf, first_guess, 
             1, zz->Proc, num_parts, 
             rcbbox->lo[dim], rcbbox->hi[dim], 
             weight, weightlo, weighthi, &norm_max,
             dotlist, rectilinear_blocks, average_cuts)
          != ZOLTAN_OK) {
          ZOLTAN_PRINT_ERROR(proc, yo,"Error returned from Zoltan_RB_find_bisector.");
          ierr = ZOLTAN_FATAL;
          goto End;
        }

        /* test for better balance */
        if ((!one_cut_dir) && 
            ((norm_best<0.) ||
             (norm_max < norm_best))) {
          norm_best = norm_max; 
          dim_best = dim;
          for (j=0; j<wgtdim; j++){
            weightlo_best[j] = weightlo[j];
            weighthi_best[j] = weighthi[j];
          }
          for (j=0; j<dotnum; j++)
            dotmark_best[j] = dotmark[j];
          valuehalf_best = valuehalf;
        }
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] Debug: cut dim=%1d, norm_max=%f, dim_best=%1d, norm_best=%f, cut value=%f\n", 
            proc, dim, norm_max, dim_best, norm_best, valuehalf);
          if (wgtflag>1)
            printf("[%1d] Debug: weightlo=(%f,%f), weighthi=(%f,%f)\n",
              proc, weightlo[0], weightlo[1],  weighthi[0], weighthi[1]);
        }
      }
      if (breakflag) break; /* if one_cut_dir is true */
    }

    if (!one_cut_dir){
      /* We have tried all cut directions. Restore best result. */
      if (dim_best<0){
        ZOLTAN_PRINT_ERROR(proc, yo, 
                           "Zoltan could not find any valid RCB cut.");
        ierr = ZOLTAN_FATAL;
        goto End;
      }
      dim = dim_best;
      for (j=0; j<wgtdim; j++){
        weightlo[j] = weightlo_best[j];
        weighthi[j] = weighthi_best[j];
      }
      for (j=0; j<dotnum; j++)
        dotmark[j] = dotmark_best[j];
      valuehalf = valuehalf_best;
      /* free temp arrays */
      ZOLTAN_FREE(&dotmark0);
      ZOLTAN_FREE(&dotmark_best);

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] Debug: BEST dim=%1d, cut value=%f, norm_max=%f \n", 
          proc, dim, valuehalf, norm_max);
        if (wgtflag>1)
          printf("[%1d] Debug: BEST weightlo=(%f,%f), weighthi=(%f,%f)\n",
            proc, weightlo[0], weightlo[1],  weighthi[0], weighthi[1]);
      }
    }

    /* By now, we have the right values for the best cut. */
    if (gen_tree) {
      treept[partmid].dim = dim;
      treept[partmid].cut = valuehalf;
      treept[partmid].parent = old_set ? -(root+1) : root+1;
      /* The following two will get overwritten when the information is
         assembled if this is not a terminal cut */
      treept[partmid].left_leaf = -partlower;
      treept[partmid].right_leaf = -partmid;
    }

    root = partmid;

    /* Create new dindx, grouping set 0 and set 1 dots together */
    for (set0 = 0, set1 = dotnum, i = 0; i < dotnum; i++) {
      if (dotmark[i] == 0) 
        tmpdindx[set0++] = dindx[i];
      else
        tmpdindx[--set1] = dindx[i];
    }
    memcpy(dindx, tmpdindx, dotnum * sizeof(int));

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
      end_time = Zoltan_Time(zz->Timer);
      timers[2] += end_time - start_time;
    }

    /* If set 0 has at least one part and at least one dot,
     * call serial_rcb for set 0 */
    new_nparts = partmid - partlower;
    if (new_nparts > 0 && set1 != 0) {
      memcpy(&tmpbox, rcbbox, sizeof(struct rcb_box));
      tmpbox.hi[dim] = valuehalf;
      ierr = serial_rcb(zz, dotpt, dotmark, dotlist, 0, root,
                        &tmpbox, weightlo, set0, new_nparts,
                        &(dindx[0]), &(tmpdindx[0]), partlower, 
                        proc, wgtflag, lock_direction, reuse, stats, gen_tree, 
                        preset_dir, 
                        rectilinear_blocks, obj_wgt_comp, mcnorm, 
                        recompute_box,
                        box_op, box_type, average_cuts, 
                        counters, treept, dim_spec, level,
                        coord, wgts, part_sizes, wgtscale, ndim, pivot_choice,
                        max_aspect_ratio, timers);
      if (ierr < 0) {
        goto End;
      }
    }

    /* If set 1 has at least one part and at least one dot,
     * call serial_rcb for set 1 */
    new_nparts = partlower + num_parts - partmid;
    if (new_nparts > 0 && set0 != dotnum) {
      memcpy(&tmpbox, rcbbox, sizeof(struct rcb_box));
      tmpbox.lo[dim] = valuehalf;
      ierr = serial_rcb(zz, dotpt, dotmark, dotlist, 1, root,
                        &tmpbox, weighthi, dotnum-set0, new_nparts,
                        &(dindx[set1]), &(tmpdindx[set1]), partmid, 
                        proc, wgtflag, lock_direction, reuse, stats, gen_tree, 
                        preset_dir, rectilinear_blocks, obj_wgt_comp, mcnorm,
                        recompute_box,
                        box_op, box_type, average_cuts,
                        counters, treept, dim_spec, level,
                        coord, wgts, part_sizes, wgtscale, ndim, pivot_choice,
                        max_aspect_ratio, timers);
      if (ierr < 0) {
        goto End;
      }
    }
  }
End:

  return ierr;
}

/*****************************************************************************/

static void compute_RCB_box(
  struct rcb_box *rcbbox,    /* bounding box of RCB sub-domain */
  int dotnum,                /* # of dots on this processor */
  struct Dot_Struct *dotpt,  /* temporary pointer to rcb->Dots. */
  int *dindx,                /* Index into dotpt used by serial RCB; 
                                NULL for parallel rcb. */
  MPI_Op box_op,    
  MPI_Datatype box_type,
  MPI_Comm comm,
  int nproc,                 /* number of processors */
  int rank,                  /* rank */
  int proc,                  /* global processor number */
  int Tflops_Special
)
{
/* Compute sub-domain bounding box to entire domain within communicator */
int i, k;
struct rcb_box boxtmp;

  rcbbox->lo[0] = rcbbox->lo[1] = rcbbox->lo[2] = DBL_MAX;
  rcbbox->hi[0] = rcbbox->hi[1] = rcbbox->hi[2] = -DBL_MAX;
  boxtmp.lo[0] = boxtmp.lo[1] = boxtmp.lo[2] = DBL_MAX;
  boxtmp.hi[0] = boxtmp.hi[1] = boxtmp.hi[2] = -DBL_MAX;

  if (dotpt->Z ==  NULL){
    boxtmp.hi[2] = boxtmp.lo[2] = 0;
    if (dotpt->Y ==  NULL){
      boxtmp.hi[1] = boxtmp.lo[1] = 0;
    }
  }
  
  for (i = 0; i < dotnum; i++) {
    k = (dindx ? dindx[i] : i);

    if (dotpt->X[k] < boxtmp.lo[0])
      boxtmp.lo[0] = dotpt->X[k];
    if (dotpt->X[k] > boxtmp.hi[0])
      boxtmp.hi[0] = dotpt->X[k];

    if (dotpt->Y){
      if (dotpt->Y[k] < boxtmp.lo[1])
        boxtmp.lo[1] = dotpt->Y[k];
      if (dotpt->Y[k] > boxtmp.hi[1])
        boxtmp.hi[1] = dotpt->Y[k];
      if (dotpt->Z){
        if (dotpt->Z[k] < boxtmp.lo[2])
          boxtmp.lo[2] = dotpt->Z[k];
        if (dotpt->Z[k] > boxtmp.hi[2])
          boxtmp.hi[2] = dotpt->Z[k];
      }
    }
  }

  if (!Tflops_Special)
    MPI_Allreduce(&boxtmp,rcbbox,1,box_type,box_op,comm);
  else {
    i = 1;
    Zoltan_RB_reduce(nproc, rank, proc, (void *) &boxtmp, (void *) rcbbox,
                     sizeof(struct rcb_box), &i, box_type, comm,
                     Zoltan_RCB_box_merge);
  }
}
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
