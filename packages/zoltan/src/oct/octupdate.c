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

#include "zz_const.h"
#include "octree_const.h"
#include "costs_const.h"
#include "oct_util_const.h"
#include "dfs_const.h"
#include "octupdate.h"
#include "octupdate_const.h"
#include "migreg_const.h"
#include "migoct_const.h"
#include "migtags_const.h"
#include "params_const.h"
#include <float.h>
#define POW(a,b) pow((double)(a),(double)(b))

/*test*/
/*extern void getMaxBounds(void *, double *, double *); */

/**************************** PROTOTYPES ************************************/
static int Zoltan_Oct_nUniqueRegions(OCT_Global_Info *OCT_info, pOctant oct);
static int Zoltan_Oct_CompareCoords(int dim, COORD pt1, COORD pt2);
static void initialize_region(ZZ *, pRegion *, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
			      int, float, int, double *);
static int lb_oct_init(ZZ *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **,
		       int **, int, int, int, int, int, int, float *); 
static void Zoltan_Oct_gen_tree_from_input_data(ZZ *zz, int, int *c1, int *c2,
						int *c3, float *c0, 
						int createpartree);
static pOctant Zoltan_Oct_global_insert(ZZ *, pRegion region);
static void Zoltan_Oct_terminal_refine(ZZ *, pOctant oct,int count);



/****************************************************************************/
/* NOTE: be careful later about region lists for nonterminal octants */

static int oct_nref=0;                             /* number of refinements */
static int oct_ncoarse=0;                          /* number of coarsenings */
#ifdef LGG_MIGOCT
static int IDcount = 0;                            /* renumbering of octs */
#endif
#define MAXOCTREGIONS_DEF 40
#define MINOCTREGIONS_DEF 10
static int MAXOCTREGIONS = MAXOCTREGIONS_DEF;
static int MINOCTREGIONS = MINOCTREGIONS_DEF;  /* min # of regions per oct */


/****************************************************************************/
/* parameters for the octpart method.  Used in  */
/* Zoltan_Oct_Set_Param and Zoltan_Octpart          */
static PARAM_VARS OCT_params[] = {
  { "OCT_DIM",          NULL, "INT", 0 },
  { "OCT_METHOD",       NULL, "INT", 0 },
  { "OCT_MAXOBJECTS",   NULL, "INT", 0 },
  { "OCT_MINOBJECTS",   NULL, "INT", 0 },
  { "OCT_OUTPUT_LEVEL", NULL, "INT", 0 },
  {  NULL,              NULL,  NULL, 0 }};

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Oct_Set_Param(
  char *name,                 /* name of variable */
  char *val                   /* value of variable */
)
{
int status;
PARAM_UTYPE result;           /* value returned from Zoltan_Check_Param */
int index;                    /* index returned from Zoltan_Check_Param */

  status = Zoltan_Check_Param(name, val, OCT_params, &result, &index);
  return(status);
}

/****************************************************************************/
int Zoltan_Octpart(
  ZZ *zz,                       /* The Zoltan structure with info for
                                   the OCTPART balancer.                    */
  float *part_sizes,            /* Input:  Array of size zz->Num_Global_Parts
                                   containing the percentage of work to be
                                   assigned to each partition.              */
  int *num_import,              /* Returned value: Number of non-local 
                                   objects assigned to this
                                   processor in the new decomposition.      */
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value: array of global IDs for
                                   non-local objects in this processor's new
                                   decomposition.                           */
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                   non-local objects in this processor's new
                                   decomposition.                           */
  int **import_procs,           /* Returned value:  array of processor IDs for
                                   processors owning the non-local objects in
                                   this processor's new decomposition.      */
  int **import_to_part,         /* Returned value:  array of partitions to
                                   which imported objects should be assigned.
                                   KDDKDD  Currently unused.  */
  int *num_export,              /* Not computed; return -1. */
  ZOLTAN_ID_PTR *export_global_ids, /* Not computed. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Not computed. */
  int **export_procs,           /* Not computed. */
  int **export_to_part          /* Not computed. */
) 
{
char *yo = "Zoltan_Octpart";
int oct_dim = 3;              /* Dimension of method to be used (2D or 3D)  */
int oct_method = 2;           /* Flag specifying curve to be used.          */
int oct_maxoctregions=MAXOCTREGIONS_DEF; /* max # of objs in leaves         */
int oct_minoctregions=MINOCTREGIONS_DEF; /* min # of objs in leaves         */
int oct_output_level = 0;     /* Flag specifying amount of output.          */
int oct_wgtflag = 0;          /* Flag specifying use of object weights.     */
int error = FALSE;            /* error flag                                 */


  Zoltan_Bind_Param(OCT_params, "OCT_DIM", (void *)&oct_dim);
  Zoltan_Bind_Param(OCT_params, "OCT_METHOD", (void *)&oct_method);
  Zoltan_Bind_Param(OCT_params, "OCT_MAXOBJECTS", (void *)&oct_maxoctregions);
  Zoltan_Bind_Param(OCT_params, "OCT_MINOBJECTS", (void *)&oct_minoctregions);
  Zoltan_Bind_Param(OCT_params, "OCT_OUTPUT_LEVEL",(void *)&oct_output_level);

  Zoltan_Assign_Param_Vals(zz->Params, OCT_params, zz->Debug_Level, zz->Proc, 
			   zz->Debug_Proc);

  /* Set oct_wgtflag based on the "key" parameter Obj_Weight_Dim */
  oct_wgtflag = (zz->Obj_Weight_Dim > 0);

  /* Initialization in case of early exit */
  *num_import = -1;  /* We don't compute any import data */
  *num_export = -1;

  /* Error checking for unimplemented features */
  if (zz->LB.PartDist) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
            "# partitions != # processors not yet implemented in OCTPART.  "
            "Try a different LB_METHOD.");
    error = TRUE;
  }
  /*if (!zz->LB.Uniform_Parts) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
            "Non-uniform partition sizes not yet implemented in OCTPART.  "
            "Try a different LB_METHOD.");
    error = TRUE;
  }*/

  /* Error checking for parameters */
  if (oct_dim < 2 || oct_dim > 3) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "OCT_DIM must be 2 or 3. ");
    error = TRUE;
  }
  if (oct_method < 0 || oct_method > 2) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "OCT_METHOD must be 0, 1, or 2");
    error = TRUE;
  }
  if (oct_maxoctregions < 1) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "OCT_MAXOBJECTS must be greater than 0");
    error = TRUE;
  }
  if (oct_minoctregions < 1) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "OCT_MINOBJECTS must be greater than 0");
    error = TRUE;
  }
  if (oct_minoctregions > oct_maxoctregions) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "OCT_MINOBJECTS must be less than "
		       "OCT_MAXOBJECTS");
    error = TRUE;
  }
  if (oct_output_level < 0 || oct_output_level > 3) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,"OCT_OUTPUT_LEVEL must be 0, 1, 2, or 3");
    error = TRUE;
  }

  if (error)
    return(ZOLTAN_FATAL);
  else {
    error = lb_oct_init(zz, num_export, export_global_ids, export_local_ids, 
			export_procs, export_to_part, 
			oct_dim, oct_method, oct_maxoctregions, 
			oct_minoctregions, oct_output_level, oct_wgtflag,
			part_sizes);
    return(error);
  }
}

/****************************************************************************/
/*
 * void lb_oct_init();
 *
 * initialize the calls needed to start the octree load balancing rounties
 */
static int lb_oct_init(
  ZZ *zz,                       /* The Zoltan structure with info for
                                   the OCTPART balancer.                    */
  int *num_export,              /* Number of non-local objs assigned to this
                                   processor in the new decomposition.      */
  ZOLTAN_ID_PTR *export_global_ids, /* Returned value: array of global IDs for
                                   non-local objects in this processor's new
                                   decomposition.                           */
  ZOLTAN_ID_PTR *export_local_ids,  /* Returned value:  array of local IDs for
                                   non-local objects in this processor's new
                                   decomposition.                           */
  int **export_procs,           /* Returned value:  array of processor IDs for
                                   processors owning the non-local objects in
                                   this processor's new decomposition.      */
  int **export_to_part,         /* Returned value:  array of partitions to 
                                   which objects are imported.
                                   KDDKDD Assume #parts==#procs.            */
  int oct_dim,                  /* Dimension of method (2D or 3D)           */
  int oct_method,               /* Flag specifying curve to be used.        */
  int oct_maxoctregions,        /* max # of objects in leaves of octree.    */
  int oct_minoctregions,        /* min # of objects in leaves of octree.    */
  int oct_output_level,         /* Flag specifying amount of output.        */
  int oct_wgtflag,              /* Flag specifying use of object weights.   */
  float *part_sizes             /* Array of size zz->Num_Global_Parts
                                   containing the percentage of work to be
                                   assigned to each partition.              */
) 
{
  char *yo = "lb_oct_init";
  OCT_Global_Info *OCT_info;
  int nsentags;                    /* number of tags being sent */
  pRegion export_regs;             /* */
  int nrectags;                    /* number of tags received */
  int kk;
  double time1,time2;              /* timers */
  double timestart,timestop;       /* timers */
  double timers[4];                /* diagnostic timers 
			              0 = start-up time before recursion
				      1 = time before median iterations
				      2 = time in median iterations
				      3 = communication time */
  int    counters[6];              /* diagnostic counts
			              0 = # of median iterations
				      1 = # of objects sent
				      2 = # of objects received
				      3 = most objects this proc ever owns
				      */
  float  c[4];
  int createpartree = 0;
  /*int num_gid_entries = zz->Num_GID;*/
  /*int num_lid_entries = zz->Num_LID;*/
  
  ZOLTAN_TRACE_ENTER(zz, yo);

  MPI_Barrier(zz->Communicator);
  timestart = MPI_Wtime();

  /* initialize timers and counters */
  counters[0] = 0;
  counters[1] = 0;
  counters[2] = 0;
  counters[3] = 0;
  counters[4] = 0;
  counters[5] = 0;
  c[0] = 0;
  c[1] = 0;
  c[2] = 0;
  c[3] = 0;
  timers[1] = 0.0;
  timers[2] = 0.0;
  timers[3] = 0.0;

  nsentags = nrectags = 0;

  if(zz->LB.Data_Structure == NULL) {
    OCT_info = Zoltan_Oct_POct_init(zz, zz->Proc, oct_dim);
    Zoltan_Oct_set_method(OCT_info, oct_method);
    Zoltan_Oct_set_maxregions(oct_maxoctregions);
    Zoltan_Oct_set_minregions(oct_minoctregions);
    createpartree = 1;
  }
  else {
    OCT_info = (OCT_Global_Info *) (zz->LB.Data_Structure);
  }

  /* create the octree structure */
  time1 = MPI_Wtime();

  ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_gen_tree_from_input_data");
  Zoltan_Oct_gen_tree_from_input_data(zz, oct_wgtflag, &counters[1],
				      &counters[2], &counters[3], &c[0], 
				      createpartree);

  time2 = MPI_Wtime();
  timers[0] = time2 - time1;                 /* time took to create octree */
  /* Zoltan_Oct_POct_printResults(OCT_info); */
  /* partition the octree structure */
  time1 = MPI_Wtime();
  ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_dfs_partition");
  /* old call to dfs_paritition: */ 
#if 0
  Zoltan_Oct_dfs_partition(zz, &counters[0], &c[1]);
#else
  /***************************
  if(zz->Proc == 0) {
    int debug_i;
    for(debug_i=0; debug_i<zz->Num_Proc; debug_i++) {
      fprintf(stdout,"Part_size[%d] = %f\n", debug_i, part_sizes[debug_i]);
    }
  }
  ****************************/
  Zoltan_Oct_dfs_partition(zz, &counters[0], &c[1], part_sizes);
#endif
  time2 = MPI_Wtime();
  timers[1] = time2 - time1;              /* time took to partition octree */

  if (oct_output_level > 2) {
    Zoltan_Oct_Plots(zz);
  }

  /* set up tags for migrations */
  time1 = MPI_Wtime();

#if 0  /* KDDKDD -- Count is never used; why is it computed? */
  {
  pRList  RootList;               /* list of all local roots */
  pOctant RootOct;                /* root octree octant */
  int count = 0; 
  RootList = Zoltan_Oct_POct_localroots(OCT_info);
  while((RootOct = RL_nextRootOctant(&RootList))) {
    while(RootOct) {
      if(Zoltan_Oct_isTerminal(RootOct)) {	
	count += Zoltan_Oct_nRegions(RootOct);
      }
      RootOct = Zoltan_Oct_POct_nextDfs(OCT_info, RootOct);
    }
  }
  }
#endif

  ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_dfs_migrate");
  Zoltan_Oct_dfs_migrate(zz, &nsentags, &export_regs, &nrectags, 
	         &c[2], &c[3], &counters[3], &counters[5]);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_fix_tags");
  if (zz->LB.Return_Lists) {
    *num_export = nrectags;
    if (nrectags > 0)
      Zoltan_Oct_fix_tags(zz, export_global_ids, export_local_ids, 
			  export_procs, export_to_part, nrectags,
			  export_regs);
  }

  time2 = MPI_Wtime();
  timers[2] = time2 - time1;               /* time took to setup migration */


#if 0  /* KDDKDD -- Count is never used; why is it computed? */
  {
  /* count the number of objects on this processor */
  pRList  RootList;               /* list of all local roots */
  pOctant RootOct;                /* root octree octant */
  int count = 0; 
  RootList = Zoltan_Oct_POct_localroots(OCT_info);
  while((RootOct = RL_nextRootOctant(&RootList))) {
    while(RootOct) {
      if(Zoltan_Oct_isTerminal(RootOct)) {	
	count += Zoltan_Oct_nRegions(RootOct);
      }
      RootOct = Zoltan_Oct_POct_nextDfs(OCT_info, RootOct);
    }
  }
  }
#endif

  counters[4] = nsentags;
  MPI_Barrier(zz->Communicator);
  timestop = MPI_Wtime();

  if (oct_output_level > 0) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_print_stats");
    Zoltan_Oct_print_stats(zz, timestop-timestart, timers, counters, c, 
                       oct_output_level);
  }

  for (kk = 0; kk < nrectags; kk++) {
    ZOLTAN_FREE(&(export_regs[kk].Global_ID));
    ZOLTAN_FREE(&(export_regs[kk].Local_ID));
  }
  ZOLTAN_FREE(&export_regs);
  ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_global_clear");
  Zoltan_Oct_global_clear(OCT_info);
  /* KDDKDD Don't understand how re-used octree will work, especially without
   * KDDKDD the Zoltan_Oct_Bounds_Geom function.  For now, we'll delete everything;
   * KDDKDD we can move back to saving some of the tree later.
   */
  Zoltan_Oct_Free_Structure(zz);
  /* KDDKDD END */

  /* Temporary return value until error codes are fully implemented. */
  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ZOLTAN_OK);
}

/****************************************************************************/
/*
 * void Zoltan_Oct_gen_tree_from_input_data()
 *
 * This function will create a root node of on each processor which will
 * then be used to create an octree with regions associated with it. The
 * tree will then be balanced and the output used to balance "mesh regions"
 * on several processors.
 */
static void Zoltan_Oct_gen_tree_from_input_data(ZZ *zz, int oct_wgtflag,
						int *c1, int *c2, int *c3, 
						float *c0, int createpartree) 
{
  char *yo = "Zoltan_Oct_gen_tree_from_input_data";
  pRList  RootList;       /* list of all local roots */
  pOctant RootOct;        /* root octree octant */
  COORD min,              /* min coord bounds of objects */
        max;              /* max coord bounds of objects */
  int num_extra;          /* number of orphaned objects */
  int num_objs;           /* total number of local objects */
  pRegion ptr,            /* pointer to iterate trough region list */
          ptr1;           /* pointer to iterate trough region list */
  pOctant root;           /* root of the partition tree */
  int     i;              /* index counter */
  int     count,          /* count for leaf nodes in partition tree */
          proc,           /* proc leaf node of parition tree belongs to */
          extra,          /* extra leaf node flag, if not evenly divisible */
          remainder;      /* remainder of node, or processors to fill */
  pOctant cursor,         /* cursor to iterate through octant list */
          cursor2,        /* another cursor to iterate through octant list */
          parent;         /* parent of an octant */
  int level,              /* number of levels of refinement */
      n,                  /* index counter */
      part;               /* partition counter */
  Map *array;             /* map of which processors own which octants */
  int hold;               /* used for calculating partition divisions */
  int ierr = 0;

#ifdef KDDKDD_NEW_BOUNDS_GEOM_QUERY_FN
  double bounds[6] = {DBL_MAX,DBL_MAX,DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX};
  COORD global_min, global_max;
#endif /* KDDKDD_NEW_BOUNDS_GEOM_QUERY_FN */
  int nroots = 0;
  /*test*/
  /* COORD gmin,gmax; */

  OCT_Global_Info *OCT_info = (OCT_Global_Info *) (zz->LB.Data_Structure);

  ZOLTAN_TRACE_ENTER(zz, yo);
  /*
   * If there are no objects on this processor, do not create a root octant.
   * The partitioner will probably assign objects to this processor
   */
  if(zz->Get_Num_Obj == NULL) {
    fprintf(stderr, "OCT %s\n\t%s\n", "Error in octree load balance:",
	    "Must register ZOLTAN_NUM_OBJ_FN function");
    abort();
  }
  *c3 = num_objs = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  if (ierr) {
    fprintf(stderr, "OCT [%d] %s: Error returned from user defined "
                    "Get_Num_Obj function.\n", zz->Proc, yo);
    exit (-1);
  }
  ptr1 = NULL;

  ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_get_bounds");
  /* Need A Function To Get The Bounds Of The Local Objects */
  Zoltan_Oct_get_bounds(zz, &ptr1, &num_objs, min, max, oct_wgtflag, c0);
  
#ifndef KDDKDD_NEW_BOUNDS_GEOM_QUERY_FN
  /* For now, don't want to add the new query function to Zoltan. */
  /* Zoltan_Oct_get_bounds appears to compute the global min and max from */
  /* the object input. */
  vector_set(OCT_info->OCT_gmin, min);
  vector_set(OCT_info->OCT_gmax, max);
#else
  /*test*/
  /*getMaxBounds(&gmin, &gmax);*/
  if(zz->Get_Bounds_Geom == NULL) {
    fprintf(stderr, "OCT %s\n\t%s\n", "Error in octree load balance:",
	    "Must register Get_Bounds_Geom function");
    abort();
  }
  zz->Get_Bounds_Geom(zz->Get_Bounds_Geom_Data, bounds, &ierr); 
  
  MPI_Allreduce(&(bounds[0]), &(global_min[0]), 3, 
		MPI_DOUBLE, MPI_MIN, zz->Communicator);
  MPI_Allreduce(&(bounds[3]), &(global_max[0]), 3,
		MPI_DOUBLE, MPI_MAX, zz->Communicator);
  vector_set(OCT_info->OCT_gmin, global_min);
  vector_set(OCT_info->OCT_gmax, global_max);
#endif
  /* 
   * the following code segment was added to create a pseudo global octree
   * needed for the partitioner. The basic idea is to regroup all the
   * regions into something close to an octree partitioning and build the
   * tree from that.
   * NOTE: This way of doing things is very costly, especially when calling
   * this for the first time on a mesh not partitioned in an octree style
   * partitioning.
   */

    level = 0;                                    /* initialize level count */

  /* 
   * if more than 1 processor, need to find what level of refinement needed
   * to initially partition bounding box among the processors 
   */

  
    if(zz->Num_Proc > 1) {
      n = zz->Num_Proc;
      if(OCT_info->OCT_dimension == 2)
	hold = 4;
      else
	hold = 8;
      remainder = hold;
      for(; remainder > 0; level++) {
	int pr = (int)POW(hold, level);
	remainder = n - pr;
      }
      level--;
    }
  ZOLTAN_TRACE_DETAIL(zz, yo, "Before createpartree");

  if(createpartree) {
    /* create the global root octant */
    root = Zoltan_Oct_POct_new(OCT_info);
    Zoltan_Oct_setbounds(root, OCT_info->OCT_gmin, OCT_info->OCT_gmax);
    /* Zoltan_Oct_setOrientation(root, 0); */
  
    /* subdivide to as many levels as calculated */
    for(i=0; i<level; i++) {
      cursor = root;
      while(cursor != NULL) {
	if(Zoltan_Oct_isTerminal(cursor)) {
	  cursor2 = Zoltan_Oct_POct_nextDfs(OCT_info, cursor);
	  Zoltan_Oct_terminal_refine(zz, cursor, 0);
	  cursor = cursor2;
	}
	else 
	  cursor = Zoltan_Oct_POct_nextDfs(OCT_info, cursor);
      }
    }
    
#if 0
    if(zz->Proc == 0)
      for(i=0; i<8; i++)
	if(Zoltan_Oct_child(root, i) == NULL)
	  fprintf(stderr,"NULL child pointer\n");
	else
	  fprintf(stderr, "child %d exists\n", i);
#endif

  ZOLTAN_TRACE_DETAIL(zz, yo, "Before create map array");
    /* this part creates the map array */
    if(OCT_info->OCT_dimension == 2) {
      hold = (int)POW(4, level);                 /* ignoring the z+ octants */
      if(hold == 0)
	hold = 1;
    }
    else
      hold = (int)POW(8, level);

    part = hold / zz->Num_Proc;          /* how many octants per partition */
    remainder = hold % zz->Num_Proc; /* extra octants, not evenly divisible */
    extra = zz->Num_Proc - remainder;/* where to start adding extra octants */
    array = (Map *) ZOLTAN_MALLOC(hold * sizeof(Map));   /* alloc map array */
    if(array == NULL) {
      fprintf(stderr, "OCT ERROR on proc %d, could not allocate array map\n",
	      zz->Proc);
      abort();
    }
    /* initialize variables */
    proc = 0;
    count = 0;
    i = 0;
    cursor = root; 
    while(cursor != NULL) {
      cursor2 = Zoltan_Oct_POct_nextDfs(OCT_info, cursor);
      if((Zoltan_Oct_isTerminal(cursor)) && (i < hold)) {
	if(proc == extra) {
	  part++;
	  extra = -1;
	}
	if(count != part) {
	  array[i].npid = proc;
	  array[i].list = RL_initRootList();
	  Zoltan_Oct_bounds(cursor, min, max);
	  vector_set(array[i].min, min);
	  vector_set(array[i].max, max);
	  count++;
	}
	else {
	  count = 1;
	  proc++;
	  array[i].npid = proc;
	  array[i].list = RL_initRootList();
	  Zoltan_Oct_bounds(cursor, min, max);
	  vector_set(array[i].min, min);
	  vector_set(array[i].max, max);
	}
	if(proc == zz->Proc) {
	  array[i].npid = -1;
          /* KDDKDD Added RL_freeList below.  The 
           * KDDKDD implementation from RPI leaked memory because the 
           * KDDKDD test cases for setting array[i].list were not mutually 
           * KDDKDD exclusive.  Freeing the list produces the result we got
           * KDDKDD before, without the memory leak.
           */
	  /* LGG --  it seems to me that this array[i].list assignment is
	   * not really necessary. It looks as though it has already been
	   * assigned with the same information from the prev if-else
	   * commented out RL_freeList(), and RL_initRootList()
	   */
          /*RL_freeList(&(array[i].list));*/
          /* KDDKDD End addition */
	  /*array[i].list = RL_initRootList();*/
	  parent = Zoltan_Oct_parent(cursor);
	  if(parent != NULL)
	    Zoltan_Oct_setchild(parent, cursor->which, NULL);
	  /* octant into local root list */
 	  Zoltan_Oct_POct_setparent(OCT_info, cursor, NULL, -1);
	  Zoltan_Oct_setMapIdx(cursor, i);
	  nroots++;
	  /* Zoltan_Oct_POct_setparent(OCT_info, cursor, NULL, zz->Proc);     
             octant into local root list */
	}
	i++;
      }
      cursor = cursor2;
    } 
    RootList = Zoltan_Oct_POct_localroots(OCT_info); 
    RootOct = RL_nextRootOctant(&RootList);
    if(RootOct != root) {
      /* KDDKDDFREE changed root to &root to allow root to be reset to NULL */
      Zoltan_Oct_POct_delTree(OCT_info,&root);
    }
    
    OCT_info->map = array;
    OCT_info->mapsize = hold;
  }

  /* 
   * attach the regions to the root... Zoltan_Oct_fix will create the octree
   * starting with the root and subdividing as needed 
   */    
  num_extra = Zoltan_Oct_fix(zz, ptr1, num_objs);
 
  ZOLTAN_TRACE_DETAIL(zz, yo, "Calling Zoltan_Oct_migreg_migrate_orphans");
  Zoltan_Oct_migreg_migrate_orphans(zz, ptr1, num_extra, level, OCT_info->map,
				    c1, c2);

  /* ZOLTAN_FREE(&array); */
  while(ptr1 != NULL) {
    ptr = ptr1->next;
    ZOLTAN_FREE(&(ptr1->Global_ID));
    ZOLTAN_FREE(&(ptr1->Local_ID));
    ZOLTAN_FREE(&ptr1);
    ptr1 = ptr;
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
}

/****************************************************************************/

static void Zoltan_Oct_get_bounds(ZZ *zz, pRegion *ptr1, int *num_objs, 
		   COORD min, COORD max, int wgtflag, float *c0) 
{
  char *yo = "Zoltan_Oct_get_bounds";
  ZOLTAN_ID_PTR obj_global_ids = NULL; 
  ZOLTAN_ID_PTR obj_local_ids = NULL;
  int *parts = NULL;   /* Input partition assignments; currently unused. */
  float *obj_wgts = NULL;
  double *geom_vec = NULL;
  float objwgt;        /* Temporary value of an object weight; used to pass
                          0. to initialize_regions when wgtflag == 0. */
  ZOLTAN_ID_PTR lid;   /* Temporary pointer to a local ID; used to pass NULL 
                          to initialize_regions when NUM_LID_ENTRIES == 0. */
  int num_dim;
  int i;
  pRegion tmp, ptr;
  COORD global_min, global_max;
  double PADDING = 0.0000001;
  int ierr = 0;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;

  /* Initialization */
  max[0] = max[1] = max[2] = -DBL_MAX;
  min[0] = min[1] = min[2] =  DBL_MAX;

  ierr = Zoltan_Get_Obj_List(zz, num_objs, &obj_global_ids, &obj_local_ids,
                             wgtflag, &obj_wgts, &parts);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                   "Error returned from user function Zoltan_Get_Obj_List.");
    exit (-1);
  }

  ierr = Zoltan_Get_Coordinates(zz, *num_objs, obj_global_ids, obj_local_ids,
                                &num_dim, &geom_vec);

  for (i = 0; i < (*num_objs); i++) {
    lid = (num_lid_entries ? obj_local_ids + i*num_lid_entries : NULL);
    objwgt = (wgtflag ? obj_wgts[i] : 0.);
    initialize_region(zz, &(ptr), obj_global_ids + i*num_gid_entries,
                      lid, wgtflag, objwgt,
                      num_dim, geom_vec + i*num_dim);
    if (i == 0) {
      tmp = ptr;
      *c0 = (float)tmp->Weight;
      vector_set(min, tmp->Coord);
      vector_set(max, tmp->Coord);
      *ptr1 = tmp;
    }
    else {
      *c0 += (float)ptr->Weight;
      /* the following is really a hack, since it has no real basis 
         in vector mathematics.... */
      if(ptr->Coord[0] < min[0])
        min[0] = ptr->Coord[0];
      if(ptr->Coord[1] < min[1])
        min[1] = ptr->Coord[1];
      if(ptr->Coord[2] < min[2])
        min[2] = ptr->Coord[2];
      if(ptr->Coord[0] > max[0])
        max[0] = ptr->Coord[0];
      if(ptr->Coord[1] > max[1])
        max[1] = ptr->Coord[1];
      if(ptr->Coord[2] > max[2])
        max[2] = ptr->Coord[2];
      tmp->next = ptr;
      tmp = tmp->next;
    }
  
    ptr = NULL;
  }
  ZOLTAN_FREE(&obj_global_ids);
  ZOLTAN_FREE(&obj_local_ids);
  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&obj_wgts);
  ZOLTAN_FREE(&geom_vec);
  
  MPI_Allreduce(&(min[0]), &(global_min[0]), 3, 
		MPI_DOUBLE, MPI_MIN, zz->Communicator);
  MPI_Allreduce(&(max[0]), &(global_max[0]), 3,
		MPI_DOUBLE, MPI_MAX, zz->Communicator);

  max[0] = global_max[0];
  max[1] = global_max[1];
  max[2] = global_max[2];
  min[0] = global_min[0];
  min[1] = global_min[1];
  min[2] = global_min[2];
  
  /* hack used for sample program since working in 2D -- */
  /* causes problems for refining the octree */
  if(max[2] == min[2])
    max[2] = 1.0;

  for(i=0; i<3; i++) {
    /* min[i] -= PADDING; */
    max[i] += PADDING;
  }

  return;
}

/****************************************************************************/
/*
 *  Function that initializes the region data structure.  It uses the 
 *  global ID, coordinates and weight provided by the application.  
 */
static void initialize_region(ZZ *zz, pRegion *ret, ZOLTAN_ID_PTR global_id,
                              ZOLTAN_ID_PTR local_id, int wgtflag, float wgt,
                              int num_dim, double *geom_vec) 
{
  pRegion reg;
  int i;
  reg = (pRegion) ZOLTAN_MALLOC(sizeof(Region));
  *ret = reg;
  reg->Global_ID = ZOLTAN_MALLOC_GID(zz);
  reg->Local_ID = ZOLTAN_MALLOC_LID(zz);
  ZOLTAN_SET_GID(zz, reg->Global_ID, global_id);
  ZOLTAN_SET_LID(zz, reg->Local_ID, local_id);
  reg->Proc = zz->Proc;
  /* reg->Proc = 0; */
  reg->Coord[0] = reg->Coord[1] = reg->Coord[2] = 0.0;
  for (i = 0; i < num_dim; i++)
    reg->Coord[i] = geom_vec[i];

#if 0
  Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
    fprintf(stderr, "Result info on %d: %d %d %d  %lf  %lf  %lf\n", zz->Proc, 
	    reg->Local_ID, reg->Global_ID, reg->Proc,
	    reg->Coord[0], reg->Coord[1], reg->Coord[2]); 
  Zoltan_Print_Sync_End(zz->Communicator, TRUE);
#endif

  if (wgtflag)
    reg->Weight = wgt;
  else
    reg->Weight = 1;

  reg->next = NULL;
}

/****************************************************************************/
/*
 * int Zoltan_Oct_fix(pMesh mesh)
 *
 * Clear all region pointers from the octree,
 * reinsert all mesh regions, 
 * and coarsen or refine the octree as needed.
 * 
 * Return the number of regions that could
 * not be inserted.
 *
 */
static int Zoltan_Oct_fix(ZZ *zz, pRegion Region_list, int num_objs) 
{
  int nreg;                               /* number of regions not inserted */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  /* initalize variables */
  nreg=0;
  oct_nref=0;
  oct_ncoarse=0;
  /* associate the objects to the octants */
  nreg=Zoltan_Oct_global_insert_object(zz, Region_list, num_objs);
  Zoltan_Oct_global_dref(zz, OCT_info); 

  return(nreg);
}

/****************************************************************************/
/*
 * int Zoltan_Oct_global_insert_object()
 *
 * Try to insert all regions in the mesh into the existing
 * local octree.  Return the number of insertion failures.
 *
 */
static int Zoltan_Oct_global_insert_object(ZZ *zz, pRegion Region_list, 
					   int num_objs) 
{
  pRegion region;                             /* region to be attached */
  int count;                                  /* count of failed insertions */
  int i;

  /* initialize variables */
  count=0;
  region = Region_list;

  /* get the next region to be inserted */
  for(i=0; i<num_objs; i++) {
#if 0
    if(region != NULL)
      printf("\n%lf   %lf   %lf\n", 
	     region->Coord[0], region->Coord[1], region->Coord[2]); 
#endif

    if (!Zoltan_Oct_global_insert(zz, region)) {
      /* obj has no octant association increment "orphan" counter */
      count++;
      region->attached = 0;
    }
    else
      region->attached = 1;
    region = region->next;
  }

  return(count);
}

/****************************************************************************/
/*
 * Zoltan_Oct_global_insert(region)
 *
 * try to insert the region into any of the local-rooted subtrees
 *
 * return the octant pointer if successful, or NULL if region's
 * centroid does not lie within the local octree.  This could
 * be due to centroid not lying within any local root, or also
 * if some local root's subtree is off-processor.
 *
 * During insertion, refinement may be performed.
 * In that case, the returned octant is an ancestor
 * of the octant to which the region is attached.
 *
 */
static pOctant Zoltan_Oct_global_insert(ZZ *zz, pRegion region) 
{
  pOctant oct;                            /* octree octant */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);

  oct = NULL;
  /* find the octant which the object lies in */
  oct=Zoltan_Oct_global_find(OCT_info,region->Coord);
  /* if octant is found, try to insert the region */
  if (oct)
    if (!Zoltan_Oct_subtree_insert(zz, oct, region)) {  /* inserting region */
      fprintf(stderr,"OCT Zoltan_Oct_global_insert: insertion failed\n");
      abort();
    }

  return(oct);
}

/****************************************************************************/
/*
 * Zoltan_Oct_subtree_insert(oct,region)
 *
 * Insert region in oct, carrying out multiple refinement if necessary
 */
int Zoltan_Oct_subtree_insert(ZZ *zz, pOctant oct, pRegion region) 
{
OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure); 
  /* if oct is not terminal, find leaf node the centroid can be attahced to */
  if (!Zoltan_Oct_isTerminal(oct)) {
    oct=Zoltan_Oct_findOctant(OCT_info, oct,region->Coord);
  }

  if (!oct)
    return(0);

  /* add the region to the octant */
  Zoltan_Oct_addRegion(zz, oct, region);
 
  /* check if octant has too many regions and needs to be refined */
  /* KDDKDD  Replaced the following to allow multiple regions with the 
   * KDDKDD  same coordinates to be placed in a single octant.
  if(Zoltan_Oct_nRegions(oct) > MAXOCTREGIONS)
   */
  if(Zoltan_Oct_nUniqueRegions(OCT_info, oct) > MAXOCTREGIONS)
    Zoltan_Oct_terminal_refine(zz, oct,0);   /* After, dest may be nonterm */

  return(1);
}

/****************************************************************************/
/*
 * pOctant Zoltan_Oct_global_find(COORD centroid)
 *
 * Return the octant in the local octree that contains the
 * given point, if it exists.  Otherwise, return NULL.
 *
 */
static pOctant Zoltan_Oct_global_find(OCT_Global_Info *OCT_info,COORD point) 
{  
  pRList  RootList;               /* list of all local roots */
  pOctant RootOct;                /* root octree octant */
  pOctant oct;                    /* octree octant */
  COORD min,                      /* minimum bounds coordinates */
        max;                      /* maximum bounds coordinates */

  /* get a list of all local roots */
  RootList = Zoltan_Oct_POct_localroots(OCT_info); 
  oct = NULL;

  /* iterate through root list to find if point lies inside root's bounds */
  while ((!oct) && (RootOct = RL_nextRootOctant(&RootList)) ) {
    Zoltan_Oct_bounds(RootOct,min,max);
    /* ATTN: Zoltan_Oct_in_box may need adjusting for the lower bounds */
    /* check if point fits inside */
    if (Zoltan_Oct_in_box(OCT_info,point,min,max)) {  
      /* find the exact octant for point */
      oct=Zoltan_Oct_findOctant(OCT_info,RootOct,point);
    }
  }

  return(oct);
}

/****************************************************************************/
/*
 * Zoltan_Oct_findOctant(oct, coord)
 *   (replaces : PO_findOctant(oct,coord))
 *  
 *
 * find the octant in a subtree containing coord (if coord
 * is not in the subtree, returns the closest octant in the subtree).
 * NOTE: return NULL if we hit an off-processor link
 *
 */
static pOctant Zoltan_Oct_findOctant(OCT_Global_Info *OCT_info,pOctant oct, 
				     COORD coord) 
{
  pOctant child;                                  /* child of an octant */
  int cnum;                                       /* child number */
  
  /* if octant is terminal, then this is the right octant */
  if (Zoltan_Oct_isTerminal(oct))
    return(oct);

  /* find closest child to coord */
  cnum = Zoltan_Oct_child_which_wrapper(OCT_info,oct,coord);
  child = Zoltan_Oct_child(oct, cnum);                    /* get that child */
  /* ATTN: are these local checks necessary? */
  /* if ( !Zoltan_Oct_POct_local(child) ) */
  if(!Zoltan_Oct_POct_local(OCT_info,oct, cnum))  /* make sure oct is local */
    return(NULL);

  /* recursivly search down the tree */
  return(Zoltan_Oct_findOctant(OCT_info,child,coord));
}

/****************************************************************************/
/*
 * Zoltan_Oct_terminal_refine(oct)
 *
 * subdivide a terminal octant and divvy up
 * its regions to the 8 children; recurse if
 * necessary to satisfy MAXOCTREGIONS
 *
 */
static void Zoltan_Oct_terminal_refine(ZZ *zz, pOctant oct,int count) 
{
  COORD min,                     /* coordinates of minimum bounds of region */
        max,                     /* coordinates of maximum bounds of region */
        origin;                  /* origin of region */
  pOctant child[8];              /* array of child octants */
  int cnum;                      /* child number */
  int i;                         /* index counter */
  pRegion region;                /* a region to be associated to an octant */
  pRegion entry;
  COORD cmin[8], cmax[8];
  OCT_Global_Info *OCT_info = (OCT_Global_Info *) (zz->LB.Data_Structure);

  for(i=0;i<3;i++)
    min[i] = max[i] = 0;
  
  /* upper limit of refinement levels */
  /* ATTN: may not be used anymore, but can be put back in if necessary */
  if (count>=20) {
    fprintf(stderr, "OCT ERROR: Zoltan_Oct_terminal_refine: bailing out at "
                    "10 levels\n");
    abort();
  }
  oct_nref++;                               /* increment refinement counter */

  /* octant should be terminal in order to be refined (subdivided) */
  if (!Zoltan_Oct_isTerminal(oct)) {
    fprintf(stderr,"OCT ref_octant: oct not terminal\n");
    abort();
  }

  /* get the bounds of an octant */
  Zoltan_Oct_bounds(oct,min,max);
  /* calculate the origin from the bounds */
  Zoltan_Oct_bounds_to_origin(min,max,origin);

  region = Zoltan_Oct_regionlist(oct);     /* Get list while still terminal */
  oct->list = NULL;  /* remove regions from octant, it won't be terminal */

  /* create the children and set their id's */

  Zoltan_Oct_child_bounds_wrapper(OCT_info,oct, cmin, cmax);
  for (i=0; i<8; i++) {
    if(OCT_info->OCT_dimension == 2) {
      /* KDDKDD 3/01 see changes to Zoltan_Oct_child_bounds_wrapper that allow this
       * KDDKDD 3/01 test to work for GRAY and HILBERT mappings.
       */
      if(cmin[i][2] > OCT_info->OCT_gmin[2]) {     /* ignore the z+ octants */
	child[i] = NULL;
	continue;
      }
    }
    
    child[i]=Zoltan_Oct_POct_new(OCT_info);          /* create a new octant */
    child[i]->dir = Zoltan_Oct_get_child_dir(OCT_info, oct->dir, i);
                  /* create a new octant */
    /* set the child->parent link */
    Zoltan_Oct_POct_setparent(OCT_info, child[i], oct, zz->Proc);
    Zoltan_Oct_setchildnum(child[i], i);       /* which child of the parent */
    Zoltan_Oct_setchild(oct, i, child[i]);    /* set the parent->child link */
#ifdef LGG_MIGOCT
    Zoltan_Oct_setID(child[i], Zoltan_Oct_nextId());    /* set child id num */
#endif /* LGG_MIGOCT */
    Zoltan_Oct_setbounds(child[i], cmin[i], cmax[i]);   /* set child bounds */
    Zoltan_Oct_setCpid(oct, i, zz->Proc);    /* set child to be a local oct */
    /* Zoltan_Oct_setOrientation(child[i], 
	      Zoltan_Oct_child_orientation(oct->orientation, oct->which));  */
  }

  /* assign newly created children to child array*/
  if(OCT_info->OCT_dimension == 3) {
    if(Zoltan_Oct_children(oct, child) != 8) {
      /* 
       * if subdivision of oct was successful, oct should have 8 children; 
       * thus a return value of 0 here is a fatal error
       */
      fprintf(stderr, "OCT ref_octant: subdivide failed, %d children.\n",
	      Zoltan_Oct_children(oct, child));
      abort();
    }
  }
  else
    if(Zoltan_Oct_children(oct, child) != 4) {
      /* 
       * if subdivision of oct was successful, oct should have 4 children; 
       * thus a return value of 0 here is a fatal error
       */
      fprintf(stderr, 
	      "OCT ref_octant:subdivide failed, %d children, expected 4\n",
	      Zoltan_Oct_children(oct, child));
      abort();
    }

  /* iterate through and find which child each region should belong to */
  while(region != NULL) {
    entry = region->next;
    cnum=Zoltan_Oct_child_which_wrapper(OCT_info,oct, region->Coord);
    /* add region to octant's regionlist */
    Zoltan_Oct_addRegion(zz, child[cnum], region);
    ZOLTAN_FREE(&(region->Global_ID));
    ZOLTAN_FREE(&(region->Local_ID));
    ZOLTAN_FREE(&region);
    region = entry;
  }

  for (i=0; i<8; i++)                                          /* Recursion */
    if(child[i] != NULL)
      /* KDDKDD  Replaced the following to allow multiple regions with the
       * KDDKDD same coordinates to be placed in the same octant.
      if (Zoltan_Oct_nRegions(child[i]) > MAXOCTREGIONS) {
       */
      if (Zoltan_Oct_nUniqueRegions(OCT_info,child[i]) > MAXOCTREGIONS) {
	Zoltan_Oct_terminal_refine(zz, child[i],count+1);
      }
}

/****************************************************************************/
/*
 * Zoltan_Oct_global_dref()
 * 
 * refine and derefine octree as necessary by number of
 * regions in each octant
 *
 */
static void Zoltan_Oct_global_dref(ZZ *zz, OCT_Global_Info *OCT_info) 
{
  pRList  RootList;                           /* list of all local roots */
  pOctant RootOct;

  RootList = Zoltan_Oct_POct_localroots(OCT_info);
  while ((RootOct = RL_nextRootOctant(&RootList))) 
    Zoltan_Oct_subtree_dref(zz, OCT_info, RootOct);
}

/****************************************************************************/
/*
 * Zoltan_Oct_subtree_dref(oct)
 *
 * Coarsen octree so that leaf octants do not have less than MINOCTREGIONS 
 * regions. Refinement takes precedence, so coarsening will not take place 
 * unless all subtrees agree on it.
 */
static int Zoltan_Oct_subtree_dref(ZZ *zz, OCT_Global_Info *OCT_info,
				   pOctant oct) 
{
  pOctant child;                        /* child of an octant */
  int coarsen;                          /* flag to indicate need to coarsen */
  int i;                                /* index counter */
  int nregions;                         /* number of regions */
  int total;                            /* total number of regions */

  /* if terminal octant, cannot coarsen on own */
  if (Zoltan_Oct_isTerminal(oct)) {
    nregions=Zoltan_Oct_nRegions(oct);

    if (nregions > (MAXOCTREGIONS*2) ) {
      fprintf(stderr, 
	      "OCT Zoltan_Oct_subtree_dref: warning: too many (%d) regions "
	     "in oct (id=%d)\n",Zoltan_Oct_nRegions(oct),Zoltan_Oct_id(oct));
      return(-1);
    }
    else
      if (nregions < MINOCTREGIONS)
	return(nregions);
      else
	return(-1);
  }

  coarsen=1;                                        /* assume to be coarsen */
  total=0;

  /* look at each child, see if they need to be coarsened */
  for (i=0; coarsen && i<8; i++) {
    child = Zoltan_Oct_child(oct,i);
    /* if child is off processor cannot coarsen */
    /* if (!Zoltan_Oct_local(child) ) */
    if(!Zoltan_Oct_POct_local(OCT_info, oct, i))
      coarsen=0;
    else {
      /* get the number of region of the child */
      nregions=Zoltan_Oct_subtree_dref(zz, OCT_info,child);
      if (nregions<0)
	coarsen=0;
      else
	total+=nregions;                          /* total the region count */
    }
  }

  /* check if octant can be coarsened */
  if (coarsen && total<MAXOCTREGIONS) {
    Zoltan_Oct_terminal_coarsen(zz, OCT_info,oct);
    
    if (total < MINOCTREGIONS)
      return(total);
    else
      return(-1);
  }

  return(-1);
}

/****************************************************************************/
/*
 * Zoltan_Oct_terminal_coarsen(oct)
 *
 * remove octant's children, accumulating regions
 * to octant
 *
 */
static void Zoltan_Oct_terminal_coarsen(ZZ *zz, OCT_Global_Info *OCT_info,
					pOctant oct) 
{
  pOctant child;                        /* child of an octant */
  pRegion region;                       /* region associated with an octant */
  int i;                                /* index counter */
  pRegion regionlist[8];                /* an array of region lists */

  oct_ncoarse++;                        /* increment coarsening counter */

  for(i=0; i<8; i++) {
    /* get the ith child of an octant */
    child = Zoltan_Oct_child(oct,i);
    
    /* cannot coarsen if child is off-processor */
    /* if(!Zoltan_Oct_POct_local(child)) X */
    /* cannot be off-processor */
    if(!Zoltan_Oct_POct_local(OCT_info, oct, i)) {
      fprintf(stderr,"OCT Zoltan_Oct_terminal_coarsen: child not local\n");
      abort();
    }
    
    if(!Zoltan_Oct_isTerminal(child)) {
      fprintf(stderr,"OCT Zoltan_Oct_terminal_coarsen: child not terminal\n");
      abort();
    }
    
    /* get each child's region list */
    regionlist[i] = Zoltan_Oct_regionlist(child);
    
    /* delete each child */
    /* KDDKDDFREE Change child to &child. */
    Zoltan_Oct_POct_free(OCT_info, &child);
    oct->child[i] = NULL;
  }
  oct->numChild = 0;
  /* 
   * copy contents of each region list into region list
   * of coarsened parent (which is now a terminal octant) 
   */
  for(i=0; i < 8; i++) {
    region = regionlist[i];
    /* go through the regionlist and add to octant */
    while(region != NULL) {
      Zoltan_Oct_addRegion(zz, oct,region);           
      region = region->next;
    }
  }
}

/****************************************************************************/
static void Zoltan_Oct_set_maxregions(int max)
{ 
  if (max < 1) {
    fprintf(stderr, "OCT Warning Zoltan_Oct_set_maxregions(): %s\n",
	    "illegal input, using default.");
    MAXOCTREGIONS = MAXOCTREGIONS_DEF;
  }
  else
    MAXOCTREGIONS = max; 
}

/****************************************************************************/
static void Zoltan_Oct_set_minregions(int min)
{ 
  if (min < 1) {
    fprintf(stderr, "OCT Warning Zoltan_Oct_set_minregions(): %s\n",
	    "illegal input, using default.");
    MINOCTREGIONS = MINOCTREGIONS_DEF;
  }
  else
    MINOCTREGIONS = min; 
}

/****************************************************************************/
#ifdef LGG_MIGOCT
void Zoltan_Oct_resetIdCount(int start_count)
{
  IDcount = start_count;
}

/****************************************************************************/

int Zoltan_Oct_nextId(void)
{
  return ++IDcount;
}

/****************************************************************************/
typedef struct
{
  pOctant ptr;
  int id;
} Rootid;

/****************************************************************************/

static int idcompare(Rootid *i, Rootid *j)
{
  return( (i->id) - (j->id) );
}
/****************************************************************************/
#endif /* LGG_MIGOCT */

/****************************************************************************/


/*
 * oct_global_clear()
 *
 * delete all regions from all octants on the local processor
 *
 */
static void Zoltan_Oct_global_clear(OCT_Global_Info * OCT_info)
{ 
  pRList  RootList;                 /* list of all local roots */
  pOctant Oct;

  /* 
   * iterate through the list of local roots 
   * traverse down the subtree 
   * delete regions associated with octant 
   */
  RootList = Zoltan_Oct_POct_localroots(OCT_info);
  while ((Oct = RL_nextRootOctant(&RootList))) {
    while(Oct) {
      if(Zoltan_Oct_isTerminal(Oct))
	Zoltan_Oct_clearRegions(Oct);
      Oct = Zoltan_Oct_POct_nextDfs(OCT_info, Oct);
    }
  }
}

#if 0
/****************************************************************************/

/* global variable from POC library */
/*
 * pOctant oct_localTree_nextDfs(pOctant oct, void **temp)
 * (previously: pOctant oct_global_nextDfs(pParoct oct, void **temp))
 *
 * do a dfs of the local octree on each processor.
 * subtrees are traversed in the order of the local roots
 * list, and each subtree is done in dfs order.
 *
 * Undefined action will take place if you modify the
 * tree while traversing.
 *
 */
pOctant oct_localTree_nextDfs(pOctant oct, void **temp)
{
  if (*temp==NULL)           /* Start of traversal */
    return(PList_next(OCT_info->OCT_rootlist,temp));


  if (oct=Zoltan_Oct_POct_nextDfs(oct, 0))
    return(oct);

  return(PList_next(OCT_info->OCT_rootlist,temp));
}
#endif  /* 0 */

/****************************************************************************/
/*
 * int Zoltan_Oct_nUniqueRegions(pOctant octant)
 * return the number of regions in the octant's list
 * KDDKDD  Return the number of regions with Unique coordinates; regions that
 * KDDKDD  have the same coordinates as some other region in the octant
 * KDDKDD  are not counted.
 * KDDKDD  This change was added so that two regions with the same coordinates
 * KDDKDD  are allowed to be in the same octant, even if having them their
 * KDDKDD  put the number of regions in the octant over MAXOCTREGIONS.
 * KDDKDD  For the default value of MAXOCTREGIONS=1, this workaround is 
 * KDDKDD  necessary when any objects have the same coordinates.
 * KDDKDD  Perhaps a better solution can be found later.
 */
int Zoltan_Oct_nUniqueRegions(OCT_Global_Info *OCT_info, pOctant oct) {
  pRegion ptr;                    /* pointer to iterate through region list */
  int count;                      /* count of number of regions in list */
  pRegion ptr2;                   /* pointer to iterate through prev regions*/
  int same_coords;

  if (oct == NULL) 
    return 0;

  if(!Zoltan_Oct_isTerminal(oct)) 
    return 0;

  count = 0;
  ptr = oct->list;
  while(ptr != NULL) {
    ptr2 = oct->list;
    same_coords = FALSE;
    while (ptr2 != ptr) {
      if (Zoltan_Oct_CompareCoords(OCT_info->OCT_dimension, 
                               ptr->Coord, ptr2->Coord)) {
        same_coords = TRUE;
        break;
      }
      else
        ptr2 = ptr2->next;
    }
    if (!same_coords) 
      count ++;
    ptr = ptr->next;
  }
  return(count);
}
/****************************************************************************/
/*
 * Routine to compare coordinates
 */
static int Zoltan_Oct_CompareCoords(int dim, COORD pt1, COORD pt2) 
{
int i;
int ret = TRUE;

  for (i = 0; i < dim; i++)
    if (pt1[i] != pt2[i]) {
      ret = FALSE;
      break;
    }
 
  return ret;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
