/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
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
#include <values.h>
#define POW(a,b) pow((double)(a),(double)(b))

/*test*/
/*extern void getMaxBounds(void *, double *, double *); */

/***************************  PROTOTYPES *************************************/
static int LB_Oct_nUniqueRegions(OCT_Global_Info *OCT_info, pOctant oct);
static int LB_Oct_CompareCoords(int dim, COORD pt1, COORD pt2);
static void initialize_region(LB *, pRegion *, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float);
static int lb_oct_init(LB *lb, int *num_import, ZOLTAN_ID_PTR *import_global_ids,
  ZOLTAN_ID_PTR *import_local_ids, int **import_procs, int oct_dim, int oct_method,
  int oct_maxoctregions, int oct_minoctregions, int oct_output_level, int oct_wgtflag); 
static void    LB_oct_gen_tree_from_input_data(LB *lb, int, int *c1, int *c2,
                                               int *c3, float *c0, int createpartree);
static pOctant LB_oct_global_insert(LB *, pRegion region);
static void    LB_oct_terminal_refine(LB *, pOctant oct,int count);



/*****************************************************************************/
/* NOTE: be careful later about region lists for nonterminal octants */

static int oct_nref=0;                             /* number of refinements */
static int oct_ncoarse=0;                          /* number of coarsenings */
#ifdef LGG_MIGOCT
static int IDcount = 0;                            /* renumbering of octants */
#endif
static int MAXOCTREGIONS = 1;
static int MINOCTREGIONS = 1;                    /* minimum number of regions per octant */

/*****************************************************************************/
/* parameters for the octpart method.  Used in  */
/* LB_Set_Octpart_Param and LB_octpart          */
static PARAM_VARS OCT_params[] = {
  { "OCT_DIM",          NULL, "INT" },
  { "OCT_METHOD",       NULL, "INT" },
  { "OCT_MAXOBJECTS",   NULL, "INT" },
  { "OCT_MINOBJECTS",   NULL, "INT" },
  { "OCT_OUTPUT_LEVEL", NULL, "INT" },
  {  NULL,              NULL,  NULL }};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Set_Octpart_Param(
  char *name,                 /* name of variable */
  char *val                   /* value of variable */
)
{
int status;
PARAM_UTYPE result;           /* value returned from LB_Check_Param */
int index;                    /* index returned from LB_Check_Param */

  status = LB_Check_Param(name, val, OCT_params, &result, &index);
  return(status);
}

/*****************************************************************************/
int LB_octpart(
  LB *lb,                       /* The load-balancing structure with info for
                                   the OCTPART balancer.                     */
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
  int *num_export,              /* Not computed; return -1. */
  ZOLTAN_ID_PTR *export_global_ids, /* Not computed. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Not computed. */
  int **export_procs            /* Not computed. */
) 
{
int oct_dim = 3;              /* Dimension of method to be used (2D or 3D)   */
int oct_method = 2;           /* Flag specifying curve to be used.           */
int oct_maxoctregions = 1;    /* max # of objects in leaves of octree.       */
int oct_minoctregions = 1;    /* min # of objects in leaves of octree.       */
int oct_output_level = 0;     /* Flag specifying amount of output.           */
int oct_wgtflag = 0;          /* Flag specifying use of object weights.      */
int error = FALSE;            /* error flag                                  */


  LB_Bind_Param(OCT_params, "OCT_DIM", (void *) &oct_dim);
  LB_Bind_Param(OCT_params, "OCT_METHOD", (void *) &oct_method);
  LB_Bind_Param(OCT_params, "OCT_MAXOBJECTS", (void *) &oct_maxoctregions);
  LB_Bind_Param(OCT_params, "OCT_MINOBJECTS", (void *) &oct_minoctregions);
  LB_Bind_Param(OCT_params, "OCT_OUTPUT_LEVEL", (void *) &oct_output_level);

  LB_Assign_Param_Vals(lb->Params, OCT_params, lb->Debug_Level, lb->Proc, 
                       lb->Debug_Proc);

  /* Set oct_wgtflag based on the "key" parameter Obj_Weight_Dim */
  oct_wgtflag = (lb->Obj_Weight_Dim > 0);

  /* Initialization in case of early exit */
  *num_import = -1;
  *num_export = -1;  /* We don't compute any export data */

  /* Error checking for parameters */
  if (oct_dim < 2 || oct_dim > 3) {
    fprintf(stderr, "OCT Error in OCTPART: OCT_DIM must be 2 or 3\n");
    error = TRUE;
  }
  if (oct_method < 0 || oct_method > 2) {
    fprintf(stderr, "OCT Error in OCTPART: OCT_METHOD must be 0, 1, or 2\n");
    error = TRUE;
  }
  if (oct_maxoctregions < 1) {
    fprintf(stderr, "OCT Error in OCTPART: OCT_MAXOBJECTS "
                    "must be greater than 0\n");
    error = TRUE;
  }
  if (oct_minoctregions < 1) {
    fprintf(stderr, "OCT Error in OCTPART: OCT_MINOCTREGIONS "
                    "must be greater than 0\n");
    error = TRUE;
  }
  if (oct_output_level < 0 || oct_output_level > 2) {
    fprintf(stderr, "OCT Error in OCTPART: OCT_OUTPUT_LEVEL must be 0, 1, or 2\n");
    error = TRUE;
  }

  if (error)
    return(ZOLTAN_FATAL);
  else
    return(lb_oct_init(lb, num_import, import_global_ids, import_local_ids, 
                       import_procs, oct_dim, oct_method, oct_maxoctregions, 
                       oct_minoctregions, oct_output_level, oct_wgtflag));
}

/*****************************************************************************/
/*
 * void oct_init(LB *load_balancing_structure, int *number_of_objects,
 *               int *numer_of_old_objects, int *number_of_non_local_objects,
 *               LB_TAG **non_local_objects)
 *
 * initialize the calls needed to start the octree load balancing rounties
 */
static int lb_oct_init(
  LB *lb,                       /* The load-balancing structure with info for
                                   the OCTPART balancer.                     */
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

  int oct_dim,                  /* Dimension of method to be used (2D or 3D) */
  int oct_method,               /* Flag specifying curve to be used.         */
  int oct_maxoctregions,        /* max # of objects in leaves of octree.     */
  int oct_minoctregions,        /* min # of objects in leaves of octree.     */
  int oct_output_level,         /* Flag specifying amount of output.         */
  int oct_wgtflag               /* Flag specifying use of object weights.    */
) 
{
  char *yo = "lb_oct_init";
  OCT_Global_Info *OCT_info;
  pRList  RootList;               /* list of all local roots */
  pOctant RootOct;                /* root octree octant */
  int nsentags;                    /* number of tags being sent */
  pRegion import_regs;             /* */
  int nrectags;                    /* number of tags received */               
  int count, kk;
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

  ZOLTAN_LB_TRACE_ENTER(lb, yo);

  MPI_Barrier(lb->Communicator);
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

  count = nsentags = nrectags = 0;

  if(lb->Data_Structure == NULL) {
    OCT_info = LB_POct_init(lb, lb->Proc, oct_dim);
    LB_set_method(OCT_info, oct_method);
    LB_oct_set_maxregions(oct_maxoctregions);
    LB_oct_set_minregions(oct_minoctregions);
    createpartree = 1;
  }
  else {
    OCT_info = (OCT_Global_Info *) (lb->Data_Structure);
  }

  /* create the octree structure */
  time1 = MPI_Wtime();

  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_oct_gen_tree_from_input_data");
  LB_oct_gen_tree_from_input_data(lb, oct_wgtflag, &counters[1], &counters[2], 
			          &counters[3], &c[0], createpartree);

  time2 = MPI_Wtime();
  timers[0] = time2 - time1;                 /* time took to create octree */
/*   LB_POct_printResults(OCT_info); */
  /* partition the octree structure */
  time1 = MPI_Wtime();
  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_dfs_partition");
  LB_dfs_partition(lb, &counters[0], &c[1]);
  time2 = MPI_Wtime();
  timers[1] = time2 - time1;              /* time took to partition octree */

  /* set up tags for migrations */
  time1 = MPI_Wtime();

  count = 0; 
  RootList = LB_POct_localroots(OCT_info);
  while((RootOct = RL_nextRootOctant(&RootList))) {
    while(RootOct) {
      if(LB_Oct_isTerminal(RootOct)) {	
	count += LB_Oct_nRegions(RootOct);
      }
      RootOct = LB_POct_nextDfs(OCT_info, RootOct);
    }
  }

  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_dfs_migrate");
  LB_dfs_migrate(lb, &nsentags, &import_regs, &nrectags, 
	         &c[2], &c[3], &counters[3], &counters[5]);


  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_fix_tags");
  if (lb->Return_Lists) {
    *num_import = nrectags;
    if (nrectags > 0)
      LB_fix_tags(lb, import_global_ids, import_local_ids, import_procs,
                  nrectags, import_regs);
  }

  time2 = MPI_Wtime();
  timers[2] = time2 - time1;               /* time took to setup migration */


  /* count the number of objects on this processor */
  count = 0; 
  RootList = LB_POct_localroots(OCT_info);
  while((RootOct = RL_nextRootOctant(&RootList))) {
    while(RootOct) {
      if(LB_Oct_isTerminal(RootOct)) {	
	count += LB_Oct_nRegions(RootOct);
      }
      RootOct = LB_POct_nextDfs(OCT_info, RootOct);
    }
  }

  counters[4] = nsentags;
  MPI_Barrier(lb->Communicator);
  timestop = MPI_Wtime();

  if (oct_output_level > 0) {
    ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_oct_print_stats");
    LB_oct_print_stats(lb, timestop-timestart, timers, counters, c, 
                       oct_output_level);
  }

  for (kk = 0; kk < nrectags; kk++) {
    ZOLTAN_FREE(&(import_regs[kk].Global_ID));
    ZOLTAN_FREE(&(import_regs[kk].Local_ID));
  }
  ZOLTAN_FREE(&import_regs);

  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_oct_global_clear");
  LB_oct_global_clear(OCT_info);
  /* KDDKDD Don't understand how re-used octree will work, especially without
   * KDDKDD the LB_Bounds_Geom function.  For now, we'll delete everything;
   * KDDKDD we can move back to saving some of the tree later.
   */
  LB_OCT_Free_Structure(lb);
  /* KDDKDD END */

  /* Temporary return value until error codes are fully implemented. */
  ZOLTAN_LB_TRACE_EXIT(lb, yo);
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*
 * void LB_oct_gen_tree_from_input_data()
 *
 * This function will create a root node of on each processor which will
 * then be used to create an octree with regions associated with it. The
 * tree will then be balanced and the output used to balance "mesh regions"
 * on several processors.
 */
static void LB_oct_gen_tree_from_input_data(LB *lb, int oct_wgtflag, int *c1, 
                                     int *c2, int *c3, float *c0, int createpartree) 
{
  char *yo = "LB_oct_gen_tree_from_input_data";
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
          proc,           /* processor leaf node of parition tree belongs to */
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
  double bounds[6] = {MAXDOUBLE,MAXDOUBLE,MAXDOUBLE,-MAXDOUBLE,-MAXDOUBLE,-MAXDOUBLE};
  COORD global_min, global_max;
#endif /* KDDKDD_NEW_BOUNDS_GEOM_QUERY_FN */
  int nroots = 0;
  /*test*/
  /* COORD gmin,gmax; */

  OCT_Global_Info *OCT_info = (OCT_Global_Info *) (lb->Data_Structure);

  ZOLTAN_LB_TRACE_ENTER(lb, yo);
  /*
   * If there are no objects on this processor, do not create a root octant.
   * The partitioner will probably assign objects to this processor
   */
  if(lb->Get_Num_Obj == NULL) {
    fprintf(stderr, "OCT %s\n\t%s\n", "Error in octree load balance:",
	    "Must register Get_Num_Local_Objects function");
    abort();
  }
  *c3 = num_objs = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr) {
    fprintf(stderr, "OCT [%d] %s: Error returned from user defined "
                    "Get_Num_Obj function.\n", lb->Proc, yo);
    exit (-1);
  }
  ptr1 = NULL;

  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_get_bounds");
  /* Need A Function To Get The Bounds Of The Local Objects */
  LB_get_bounds(lb, &ptr1, &num_objs, min, max, oct_wgtflag, c0);
  
#ifndef KDDKDD_NEW_BOUNDS_GEOM_QUERY_FN
  /* For now, don't want to add the new query function to Zoltan. */
  /* LB_get_bounds appears to compute the global min and max from */
  /* the object input. */
  vector_set(OCT_info->OCT_gmin, min);
  vector_set(OCT_info->OCT_gmax, max);
#else
  /*test*/
  /*getMaxBounds(&gmin, &gmax);*/
  if(lb->Get_Bounds_Geom == NULL) {
    fprintf(stderr, "OCT %s\n\t%s\n", "Error in octree load balance:",
	    "Must register Get_Bounds_Geom function");
    abort();
  }
  lb->Get_Bounds_Geom(lb->Get_Bounds_Geom_Data, bounds, &ierr); 
  
  MPI_Allreduce(&(bounds[0]), &(global_min[0]), 3, 
		MPI_DOUBLE, MPI_MIN, lb->Communicator);
  MPI_Allreduce(&(bounds[3]), &(global_max[0]), 3,
		MPI_DOUBLE, MPI_MAX, lb->Communicator);
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

    level = 0;                                     /* initialize level count */

  /* 
   * if more than 1 processor, need to find what level of refinement needed
   * to initially partition bounding box among the processors 
   */

  
    if(lb->Num_Proc > 1) {
      n = lb->Num_Proc;
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
  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Before createpartree");

  if(createpartree) {
    /* create the global root octant */
    root = LB_POct_new(OCT_info);
    LB_Oct_setbounds(root, OCT_info->OCT_gmin, OCT_info->OCT_gmax);
    /* LB_Oct_setOrientation(root, 0); */
  
    /* subdivide to as many levels as calculated */
    for(i=0; i<level; i++) {
      cursor = root;
      while(cursor != NULL) {
	if(LB_Oct_isTerminal(cursor)) {
	  cursor2 = LB_POct_nextDfs(OCT_info, cursor);
	  LB_oct_terminal_refine(lb, cursor, 0);
	  cursor = cursor2;
	}
	else 
	  cursor = LB_POct_nextDfs(OCT_info, cursor);
      }
    }
    
#if 0
    if(lb->Proc == 0)
      for(i=0; i<8; i++)
	if(LB_Oct_child(root, i) == NULL)
	  fprintf(stderr,"NULL child pointer\n");
	else
	  fprintf(stderr, "child %d exists\n", i);
#endif

  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Before create map array");
    /* this part creates the map array */
    if(OCT_info->OCT_dimension == 2) {
      hold = (int)POW(4, level);                  /* ignoring the z+ octants */
      if(hold == 0)
	hold = 1;
    }
    else
      hold = (int)POW(8, level);

    part = hold / lb->Num_Proc;          /* how many octants per partition */
    remainder = hold % lb->Num_Proc; /* extra octants, not evenly divisible */
    extra = lb->Num_Proc - remainder;/* where to start adding extra octants */
    array = (Map *) ZOLTAN_MALLOC(hold * sizeof(Map));       /* allocate map array */
    if(array == NULL) {
      fprintf(stderr, "OCT ERROR on proc %d, could not allocate array map\n",
	      lb->Proc);
      abort();
    }
    /* initialize variables */
    proc = 0;
    count = 0;
    i = 0;
    cursor = root; 
    while(cursor != NULL) {
      cursor2 = LB_POct_nextDfs(OCT_info, cursor);
      if((LB_Oct_isTerminal(cursor)) && (i < hold)) {
	if(proc == extra) {
	  part++;
	  extra = -1;
	}
	if(count != part) {
	  array[i].npid = proc;
	  array[i].list = RL_initRootList();
	  LB_Oct_bounds(cursor, min, max);
	  vector_set(array[i].min, min);
	  vector_set(array[i].max, max);
	  count++;
	}
	else {
	  count = 1;
	  proc++;
	  array[i].npid = proc;
	  array[i].list = RL_initRootList();
	  LB_Oct_bounds(cursor, min, max);
	  vector_set(array[i].min, min);
	  vector_set(array[i].max, max);
	}
	if(proc == lb->Proc) {
	  array[i].npid = -1;
          /* KDDKDD Added RL_freeList below.  The 
           * KDDKDD implementation from RPI leaked memory because the 
           * KDDKDD test cases for setting array[i].list were not mutually 
           * KDDKDD exclusive.  Freeing the list produces the result we got
           * KDDKDD before, without the memory leak.
           */
          RL_freeList(&(array[i].list));
          /* KDDKDD End addition */
	  array[i].list = RL_initRootList();
	  parent = LB_Oct_parent(cursor);
	  if(parent != NULL)
	    LB_Oct_setchild(parent, cursor->which, NULL);
 	  LB_POct_setparent(OCT_info, cursor, NULL, -1);    /* octant into local root list */
	  LB_Oct_setMapIdx(cursor, i);
	  nroots++;
	  /*	  LB_POct_setparent(OCT_info, cursor, NULL, lb->Proc);     octant into local root list */
	}
	i++;
      }
      cursor = cursor2;
    } 
    RootList = LB_POct_localroots(OCT_info); 
    RootOct = RL_nextRootOctant(&RootList);
    if(RootOct != root) {
      /* KDDKDDFREE changed root to &root to allow root to be reset to NULL. */
      LB_POct_delTree(OCT_info,&root);
    }
    
    OCT_info->map = array;
    OCT_info->mapsize = hold;
  }

  /* 
   * attach the regions to the root... LB_oct_fix will create the octree
   * starting with the root and subdividing as needed 
   */    
  num_extra = LB_oct_fix(lb, ptr1, num_objs);
 
  ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Calling LB_migreg_migrate_orphans");
  LB_migreg_migrate_orphans(lb, ptr1, num_extra, level, OCT_info->map, c1, c2);

/*   ZOLTAN_FREE(&array); */
  while(ptr1 != NULL) {
    ptr = ptr1->next;
    ZOLTAN_FREE(&(ptr1->Global_ID));
    ZOLTAN_FREE(&(ptr1->Local_ID));
    ZOLTAN_FREE(&ptr1);
    ptr1 = ptr;
  }
  ZOLTAN_LB_TRACE_EXIT(lb, yo);
}

/*****************************************************************************/

static void LB_get_bounds(LB *lb, pRegion *ptr1, int *num_objs, 
		   COORD min, COORD max, int wgtflag, float *c0) 
{
  char *yo = "LB_get_bounds";
  ZOLTAN_ID_PTR obj_global_ids = NULL; 
  ZOLTAN_ID_PTR obj_local_ids = NULL;
  ZOLTAN_ID_PTR lid;       /* Temporary pointer to a local ID; used to pass NULL 
                          to query functions when NUM_LID_ENTRIES == 0. */
  ZOLTAN_ID_PTR next_lid;  /* Temporary pointer to a local ID; used to pass NULL 
                          to query functions when NUM_LID_ENTRIES == 0. */
  float *obj_wgts = NULL;
  int i, found;
  pRegion tmp, ptr;
  COORD global_min, global_max;
  double PADDING = 0.0000001;
  int ierr = 0;
  int num_gid_entries = lb->Num_GID;
  int num_lid_entries = lb->Num_LID;

  /* Initialization */
  max[0] = max[1] = max[2] = -MAXDOUBLE;
  min[0] = min[1] = min[2] =  MAXDOUBLE;

  *num_objs = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr) {
    fprintf(stderr, "OCT [%d] %s: Error returned from user defined "
                    "Get_Num_Obj function.\n", lb->Proc, yo);
    exit (-1);
  }

  if (*num_objs > 0) {
    obj_global_ids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb,(*num_objs));
    obj_local_ids  = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb,(*num_objs));
    obj_wgts       = (float *) ZOLTAN_MALLOC((*num_objs) * sizeof(float));
    if (!obj_global_ids || (num_lid_entries && !obj_local_ids) || !obj_wgts) {
      fprintf(stderr, "OCT [%d] Error from %s: Insufficient memory\n",lb->Proc,yo);
      exit(-1);
    }
    if (wgtflag == 0)
      for (i = 0; i < *num_objs; i++) obj_wgts[i] = 0.;


    if(lb->Get_Obj_List == NULL &&
      (lb->Get_First_Obj == NULL || lb->Get_Next_Obj == NULL)) {
      fprintf(stderr, "OCT Error in octree load balance:  user must declare " 
              "function Get_Obj_List or Get_First_Obj/Get_Next_Obj.");
      abort();
    }

    lid = (num_lid_entries ? &(obj_local_ids[0]) : NULL);
    if (lb->Get_Obj_List != NULL) {
      lb->Get_Obj_List(lb->Get_Obj_List_Data,
                       num_gid_entries, num_lid_entries,
                       obj_global_ids, obj_local_ids,
                       wgtflag, obj_wgts, &ierr);
      found = TRUE;
    }
    else {
      found = lb->Get_First_Obj(lb->Get_First_Obj_Data, 
                                num_gid_entries, num_lid_entries,
                                &(obj_global_ids[0]), lid, 
                                wgtflag, &(obj_wgts[0]), &ierr);
    }
    if (ierr) {
      fprintf(stderr, "OCT [%d] %s: Error returned from user defined "
                      "Get_Obj_List/Get_First_Obj function.\n", lb->Proc, yo);
      exit (-1);
    }

    if(found) {
      initialize_region(lb, &tmp, &(obj_global_ids[0]), lid,
                        wgtflag, obj_wgts[0]);
      *c0 = (float)tmp->Weight;
      vector_set(min, tmp->Coord);
      vector_set(max, tmp->Coord);
    }
    *ptr1 = tmp;
    for (i = 1; i < (*num_objs); i++) {
      if (num_lid_entries) {
        lid = &(obj_local_ids[(i-1)*num_lid_entries]);
        next_lid = &(obj_local_ids[i*num_lid_entries]);
      }
      else
        lid = next_lid = NULL;
      if (lb->Get_Obj_List == NULL) {
        found = lb->Get_Next_Obj(lb->Get_Next_Obj_Data, 
                                 num_gid_entries, num_lid_entries,
                                 &(obj_global_ids[(i-1)*num_gid_entries]),
                                 lid,
                                 &(obj_global_ids[i*num_gid_entries]),
                                 next_lid,
                                 wgtflag, &(obj_wgts[i]), &ierr);
        if (ierr) {
          fprintf(stderr, "OCT [%d] %s: Error returned from user defined "
                          "Get_Next_Obj function.\n", lb->Proc, yo);
          exit (-1);
        }
      }
      if (!found) {
        fprintf(stderr, "OCT Error in octree load balance:  number of objects "
               "declared by LB_NUM_OBJ_FN %d != number obtained by "
               "GET_NEXT_OBJ %d\n", *num_objs, i);
        exit(-1);
      }
      initialize_region(lb, &(ptr), &(obj_global_ids[i*num_gid_entries]), 
                        next_lid, wgtflag, obj_wgts[i]);
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
  
      ptr = NULL;
    }
    ZOLTAN_FREE(&obj_global_ids);
    ZOLTAN_FREE(&obj_local_ids);
    ZOLTAN_FREE(&obj_wgts);
  }
  
  MPI_Allreduce(&(min[0]), &(global_min[0]), 3, 
		MPI_DOUBLE, MPI_MIN, lb->Communicator);
  MPI_Allreduce(&(max[0]), &(global_max[0]), 3,
		MPI_DOUBLE, MPI_MAX, lb->Communicator);

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
static void initialize_region(LB *lb, pRegion *ret, ZOLTAN_ID_PTR global_id,
                              ZOLTAN_ID_PTR local_id, int wgtflag, float wgt) 
{
  pRegion reg;
  int ierr = 0;
  char *yo = "initialize_region";
  reg = (pRegion) ZOLTAN_MALLOC(sizeof(Region));
  *ret = reg;
  reg->Global_ID = ZOLTAN_ZOLTAN_MALLOC_GID(lb);
  reg->Local_ID = ZOLTAN_ZOLTAN_MALLOC_LID(lb);
  ZOLTAN_LB_SET_GID(lb, reg->Global_ID, global_id);
  ZOLTAN_LB_SET_LID(lb, reg->Local_ID, local_id);
  reg->Proc = lb->Proc;
  /* reg->Proc = 0; */
  reg->Coord[0] = reg->Coord[1] = reg->Coord[2] = 0.0;
  lb->Get_Geom(lb->Get_Geom_Data, lb->Num_GID, lb->Num_LID,
               global_id, local_id, reg->Coord, &ierr);
  if (ierr) {
    fprintf(stderr, "OCT [%d] %s: Error returned from user defined "
                    "Get_Geom function.\n", lb->Proc, yo);
    exit (-1);
  }

#if 0
  LB_Print_Sync_Start(lb->Communicator, TRUE);
    fprintf(stderr, "Result info on %d: %d %d %d  %lf  %lf  %lf\n", lb->Proc, 
	    reg->Local_ID, reg->Global_ID, reg->Proc,
	    reg->Coord[0], reg->Coord[1], reg->Coord[2]); 
  LB_Print_Sync_End(lb->Communicator, TRUE);
#endif

  if (wgtflag)
    reg->Weight = wgt;
  else
    reg->Weight = 1;

  if (ierr) {
    fprintf(stderr, "OCT [%d] %s: Error returned from user defined "
                    "Get_Obj_Weight function.\n", lb->Proc, yo);
    exit (-1);
  }

  reg->next = NULL;
}

/*****************************************************************************/
/*
 * int LB_oct_fix(pMesh mesh)
 *
 * Clear all region pointers from the octree,
 * reinsert all mesh regions, 
 * and coarsen or refine the octree as needed.
 * 
 * Return the number of regions that could
 * not be inserted.
 *
 */
static int LB_oct_fix(LB *lb, pRegion Region_list, int num_objs) 
{
  int nreg;                                /* number of regions not inserted */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(lb->Data_Structure);

  /* initalize variables */
  nreg=0;
  oct_nref=0;
  oct_ncoarse=0;
  /* associate the objects to the octants */
  nreg=LB_oct_global_insert_object(lb, Region_list, num_objs);
  LB_oct_global_dref(lb, OCT_info); 

  return(nreg);
}

/*****************************************************************************/
/*
 * int LB_oct_global_insert_object()
 *
 * Try to insert all regions in the mesh into the existing
 * local octree.  Return the number of insertion failures.
 *
 */
static int LB_oct_global_insert_object(LB *lb, pRegion Region_list, int num_objs) 
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

    if (!LB_oct_global_insert(lb, region)) {
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

/*****************************************************************************/
/*
 * LB_oct_global_insert(region)
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
static pOctant LB_oct_global_insert(LB *lb, pRegion region) 
{
  pOctant oct;                            /* octree octant */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(lb->Data_Structure);

  oct = NULL;
  /* find the octant which the object lies in */
  oct=LB_oct_global_find(OCT_info,region->Coord);
  /* if octant is found, try to insert the region */
  if (oct)
    if (!LB_oct_subtree_insert(lb, oct, region))         /* inserting region */
      {
	fprintf(stderr,"OCT LB_oct_global_insert: insertion failed\n");
	abort();
      }

  return(oct);
}

/*****************************************************************************/
/*
 * LB_oct_subtree_insert(oct,region)
 *
 * Insert region in oct, carrying out multiple refinement if necessary
 */
int LB_oct_subtree_insert(LB *lb, pOctant oct, pRegion region) 
{
OCT_Global_Info *OCT_info = (OCT_Global_Info *)(lb->Data_Structure); 
  /* if oct is not terminal, find leaf node the centroid can be attahced to */
  if (!LB_Oct_isTerminal(oct)) {
    oct=LB_oct_findOctant(OCT_info, oct,region->Coord);
  }

  if (!oct)
    return(0);

  /* add the region to the octant */
  LB_Oct_addRegion(lb, oct, region);
 
  /* check if octant has too many regions and needs to be refined */
  /* KDDKDD  Replaced the following to allow multiple regions with the 
   * KDDKDD  same coordinates to be placed in a single octant.
  if(LB_Oct_nRegions(oct) > MAXOCTREGIONS)
   */
  if(LB_Oct_nUniqueRegions(OCT_info, oct) > MAXOCTREGIONS)
    LB_oct_terminal_refine(lb, oct,0);    /* After this, dest may be nonterm */

  return(1);
}

/*****************************************************************************/
/*
 * pOctant LB_oct_global_find(COORD centroid)
 *
 * Return the octant in the local octree that contains the
 * given point, if it exists.  Otherwise, return NULL.
 *
 */
static pOctant LB_oct_global_find(OCT_Global_Info *OCT_info,COORD point) 
{  
  pRList  RootList;               /* list of all local roots */
  pOctant RootOct;                /* root octree octant */
  pOctant oct;                    /* octree octant */
  COORD min,                      /* minimum bounds coordinates */
        max;                      /* maximum bounds coordinates */

  RootList = LB_POct_localroots(OCT_info);                   /* get a list of all local roots */
  oct = NULL;

  /* iterate through root list to find if point lies inside root's bounds */
  while ((!oct) && (RootOct = RL_nextRootOctant(&RootList)) ) {
    LB_Oct_bounds(RootOct,min,max);
    /* ATTN: LB_in_box may need adjusting for the lower bounds */
    /* check if point fits inside */
    if (LB_in_box(OCT_info,point,min,max)) {  
      /* find the exact octant for point */
      oct=LB_oct_findOctant(OCT_info,RootOct,point);
    }
  }

  return(oct);
}

/*****************************************************************************/
/*
 * LB_oct_findOctant(oct, coord)
 *   (replaces : PO_findOctant(oct,coord))
 *  
 *
 * find the octant in a subtree containing coord (if coord
 * is not in the subtree, returns the closest octant in the subtree).
 * NOTE: return NULL if we hit an off-processor link
 *
 */
static pOctant LB_oct_findOctant(OCT_Global_Info *OCT_info,pOctant oct, COORD coord) 
{
  pOctant child;                                  /* child of an octant */
  int cnum;                                       /* child number */
  
  /* if octant is terminal, then this is the right octant */
  if (LB_Oct_isTerminal(oct))
    return(oct);

  cnum = LB_child_which_wrapper(OCT_info,oct,coord);   /* find closest child to coord */
  child = LB_Oct_child(oct, cnum);                            /* get that child */
  /* ATTN: are these local checks necessary? */
  /* if ( !LB_POct_local(child) ) */
  if(!LB_POct_local(OCT_info,oct, cnum))                     /* make sure octant is local */
    return(NULL);

  return(LB_oct_findOctant(OCT_info,child,coord)); /* recursivly search down the tree */
}

/*****************************************************************************/
/*
 * LB_oct_terminal_refine(oct)
 *
 * subdivide a terminal octant and divvy up
 * its regions to the 8 children; recurse if
 * necessary to satisfy MAXOCTREGIONS
 *
 */
static void LB_oct_terminal_refine(LB *lb, pOctant oct,int count) 
{
  COORD min,                      /* coordinates of minimum bounds of region */
        max,                      /* coordinates of maximum bounds of region */
        origin;                   /* origin of region */
  pOctant child[8];               /* array of child octants */
  int cnum;                       /* child number */
  int i;                       /* index counter */
  pRegion region;                 /* a region to be associated to an octant */
  pRegion entry;
  COORD cmin[8], cmax[8];
  OCT_Global_Info *OCT_info = (OCT_Global_Info *) (lb->Data_Structure);

  for(i=0;i<3;i++)
    min[i] = max[i] = 0;
  
  /* upper limit of refinement levels */
  /* ATTN: may not be used anymore, but can be put back in if necessary */
  if (count>=20) {
    fprintf(stderr, "OCT ERROR: LB_oct_terminal_refine: bailing out at "
                    "10 levels\n");
    abort();
  }
  oct_nref++;                                /* increment refinement counter */

  /* octant should be terminal in order to be refined (subdivided) */
  if (!LB_Oct_isTerminal(oct)) {
    fprintf(stderr,"OCT ref_octant: oct not terminal\n");
    abort();
  }

  /* get the bounds of an octant */
  LB_Oct_bounds(oct,min,max);
  /* calculate the origin from the bounds */
  LB_bounds_to_origin(min,max,origin);

  region = LB_Oct_regionlist(oct);             /* Get list while still terminal */
  oct->list = NULL;  /* remove regions from octant, it won't be terminal */

  /* create the children and set their id's */

  LB_child_bounds_wrapper(OCT_info,oct, cmin, cmax);
  for (i=0; i<8; i++) {
    if(OCT_info->OCT_dimension == 2) {
      /* KDDKDD 3/01 see changes to LB_child_bounds_wrapper that allow this
       * KDDKDD 3/01 test to work for GRAY and HILBERT mappings.
       */
      if(cmin[i][2] > OCT_info->OCT_gmin[2]) {                /* ignore the z+ octants */
	child[i] = NULL;
	continue;
      }
    }
    
    child[i]=LB_POct_new(OCT_info);                       /* create a new octant */
    child[i]->dir = LB_get_child_dir(OCT_info, oct->dir, i);
                  /* create a new octant */
    LB_POct_setparent(OCT_info, child[i], oct, lb->Proc);   /* set the child->parent link */
    LB_Oct_setchildnum(child[i], i);              /* which child of the parent */
    LB_Oct_setchild(oct, i, child[i]);            /* set the parent->child link */
#ifdef LGG_MIGOCT
    LB_Oct_setID(child[i], LB_oct_nextId());             /* set child id number */
#endif /* LGG_MIGOCT */
    LB_Oct_setbounds(child[i], cmin[i], cmax[i]);           /* set child bounds */
    LB_Oct_setCpid(oct, i, lb->Proc);      /* set child to be a local octant */
    /*    LB_Oct_setOrientation(child[i], 
		       LB_child_orientation(oct->orientation, oct->which));
		       */
  }

  /* assign newly created children to child array*/
  if(OCT_info->OCT_dimension == 3) {
    if(LB_Oct_children(oct, child) != 8) {
      /* 
       * if subdivision of oct was successful, oct should have 8 children; 
       * thus a return value of 0 here is a fatal error
       */
      fprintf(stderr, "OCT ref_octant: subdivide failed, %d children.\n",
	      LB_Oct_children(oct, child));
      abort();
    }
  }
  else
    if(LB_Oct_children(oct, child) != 4) {
      /* 
       * if subdivision of oct was successful, oct should have 4 children; 
       * thus a return value of 0 here is a fatal error
       */
      fprintf(stderr, "OCT ref_octant:subdivide failed, %d children, expected 4\n",
	      LB_Oct_children(oct, child));
      abort();
    }

  /* iterate through and find which child each region should belong to */
  while(region != NULL) {
    entry = region->next;
    cnum=LB_child_which_wrapper(OCT_info,oct, region->Coord);
    /* add region to octant's regionlist */
    LB_Oct_addRegion(lb, child[cnum], region);
    ZOLTAN_FREE(&(region->Global_ID));
    ZOLTAN_FREE(&(region->Local_ID));
    ZOLTAN_FREE(&region);
    region = entry;
  }

  for (i=0; i<8; i++)                                           /* Recursion */
    if(child[i] != NULL)
      /* KDDKDD  Replaced the following to allow multiple regions with the
       * KDDKDD same coordinates to be placed in the same octant.
      if (LB_Oct_nRegions(child[i]) > MAXOCTREGIONS) {
       */
      if (LB_Oct_nUniqueRegions(OCT_info,child[i]) > MAXOCTREGIONS) {
	LB_oct_terminal_refine(lb, child[i],count+1);
      }
}

/*****************************************************************************/
/*
 * LB_oct_global_dref()
 * 
 * refine and derefine octree as necessary by number of
 * regions in each octant
 *
 */
static void LB_oct_global_dref(LB *lb, OCT_Global_Info *OCT_info) 
{
  pRList  RootList;                           /* list of all local roots */
  pOctant RootOct;

  RootList = LB_POct_localroots(OCT_info);
  while ((RootOct = RL_nextRootOctant(&RootList))) 
    LB_oct_subtree_dref(lb, OCT_info, RootOct);
}

/*****************************************************************************/
/*
 * LB_oct_subtree_dref(oct)
 *
 * Coarsen octree so that leaf octants do not have less than MINOCTREGIONS 
 * regions. Refinement takes precedence, so coarsening will not take place 
 * unless all subtrees agree on it.
 */
static int LB_oct_subtree_dref(LB *lb, OCT_Global_Info *OCT_info,pOctant oct) 
{
  pOctant child;                         /* child of an octant */
  int coarsen;                           /* flag to indicate need to coarsen */
  int i;                                 /* index counter */
  int nregions;                          /* number of regions */
  int total;                             /* total number of regions */

  /* if terminal octant, cannot coarsen on own */
  if (LB_Oct_isTerminal(oct)) {
    nregions=LB_Oct_nRegions(oct);

    if (nregions > (MAXOCTREGIONS*2) ) {
      fprintf(stderr, "OCT LB_oct_subtree_dref: warning: too many (%d) regions "
	     "in oct (id=%d)\n",LB_Oct_nRegions(oct),LB_Oct_id(oct));
      return(-1);
    }
    else
      if (nregions < MINOCTREGIONS)
	return(nregions);
      else
	return(-1);
  }

  coarsen=1;                                         /* assume to be coarsen */
  total=0;

  /* look at each child, see if they need to be coarsened */
  for (i=0; coarsen && i<8; i++) {
    child = LB_Oct_child(oct,i);
    /* if child is off processor cannot coarsen */
    /* if (!LB_Oct_local(child) ) */
    if(!LB_POct_local(OCT_info, oct, i))
      coarsen=0;
    else {
      /* get the number of region of the child */
      nregions=LB_oct_subtree_dref(lb, OCT_info,child);
      if (nregions<0)
	coarsen=0;
      else
	total+=nregions;                           /* total the region count */
    }
  }

  /* check if octant can be coarsened */
  if (coarsen && total<MAXOCTREGIONS) {
    LB_oct_terminal_coarsen(lb, OCT_info,oct);
    
    if (total < MINOCTREGIONS)
      return(total);
    else
      return(-1);
  }

  return(-1);
}

/*****************************************************************************/
/*
 * LB_oct_terminal_coarsen(oct)
 *
 * remove octant's children, accumulating regions
 * to octant
 *
 */
static void LB_oct_terminal_coarsen(LB *lb, OCT_Global_Info *OCT_info,pOctant oct) 
{
  pOctant child;                         /* child of an octant */
  pRegion region;                        /* region associated with an octant */
  int i;                                 /* index counter */
  pRegion regionlist[8];                 /* an array of region lists */

  oct_ncoarse++;                             /* increment coarsening counter */

  for(i=0; i<8; i++) {
    /* get the ith child of an octant */
    child = LB_Oct_child(oct,i);
    
    /* cannot coarsen if child is off-processor */
    /* if(!LB_POct_local(child)) X */
    if(!LB_POct_local(OCT_info, oct, i)) {          /* cannot be off-processor */
      fprintf(stderr,"OCT LB_oct_terminal_coarsen: child not local\n");
      abort();
    }
    
    if(!LB_Oct_isTerminal(child)) {
      fprintf(stderr,"OCT LB_oct_terminal_coarsen: child not terminal\n");
      abort();
    }
    
    /* get each child's region list */
    regionlist[i] = LB_Oct_regionlist(child);
    
    /* delete each child */
    /* KDDKDDFREE Change child to &child. */
    LB_POct_free(OCT_info, &child);
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
      LB_Oct_addRegion(lb, oct,region);           
      region = region->next;
    }
  }
}

/*****************************************************************************/
static void LB_oct_set_maxregions(int max)
{ 
  if (max < 1) {
    fprintf(stderr, "OCT Warning LB_oct_set_maxregions(): %s\n",
	    "illegal input, using default.");
    MAXOCTREGIONS = 1;
  }
  else
    MAXOCTREGIONS = max; 
}

/*****************************************************************************/
static void LB_oct_set_minregions(int min)
{ 
  if (min < 1) {
    fprintf(stderr, "OCT Warning LB_oct_set_minregions(): %s\n",
	    "illegal input, using default.");
    MINOCTREGIONS = 1;
  }
  else
    MINOCTREGIONS = min; 
}

/*****************************************************************************/
#ifdef LGG_MIGOCT
void LB_oct_resetIdCount(int start_count)
{
  IDcount = start_count;
}

/*****************************************************************************/

int LB_oct_nextId(void)
{
  return ++IDcount;
}

/*****************************************************************************/
typedef struct
{
  pOctant ptr;
  int id;
} Rootid;

/*****************************************************************************/

static int idcompare(Rootid *i, Rootid *j)
{
  return( (i->id) - (j->id) );
}
/*****************************************************************************/
#endif /* LGG_MIGOCT */

/*****************************************************************************/


/*
 * oct_global_clear()
 *
 * delete all regions from all octants on the local processor
 *
 */
static void LB_oct_global_clear(OCT_Global_Info * OCT_info)
{ 
  pRList  RootList;                 /* list of all local roots */
  pOctant Oct;

  /* 
   * iterate through the list of local roots 
   * traverse down the subtree 
   * delete regions associated with octant 
   */
  RootList = LB_POct_localroots(OCT_info);
  while ((Oct = RL_nextRootOctant(&RootList))) {
    while(Oct) {
      if(LB_Oct_isTerminal(Oct))
	LB_Oct_clearRegions(Oct);
      Oct = LB_POct_nextDfs(OCT_info, Oct);
    }
  }
}

#if 0
/*****************************************************************************/

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


  if (oct=LB_POct_nextDfs(oct, 0))
    return(oct);

  return(PList_next(OCT_info->OCT_rootlist,temp));
}
#endif  /* 0 */

/*****************************************************************************/
/*
 * int LB_Oct_nUniqueRegions(pOctant octant)
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
int LB_Oct_nUniqueRegions(OCT_Global_Info *OCT_info, pOctant oct) {
  pRegion ptr;                     /* pointer to iterate through region list */
  int count;                       /* count of number of regions in list */
  pRegion ptr2;                    /* pointer to iterate through prev regions*/
  int same_coords;

  if (oct == NULL) 
    return 0;

  if(!LB_Oct_isTerminal(oct)) 
    return 0;

  count = 0;
  ptr = oct->list;
  while(ptr != NULL) {
    ptr2 = oct->list;
    same_coords = FALSE;
    while (ptr2 != ptr) {
      if (LB_Oct_CompareCoords(OCT_info->OCT_dimension, 
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
/*****************************************************************************/
/*
 * Routine to compare coordinates
 */
static int LB_Oct_CompareCoords(int dim, COORD pt1, COORD pt2) 
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
