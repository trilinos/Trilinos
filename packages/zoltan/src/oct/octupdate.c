/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_octupdatec_id = "$Id$";
#endif

#include "lb_const.h"
#include "octant_const.h"
#include "costs_const.h"
#include "util_const.h"
#include "dfs_const.h"
#include "octupdate.h"
#include "octupdate_const.h"
#include "migreg_const.h"
#include "migoct_const.h"
#include "all_allo_const.h"

/***************************  PROTOTYPES *************************************/
static void initialize_region(LB *, pRegion *, LB_GID, LB_LID, int, float);

/*****************************************************************************/
/* NOTE: be careful later about region lists for nonterminal octants */

static int oct_nref=0;                              /* number of refinements */
static int oct_ncoarse=0;                           /* number of coarsenings */
static int IDcount = 0;                            /* renumbering of octants */
static int MAXOCTREGIONS = 1;


/*****************************************************************************/
/*
 * void oct_init(LB *load_balancing_structure, int *number_of_objects,
 *               int *numer_of_old_objects, int *number_of_non_local_objects,
 *               LB_TAG **non_local_objects)
 *
 * initialize the calls needed to start the octree load balancing rounties
 */
int lb_oct_init(
  LB *lb,                     /* The load-balancing structure with info for
                                 the OCTPART balancer.                       */
  int *num_import,            /* Number of non-local objects assigned to this
                                 processor in the new decomposition.         */
  LB_GID **import_global_ids, /* Returned value:  array of global IDs for
                                 non-local objects in this processor's new
                                 decomposition.                              */
  LB_LID **import_local_ids,  /* Returned value:  array of local IDs for
                                 non-local objects in this processor's new
                                 decomposition.                              */
  int **import_procs          /* Returned value:  array of processor IDs for
                                 processors owning the non-local objects in
                                 this processor's new decomposition.         */
) 
{
  LB_TAG *export_tags;             /* array of LB_TAGS being exported */
  pRegion export_regs;             /* */
  int nsentags;                    /* number of tags being sent */
  LB_TAG *import_tags;             /* array of LB_TAGS being imported */
  pRegion import_regs;             /* */
  int nrectags;                    /* number of tags received */
  pOctant ptr;                     /* pointer to an octant */
  pRList root;                     
  pRList root2;                     
  int count, i;
  double time1,time2,time3,time4;  /* timers */
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

  export_tags = NULL;
  import_tags = NULL;
  count = nsentags = nrectags = 0;

  if(lb->Params != NULL) {
/* BAH: need to put new params stuff into octupdate */
/*
    if(lb->Params[0] == LB_PARAMS_INIT_VALUE)
      POC_init(lb->Proc, 3);
    else
      POC_init(lb->Proc, lb->Params[0]);

    if(lb->Params[1] == LB_PARAMS_INIT_VALUE)
      LB_set_method(0);
    else
      LB_set_method(lb->Params[1]);

    if(lb->Params[2] == LB_PARAMS_INIT_VALUE)
      LB_oct_set_maxregions(1);
    else
      LB_oct_set_maxregions(lb->Params[2]);
*/
  }
  else {
    POC_init(lb->Proc, 3);
    LB_set_method(0);
    LB_oct_set_maxregions(1);
  }

  /* create the octree structure */
  time1 = MPI_Wtime();

  LB_oct_gen_tree_from_input_data(lb, &counters[1], &counters[2], 
			          &counters[3], &c[0]);
  time2 = MPI_Wtime();
  timers[0] = time2 - time1;                 /* time took to create octree */
  
  /* partition the octree structure */
  time1 = MPI_Wtime();
  LB_dfs_partition(lb, &counters[0], &c[1]);
  time2 = MPI_Wtime();
  timers[1] = time2 - time1;              /* time took to partition octree */

  /* intermediate result print out */
/* BAH: need to put new params stuff into octupdate */
/*if(lb->Params[4] != LB_PARAMS_INIT_VALUE)*/          /* WARNING BIG OUTPUT!! */
  if (0) {
    for(i=0; i<lb->Num_Proc; i++) {
      if(lb->Proc == i) {
	POC_printResults();
	printf("\n\n");
      }
      MPI_Barrier(lb->Communicator);
    }
  }	

  /* set up tags for migrations */
  time1 = MPI_Wtime();
  LB_dfs_migrate(lb, &export_regs, &nsentags, &import_regs, &nrectags, 
	         &c[2], &c[3], &counters[3], &counters[5]);

  *num_import = nrectags;
  LB_fix_tags(import_global_ids, import_local_ids, import_procs, nrectags, 
	      import_regs);

  time2 = MPI_Wtime();
  timers[2] = time2 - time1;               /* time took to setup migration */


  /* count the number of objects on this processor */
  root = OCT_rootlist;
  while(root != NULL) {
    ptr = root->oct;
    while(ptr != NULL) {
      if(POC_isTerminal(ptr))
	count += POC_nRegions(ptr);
      ptr = POC_nextDfs(ptr);
    }
    root = root->next;
  }

  counters[4] = nsentags;
  MPI_Barrier(lb->Communicator);
  timestop = MPI_Wtime();

  if(lb->Params != NULL) {
/* BAH: need to put new params stuff into octupdate */
/*
    if(lb->Params[3] == LB_PARAMS_INIT_VALUE)
      LB_print_stats(lb, timestop - timestart, timers, counters, c, 1);
    else
      if(lb->Params[3] != 0)
	LB_print_stats(lb, timestop-timestart, timers, counters, c,
                       lb->Params[3]);
*/
  }
  else
    LB_print_stats(lb, timestop - timestart, timers, counters, c, 1);

  LB_Free((void **) &export_regs);
  LB_Free((void **) &import_regs);
  LB_Free((void **) &export_tags);
  root = POC_localroots();
  while(root != NULL) {
    root2 = root->next;
    POC_delTree(root->oct);
    root = root2;
  }

  /* Temporary return value until error codes are fully implemented. */
  return(DLB_OK);
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
void LB_oct_gen_tree_from_input_data(LB *lb, int *c1, int *c2, 
				     int *c3, float *c0) 
{
  char *yo = "LB_oct_gen_tree_from_input_data";
  COORD min,              /* min coord bounds of objects */
        max;              /* max coord bounds of objects */
  int num_extra;          /* number of orphaned objects */
  int num_objs;           /* total number of local objects */
  pRegion Region_list;    /* array of objects to be inserted */
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

  /*
   * If there are no objects on this processor, do not create a root octant.
   * The partitioner will probably assign objects to this processor
   */
  if(lb->Get_Num_Obj == NULL) {
    fprintf(stderr, "%s\n\t%s\n", "Error in octree load balance:",
	    "Must register Get_Num_Local_Objects function");
    abort();
  }
  *c3 = num_objs = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr) {
    fprintf(stderr, "[%d] %s: Error returned from user defined "
                    "Get_Num_Obj function.\n", lb->Proc, yo);
    exit (-1);
  }
  Region_list = NULL;
  ptr1 = NULL;
  if(num_objs > 0) {
    /* Need A Function To Get The Bounds Of The Local Objects */
    LB_get_bounds(lb, &ptr1, &num_objs, min, max, c0);

    vector_set(OCT_gmin, min);
    vector_set(OCT_gmax, max);
    
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
      if(OCT_dimension == 2)
	hold = 4;
      else
	hold = 8;
      remainder = hold;
      for(; remainder >=hold; level++) {
	int pr;

	pr = (int)pow(hold, level);
	remainder = n - pr;
	/*
	f = 0printf(stderr,"n = %d, pow = %d, level = %d, remainder = %d\n", 
		n, pr, level, remainder);
	*/
      }
      if(remainder == 0)
	level--;
    }

    /* create the global root octant */
    root = POC_new();
    POC_setbounds(root, OCT_gmin, OCT_gmax);
    /* POC_setOrientation(root, 0); */
    
    /* subdivide to as many levels as calculated */
    for(i=0; i<level; i++) {
      cursor = root;
      while(cursor != NULL) {
	if(POC_isTerminal(cursor)) {
	  cursor2 = POC_nextDfs(cursor);
	  LB_oct_terminal_refine(lb, cursor, 0);
	  cursor = cursor2;
	}
	else 
	  cursor = POC_nextDfs(cursor);
      }
    }
    
#if 0
    if(lb->Proc == 0)
      for(i=0; i<8; i++)
	if(POC_child(root, i) == NULL)
	  fprintf(stderr,"NULL child pointer\n");
	else
	  fprintf(stderr, "child %d exists\n", i);
#endif

    /* this part creates the map array */
    if(OCT_dimension == 2) {
      hold = (int)pow(4, level);                  /* ignoring the z+ octants */
      if(hold == 0)
	hold = 1;
    }
    else
      hold = (int)pow(8, level);

    part = hold / lb->Num_Proc;          /* how many octants per partition */
    remainder = hold % lb->Num_Proc; /* extra octants, not evenly divisible */
    extra = lb->Num_Proc - remainder;/* where to start adding extra octants */
    array = (Map *) LB_Array_Alloc(__FILE__, __LINE__, 1, hold,
                                   sizeof(Map));       /* allocate map array */
    if(array == NULL) {
      fprintf(stderr, "ERROR on proc %d, could not allocate array map\n",
	      lb->Proc);
      abort();
    }
    /* initialize variables */
    proc = 0;
    count = 0;
    i = 0;
    cursor = root;
    while(cursor != NULL) {
      cursor2 = POC_nextDfs(cursor);
      if((POC_isTerminal(cursor)) && (i < hold)) {
	if(proc == extra) {
	  part++;
	  extra = -1;
	}
	if(count != part) {
	  array[i].npid = proc;
	  POC_bounds(cursor, min, max);
	  vector_set(array[i].min, min);
	  vector_set(array[i].max, max);
	  count++;
	}
	else {
	  count = 1;
	  proc++;
	  array[i].npid = proc;
	  POC_bounds(cursor, min, max);
	  vector_set(array[i].min, min);
	  vector_set(array[i].max, max);
	}
	if(proc == lb->Proc) {
	  array[i].npid = -1;
	  parent = POC_parent(cursor);
	  if(parent != NULL)
	    POC_setchild(parent, cursor->which, NULL);
	  POC_setparent(cursor, NULL, -1);    /* octant into local root list */
	}
	i++;
      }
      cursor = cursor2;
    } 

    if(OCT_rootlist->oct != root)
      POC_delTree(root);
  }

  /*
   *  if(lb->Proc == 0)
   *    for(i=0; i<hold; i++)
   *      fprintf(stderr,"(%d) %lf %lf %lf, %lf %lf %lf\n", array[i].npid,
   *	      array[i].min[0], array[i].min[1], array[i].min[2],
   *	      array[i].max[0], array[i].max[1], array[i].max[2]);
   *  msg_sync();
   */

  /* 
   * attach the regions to the root... LB_oct_fix will create the octree
   * starting with the root and subdividing as needed 
   */  
  num_extra = LB_oct_fix(lb, ptr1, num_objs);

/* 
 * fprintf(stderr,"(%d) number of extra regions %d\n", lb->Proc, num_extra);
 */
  LB_migreg_migrate_orphans(lb, ptr1, num_extra, level, array, c1, c2);
  
  LB_Free((void **) &array);
  while(ptr1 != NULL) {
    ptr = ptr1->next;
    LB_Free((void **) &ptr1);
    ptr1 = ptr;
  }
}

/*****************************************************************************/

void LB_get_bounds(LB *lb, pRegion *ptr1, int *num_objs, 
		   COORD min, COORD max, float *c0) 
{
  char *yo = "LB_get_bounds";
  int max_num_objs;
  LB_GID *obj_global_ids; 
  LB_LID *obj_local_ids;
  float *obj_wgts;
  int i, found, wgtflag = 0;
  pRegion tmp, ptr;
  COORD global_min, global_max;
  double x;
  double PADDING = 0.0000001;
  int ierr = 0;

  /* ATTN: wgtflag should be determined by a user-defined option,
     probably set by LB_Set_Param. For now, read as input parameter. */
  if (lb->Params != NULL) {
/* BAH: need to put new params stuff into octupdate */
/*
    if (lb->Params[LB_PARAMS_MAX_SIZE-1] == LB_PARAMS_INIT_VALUE)
      wgtflag = 0;
    else
      wgtflag = lb->Params[LB_PARAMS_MAX_SIZE-1];
*/
  }

  *num_objs = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr) {
    fprintf(stderr, "[%d] %s: Error returned from user defined "
                    "Get_Num_Obj function.\n", lb->Proc, yo);
    exit (-1);
  }

  /* ATTN: an arbitrary choice, is this necessary? */
  max_num_objs = 2 * (*num_objs); 

  obj_global_ids = (LB_GID *) LB_Array_Alloc(__FILE__, __LINE__,
                                            1, *num_objs, sizeof(LB_GID));
  obj_local_ids  = (LB_LID *) LB_Array_Alloc(__FILE__, __LINE__,
                                            1, *num_objs, sizeof(LB_LID));
  obj_wgts       = (float *) LB_Array_Alloc(__FILE__, __LINE__,
                                            1, *num_objs, sizeof(float));
  if (!obj_global_ids || !obj_local_ids || !obj_wgts) {
    fprintf(stderr, "[%d] Error from %s: Insufficient memory\n", lb->Proc, yo);
    exit(-1);
  }


  if(lb->Get_Obj_List == NULL &&
    (lb->Get_First_Obj == NULL || lb->Get_Next_Obj == NULL)) {
    fprintf(stderr, "Error in octree load balance:  user must declare " 
            "function Get_Obj_List or Get_First_Obj/Get_Next_Obj.");
    abort();
  }

  if (lb->Get_Obj_List != NULL) {
    lb->Get_Obj_List(lb->Get_Obj_List_Data, obj_global_ids, obj_local_ids,
                     wgtflag, obj_wgts, &ierr);
    found = TRUE;
  }
  else {
    found = lb->Get_First_Obj(lb->Get_First_Obj_Data, &(obj_global_ids[0]),
                              &(obj_local_ids[0]), wgtflag, &(obj_wgts[0]), 
                              &ierr);
  }
  if (ierr) {
    fprintf(stderr, "[%d] %s: Error returned from user defined "
                    "Get_Obj_List/Get_First_Obj function.\n", lb->Proc, yo);
    exit (-1);
  }

  if(*num_objs > 0 && found) {
    initialize_region(lb, &tmp, obj_global_ids[0], obj_local_ids[0],
                      wgtflag, obj_wgts[0]);
    *c0 = (float)tmp->Weight;
    vector_set(min, tmp->Coord);
    vector_set(max, tmp->Coord);
  }
  *ptr1 = tmp;
  for (i = 1; i < (*num_objs); i++) {
    if (lb->Get_Obj_List == NULL) {
      found = lb->Get_Next_Obj(lb->Get_Next_Obj_Data, obj_global_ids[i-1],
                               obj_local_ids[i-1], &(obj_global_ids[i]),
                               &(obj_local_ids[i]), wgtflag, &(obj_wgts[i]),
                               &ierr);
      if (ierr) {
        fprintf(stderr, "[%d] %s: Error returned from user defined "
                        "Get_Next_Obj function.\n", lb->Proc, yo);
        exit (-1);
      }
    }
    if (!found) {
      fprintf(stderr, "Error in octree load balance:  number of objects "
             "declared by LB_NUM_OBJ_FN %d != number obtained by "
             "GET_NEXT_OBJ %d\n", *num_objs, i);
      exit(-1);
    }
    initialize_region(lb, &(ptr), obj_global_ids[i], obj_local_ids[i],
                      wgtflag, obj_wgts[i]);
    *c0 += (float)tmp->Weight;
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
  LB_Free((void **) &obj_global_ids);
  LB_Free((void **) &obj_local_ids);
  LB_Free((void **) &obj_wgts);
  
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
static void initialize_region(LB *lb, pRegion *ret, LB_GID global_id,
                              LB_LID local_id, int wgtflag, float wgt) 
{
  pRegion reg;
  int ierr = 0;
  char *yo = "initialize_region";

  reg = (pRegion) LB_MALLOC(sizeof(Region));
  *ret = reg;
  LB_SET_GID(reg->Tag.Global_ID, global_id);
  LB_SET_LID(reg->Tag.Local_ID, local_id);
  reg->Tag.Proc = lb->Proc;
  /* reg->Proc = 0; */
  reg->Coord[0] = reg->Coord[1] = reg->Coord[2] = 0.0;
  lb->Get_Geom(lb->Get_Geom_Data, global_id, local_id, reg->Coord, &ierr);
  if (ierr) {
    fprintf(stderr, "[%d] %s: Error returned from user defined "
                    "Get_Geom function.\n", lb->Proc, yo);
    exit (-1);
  }

#if 0
  LB_print_sync_start(lb, TRUE);
    fprintf(stderr, "Result info on %d: %d %d %d  %lf  %lf  %lf\n", lb->Proc, 
	    reg->Tag.Local_ID, reg->Tag.Global_ID, reg->Tag.Proc,
	    reg->Coord[0], reg->Coord[1], reg->Coord[2]); 
  LB_print_sync_end(lb, TRUE);
#endif

  if (wgtflag)
    reg->Weight = wgt;
  else
    reg->Weight = 1;

  if (ierr) {
    fprintf(stderr, "[%d] %s: Error returned from user defined "
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
int LB_oct_fix(LB *lb, pRegion Region_list, int num_objs) 
{
  int nreg;                                /* number of regions not inserted */
  double pct;                              /* percentage of bad regions */
  int mregions;                            /* total number of regions */
  int i;

  /* initalize variables */
  nreg=0;
  oct_nref=0;
  oct_ncoarse=0;
  
  /* associate the objects to the octants */
  nreg=LB_oct_global_insert_object(lb, Region_list, num_objs);
  LB_oct_global_dref(); 

  if (num_objs)
    pct=100.0*nreg/num_objs;
  else
    pct=0;

  /*
   *  PRINT_IN_ORDER()
   *    printf("refs: %5d  drefs: %5d  bad regions: %5d  = %4.1f%%\n",
   *	   oct_nref,oct_ncoarse,nreg,pct);
   */

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
int LB_oct_global_insert_object(LB *lb, pRegion Region_list, int num_objs) 
{
  pRegion region;                             /* region to be attached */
  int count;                                  /* count of failed insertions */
  int i;

  /* initialize variables */
  count=0;
  region = Region_list;

  /* get the next region to be inserted */
  for(i=0; i<num_objs; i++) {
    if(region != NULL)

#if 0
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
pOctant LB_oct_global_insert(LB *lb, pRegion region) 
{
  pOctant oct;                            /* octree octant */

  oct = NULL;
  /* find the octant which the object lies in */
  oct=LB_oct_global_find(region->Coord);

  /* if octant is found, try to insert the region */
  if (oct)
    if (!LB_oct_subtree_insert(lb, oct, region))         /* inserting region */
      {
	fprintf(stderr,"LB_oct_global_insert: insertion failed\n");
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
  COORD centroid;                              /* coordintes of centroid */
  pRegion entry;

  /* if oct is not terminal, find leaf node the centroid can be attahced to */
  if (!POC_isTerminal(oct)) {
    oct=LB_oct_findOctant(oct,region->Coord);
  }

  if (!oct)
    return(0);

  /* add the region to the octant */
  POC_addRegion(oct, region);

  /* check if octant has too many regions and needs to be refined */
  if(POC_nRegions(oct) > MAXOCTREGIONS)
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
pOctant LB_oct_global_find(COORD point) 
{
  pRList ptr;                      /* ptr used for iterating local root list */
  pOctant root;                    /* root of a subtree */
  pOctant oct;                     /* octree octant */
  COORD min,                      /* minimum bounds coordinates */
        max;                      /* maximum bounds coordinates */

  ptr = POC_localroots();                   /* get a list of all local roots */
  oct = NULL;

  /* iterate through root list to find if point lies inside root's bounds */
  while ((!oct) && (ptr != NULL) ) {
    POC_bounds(ptr->oct,min,max);

    /* ATTN: LB_in_box may need adjusting for the lower bounds */
    if (LB_in_box(point,min,max)) {     /* check if point fits inside */
      /* find the exact octant for point */
      oct=LB_oct_findOctant(ptr->oct,point);
    }
    ptr = ptr->next;
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
pOctant LB_oct_findOctant(pOctant oct, COORD coord) 
{
  pOctant child;                                  /* child of an octant */
  int cnum;                                       /* child number */
  
  /* if octant is terminal, then this is the right octant */
  if (POC_isTerminal(oct))
    return(oct);

  cnum = LB_child_which_wrapper(oct,coord);   /* find closest child to coord */
  child = POC_child(oct, cnum);                            /* get that child */
  /* ATTN: are these local checks necessary? */
  /* if ( !POC_local(child) ) */
  if(!POC_local(oct, cnum))                     /* make sure octant is local */
    return(NULL);

  return(LB_oct_findOctant(child,coord)); /* recursivly search down the tree */
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
void LB_oct_terminal_refine(LB *lb, pOctant oct,int count) 
{
  COORD min,                      /* coordinates of minimum bounds of region */
        max,                      /* coordinates of maximum bounds of region */
        origin,                   /* origin of region */
        centroid;                 /* coordinates of region's centroid */
  pOctant child[8];               /* array of child octants */
  int cnum;                       /* child number */
  int i, j;                       /* index counter */
  pRegion region;                 /* a region to be associated to an octant */
  pRegion entry;
  Region reg;
  COORD cmin[8], cmax[8];
  int new_order;

  for(i=0;i<3;i++)
    min[i] = max[i] = 0;

  /* upper limit of refinement levels */
  /* ATTN: may not be used anymore, but can be put back in if necessary */
  if (count>=10)
    printf("LB_oct_terminal_refine: bailing out at 10 levels\n");
  
  oct_nref++;                                /* increment refinement counter */

  /* octant should be terminal in order to be refined (subdivided) */
  if (!POC_isTerminal(oct)) {
    fprintf(stderr,"ref_octant: oct not terminal\n");
    abort();
  }

  /* get the bounds of an octant */
  POC_bounds(oct,min,max);
  /* calculate the origin from the bounds */
  LB_bounds_to_origin(min,max,origin);

  region = POC_regionlist(oct);             /* Get list while still terminal */
  oct->list = NULL;  /* remove regions from octant, it won't be terminal */

  /* create the children and set their id's */
  j = 0;

  LB_child_bounds_wrapper(oct, cmin, cmax);
  for (i=0; i<8; i++) {
    if(OCT_dimension == 2) {
      if(cmin[i][2] > OCT_gmin[2]) {                /* ignore the z+ octants */
	child[i] = NULL;
	continue;
      }
    }
    
    /*
     * fprintf(stderr,"%d (%lf %lf %lf, %lf %lf %lf)\n", i,
     *	    cmin[0], cmin[1], cmin[2], cmax[0], cmax[1], cmax[2]);
     */

    child[i]=POC_new();                               /* create a new octant */
    POC_setparent(child[i], oct, lb->Proc);   /* set the child->parent link */
    POC_setchildnum(child[i], i);               /* which child of the parent */
    POC_setchild(oct, i, child[i]);            /* set the parent->child link */
#ifdef LGG_MIGOCT
    POC_setID(child[i], LB_oct_nextId());             /* set child id number */
#endif /* LGG_MIGOCT */
    POC_setbounds(child[i], cmin[i], cmax[i]);           /* set child bounds */
    POC_setCpid(oct, i, lb->Proc);      /* set child to be a local octant */
    /*    POC_setOrientation(child[i], 
		       LB_child_orientation(oct->orientation, oct->which));
		       */
  }

  /* assign newly created children to child array*/
  if(OCT_dimension == 3) {
    if(POC_children(oct, child) != 8) {
      /* 
       * if subdivision of oct was successful, oct should have 8 children; 
       * thus a return value of 0 here is a fatal error
       */
      fprintf(stderr, "ref_octant: subdivide failed, %d children.\n",
	      POC_children(oct, child));
      abort();
    }
  }
  else
    if(POC_children(oct, child) != 4) {
      /* 
       * if subdivision of oct was successful, oct should have 4 children; 
       * thus a return value of 0 here is a fatal error
       */
      fprintf(stderr, "ref_octant:subdivide failed, %d children, expected 4\n",
	      POC_children(oct, child));
      abort();
    }

  /* iterate through and find which child each region should belong to */
  while(region != NULL) {
    entry = region->next;
    cnum=LB_child_which_wrapper(oct, region->Coord);
    /* add region to octant's regionlist */
    POC_addRegion(child[cnum], region);
    LB_Free((void **) &region);
    region = entry;
  }

  for (i=0; i<8; i++)                                           /* Recursion */
    if(child[i] != NULL)
      if (POC_nRegions(child[i]) > MAXOCTREGIONS)
	LB_oct_terminal_refine(lb, child[i],count+1);
}

/*****************************************************************************/
/*
 * LB_oct_global_dref()
 * 
 * refine and derefine octree as necessary by number of
 * regions in each octant
 *
 */
void LB_oct_global_dref() 
{
  pRList lroots;                              /* list of all the local roots */

  lroots = POC_localroots();
  while (lroots != NULL) {
    LB_oct_subtree_dref(lroots->oct);
    lroots = lroots->next;
  }
}

/*****************************************************************************/
/*
 * LB_oct_subtree_dref(oct)
 *
 * Coarsen octree so that leaf octants do not have less than MINOCTREGIONS 
 * regions. Refinement takes precedence, so coarsening will not take place 
 * unless all subtrees agree on it.
 */
int LB_oct_subtree_dref(pOctant oct) 
{
  pOctant child;                         /* child of an octant */
  int coarsen;                           /* flag to indicate need to coarsen */
  int i;                                 /* index counter */
  int nregions;                          /* number of regions */
  int total;                             /* total number of regions */

  /* if terminal octant, cannot coarsen on own */
  if (POC_isTerminal(oct)) {
    nregions=POC_nRegions(oct);

    if (nregions > (MAXOCTREGIONS*2) ) {
      printf("LB_oct_subtree_dref: warning: too many (%ld) regions "
	     "in oct %d (id=%d)\n",POC_nRegions(oct),(int)oct,POC_id(oct));
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
    child = POC_child(oct,i);
    /* if child is off processor cannot coarsen */
    /* if (!POC_local(child) ) */
    if(!POC_local(oct, i))
      coarsen=0;
    else {
      /* get the number of region of the child */
      nregions=LB_oct_subtree_dref(child);
      if (nregions<0)
	coarsen=0;
      else
	total+=nregions;                           /* total the region count */
    }
  }

  /* check if octant can be coarsened */
  if (coarsen && total<MAXOCTREGIONS) {
    LB_oct_terminal_coarsen(oct);
    
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
void LB_oct_terminal_coarsen(pOctant oct) 
{
  pOctant child;                         /* child of an octant */
  pRegion region;                        /* region associated with an octant */
  int i;                                 /* index counter */
  void *temp;                            /* temp variable for iterating */
  pRegion regionlist[8];                 /* an array of region lists */

  oct_ncoarse++;                             /* increment coarsening counter */

  for(i=0; i<8; i++) {
    /* get the ith child of an octant */
    child = POC_child(oct,i);
    
    /* cannot coarsen if child is off-processor */
    /* if(!POC_local(child)) X */
    if(!POC_local(oct, i)) {          /* cannot be off-processor */
      fprintf(stderr,"LB_oct_terminal_coarsen: child not local\n");
      abort();
    }
    
    if(!POC_isTerminal(child)) {
      fprintf(stderr,"LB_oct_terminal_coarsen: child not terminal\n");
      abort();
    }
    
    /* get each child's region list */
    regionlist[i] = POC_regionlist(child);
    
    /* delete each child */
    POC_free(child);
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
      POC_addRegion(oct,region);           
      region = region->next;
    }
  }
}

/*****************************************************************************/
void LB_oct_set_maxregions(int max)
{ 
  if (max < 1) {
    fprintf(stderr, "Warning LB_oct_set_maxregions(): %s\n",
	    "illeage input, using default.");
    MAXOCTREGIONS = 1;
  }
  else
    MAXOCTREGIONS = max; 
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
/*
 * LB_oct_roots_in_order(pOctant **roots_ret,int *nroots_ret)
 *
 * Return an array of the local roots, sorted by id
 * Caller must free this array.
 *
 */

void LB_oct_roots_in_order(pOctant **roots_ret, int *nroots_ret)
{
  int nroots;                                 /* number of roots */
  pOctant *roots = NULL;                      /* array of roots */
  Rootid *rootid = NULL;                      /* array of root id's */
  pRList lroots;                              /* list of local roots */
  pRList ptr;                                 /* ptr for iterating */
  pOctant root ;                              /* a root octant */
  int i;                                      /* index counter */
  void *temp;                                 /* temp var used for iterating */

  lroots = POC_localroots();          /* get the list of all the local roots */
  {
    /* find the total number of roots on the list */
    i=0;
    ptr = lroots;
    while(ptr != NULL) {
      i++;
      ptr = ptr->next;
    }
    nroots = i; 

    /* set the return variables */
    *nroots_ret=nroots;
    if(nroots) {
      roots=(pOctant *) LB_Array_Alloc(__FILE__, __LINE__, 1, nroots,
                                       sizeof(pOctant));
      rootid=(Rootid *) LB_Array_Alloc(__FILE__, __LINE__, 1, nroots,
                                       sizeof(Rootid));
      if((roots == NULL) || (rootid == NULL)) {
	fprintf(stderr, "LB_oct_roots_in_order: error in malloc\n");
	abort();
      }
    }
    *roots_ret=roots;
    
    i=0;
    temp=NULL;
    /* place roots in an array to be sorted */
    while (lroots != NULL) {
      rootid[i].ptr = lroots->oct;
      rootid[i].id  = POC_id(lroots->oct);
      i++;
      lroots = lroots->next;
    }
  }
  if (i!=nroots)
    abort();

  /* sort the array of roots */
  qsort(rootid,(size_t)nroots,sizeof(Rootid),(int(*)())idcompare);
  /* qsort(rootid, nroots, sizeof(Rootid), (int(*)())idcompare); */
  
  /* give data to the return variable */
  for (i=0; i<nroots; i++)
    roots[i]=rootid[i].ptr;

  LB_Free((void **) &rootid);
}

/*****************************************************************************/

/*
 * pOctant LB_oct_findId(int id)
 *
 * find the first octant with given id
 *
 */
pOctant LB_oct_findId(int id)
{
  void *temp;
  pOctant oct;

  oct=NULL;
  temp=NULL;
  while (oct=POC_nextDfs(oct)) {
    if (POC_id(oct)==id)
      return(oct);
  }

  return(NULL);
}
#endif /* LGG_MIGOCT */

/*****************************************************************************/

#if 0
/*
 * oct_global_clear()
 *
 * delete all regions from all octants on the local processor
 *
 */
void oct_global_clear()
{
  pPList lroots;                       /* list of all local roots */
  void *temp;                          /* temp variable used for iterations */
  pOctant root;                        /* root octant of a tree */
  pOctant oct;                         /* variable used to reference octants */

  lroots=POC_getRoots();                    /* get a list of all local roots */

  /* 
   * iterate through the list of roots 
   * traverse down the subtree 
   * delete regions associated with octant 
   */
  temp=NULL;
  while (root=PList_next(lroots,&temp))
    for (oct=root; oct; oct=POC_nextDfs(oct,0))
      if (POC_isTerminal(oct))
	POC_delEntities(oct);                          /* delete the regions */

  /* clean up */
  PList_delete(lroots);
}

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
    return(PList_next(OCT_rootlist,temp));


  if (oct=POC_nextDfs(oct, 0))
    return(oct);

  return(PList_next(OCT_rootlist,temp));
}
#endif  /* 0 */
