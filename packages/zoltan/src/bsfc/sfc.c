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
#include "params_const.h"
#include "timer_const.h"
#include <values.h>
#include <limits.h>
#include "hilbert_const.h"
#include "sfc.h"
#include "all_allo_const.h"

/* decent values for this partitioning scheme, the user can
   set them in the as parameters to tune for better performance */
#define BINS_PER_PROC 50 /* minimum amount of coarse bins on each processor */
#define HASHTABLE_DIVIDER 100 /* for this algorithm, need to send out the
				 weights of objects on this processor to 
				 the corresponding bins which may be on
				 other processors.  a hashtable is used so
				 that we do not have to allocate an array
				 of size BINS_PER_PROC * Num_Proc.  only 
				 an array of size 
				 BINS_PER_PROC * Num_Proc / HASHTABLE_DIVIDER
				 needs to be allocated which allows the 
				 algorithm to scale */
#define MAX_CUTS_IN_BIN 10 /* maximum amount of cuts in a coarse level bin */
#define SUBBINS_PER_BIN 20 /* amount of subbins a bin is divided into */
#define MAX_REFINEMENT_LEVEL 20 /* amount of refinement of the bins */
#define BSFC_BOUNDING_BOX_EPSILON 0.0000001 /* used to increase the bounding
					      box slightly so that no objects
					      are on the boundary of the box */

int Zoltan_BSFC_create_refinement_info(ZZ *zz, int* number_of_cuts, 
			       float* global_actual_work_allocated,
			       int wgt_dim, float* total_weight_array, 
			       float* work_percent_array, int num_vert_in_cut, 
			       BSFC_VERTEX_PTR vert_in_cut_ptr, 
			       float* wgts_in_cut_ptr,
			       float** work_prev_allocated_ptr);

int Zoltan_BSFC_create_bins(ZZ *zz, int num_local_objects, int wgt_dim,
		    BSFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[],
		    int* amount_of_bits_used, int size_of_unsigned, 
		    float* global_actual_work_allocated, 
		    float *work_percent_array, float* total_weight_array,
		    int* balanced_flag, BSFC_VERTEX_PTR *vert_in_cut_ptr,
		    float** wgts_in_cut_ptr, int* num_vert_in_cut,
		    int* number_of_cuts, int bins_per_proc, 
		    int hashtable_divider, ZOLTAN_COMM_OBJ **plan, 
		    int* num_vert_sent, int max_cuts_in_bin);

static PARAM_VARS BSFC_params[] = {
  { "BSFC_BINS_PER_PROC", NULL, "INT" },
  { "BSFC_HASHTABLE_DIVIDER", NULL, "INT" }, 
  { "BSFC_MAX_CUTS_IN_BIN", NULL, "INT" },
  { "BSFC_SUBBINS_PER_BIN", NULL, "INT" },
  { "BSFC_MAX_REFINEMENT_LEVEL", NULL, "INT" },
  { NULL, NULL, NULL } };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int Zoltan_BSFC_Set_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, BSFC_params, &result, &index);

    return(status);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/* Space filling curve (BSFC) partioning routine */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int Zoltan_BSFC(
  ZZ *zz,                       /* The Zoltan structure with info for
                                   the RCB balancer.                         */
  int *num_import,              /* Not computed.  Set to -1. */
  ZOLTAN_ID_PTR *import_global_ids, /* Not computed. */
  ZOLTAN_ID_PTR *import_local_ids,  /* Not computed. */
  int **import_procs,           /* Not computed. */
  int *num_export,              /* Number of local objects assigned to another
                                   processor in the new decomposition. */
  ZOLTAN_ID_PTR *export_global_ids, /* Returned value:  array of global IDs of
                                   local objects assigned to other processors
				   for the new decomposition.                */
  ZOLTAN_ID_PTR *export_local_ids,  /* Returned value:  array of local IDs of
                                   local objects assigned to other processors
				   for the new decomposition.                */ 
  int **export_procs            /* Returned value:  array of processor IDs for
                                   the local objects to be sent to in
                                   this processor's new decomposition.       */
)
{
  char    yo[] = "Zoltan_BSFC";
  int wgt_dim = zz->Obj_Weight_Dim;   /* dimension of weights of each object */
  int num_dims;                       /* geometric dimension */
  int ierr, i, j;                     /* local variables */
  double bounding_box[6], global_bounding_box[6]; /* local and global geometric
						     bounding boxes for the objects */
  int num_local_objects;              /* the number of objects this processor owns */
  BSFC_VERTEX_PTR sfc_vert_ptr;        /* array that stores the sfc objects */
  int num_gid_entries = zz->Num_GID;  
  int num_lid_entries = zz->Num_LID;
  ZOLTAN_ID_PTR global_ids = NULL;
  ZOLTAN_ID_PTR local_ids = NULL;
  float* objs_wgt = NULL;             /* array of objects weights */
  int size_of_unsigned;               /* minimum size of an unsigned integer,
					 used only for heterogeneous systems */
  float *global_actual_work_allocated = NULL; /* cumulative actual work allocated */
  float *total_weight_array = NULL;   /* sum of weights (length of wgt_dim) */
  float *work_percent_array = NULL;   /* the cumulative percent of work each
					 processor should ideally get */
  int balanced_flag;                  /* flag to indicate if all partitions
					 are balanced */
  BSFC_VERTEX_PTR vert_in_cut_ptr = NULL; /* array of sfc objects that are in a 
					    cut bin and get partitioned with 
					    the multi-level scheme */
  float* wgts_in_cut_ptr = NULL;      /* array of weights for sfc objects in 
					 a cut bin */
  int num_vert_in_cut;                /* number of sfc objects in the cut bin */
  int number_of_cuts = 0; /* maximum amount of cuts in a coarse bin on this processor */
  int amount_of_bits_used = 0;        /* amount of bits used in calculating the
					 bin an sfc object belongs to */
  ZOLTAN_COMM_OBJ *plan;                     /* used to put all sfc objects that were 
					 moved to a new proc back on the original proc */
  int comm_tag = 8765;                /* randomly chosen communication tag */
  int num_vert_sent;                  /* the number of sfc objects that this processor 
					 sent to other processors */
  int refinement_level_counter = 0;   /* counter to keep track of how many 
					 levels of bin refinement have been performed */
  int max_cuts_in_bin, bin_refinement_method, max_refinement_level,
    subbins_per_bin, hashtable_divider, bins_per_proc; /* tuning parameters */
  double* coords; /* array for objects coordinates */

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* set the of parameters */
  Zoltan_Bind_Param(BSFC_params,"BSFC_BINS_PER_PROC",(void*) &bins_per_proc);
  Zoltan_Bind_Param(BSFC_params,"BSFC_HASHTABLE_DIVIDER",(void*) &hashtable_divider); 
  Zoltan_Bind_Param(BSFC_params,"BSFC_MAX_CUTS_IN_BIN",(void*) &max_cuts_in_bin);
  Zoltan_Bind_Param(BSFC_params,"BSFC_SUBBINS_PER_BIN",(void*) &subbins_per_bin); 
  Zoltan_Bind_Param(BSFC_params,"BSFC_MAX_REFINEMENT_LEVEL",(void*) &max_refinement_level);
  bins_per_proc = BINS_PER_PROC;
  hashtable_divider = HASHTABLE_DIVIDER; 
  max_cuts_in_bin = MAX_CUTS_IN_BIN;
  subbins_per_bin = SUBBINS_PER_BIN;
  max_refinement_level = MAX_REFINEMENT_LEVEL;

  Zoltan_Assign_Param_Vals(zz->Params, BSFC_params, zz->Debug_Level, zz->Proc,
		       zz->Debug_Proc);

  /* make sure that all parameters have feasible values */
  if(bins_per_proc <= 0) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, 
		  "BSFC_BINS_PER_PROC parameter must be greater than 0.");
    bins_per_proc = BINS_PER_PROC;
  }
  if(hashtable_divider <= 0) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, 
		  "BSFC_HASH_TABLE_DIVIDER parameter must be greater than 0.");
    hashtable_divider = HASHTABLE_DIVIDER;
  }
  if(max_cuts_in_bin <= 0) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, 
		  "BSFC_MAX_CUTS_IN_BIN parameter must be greater than 0.");
    max_cuts_in_bin = MAX_CUTS_IN_BIN;
  }
  if(subbins_per_bin <= 1) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, 
		  "BSFC_SUBBINS_PER_BIN parameter must be greater than 1.");
    subbins_per_bin = BINS_PER_PROC;
  }

  /* Initializations in case of early exit. */
  *num_import = -1;  /* We don't compute the import map. */
  *num_export = -1;  

  /* get the dimension of the problem and make sure that it is either 2 or 3 */
  num_dims = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &ierr);
  if(ierr != 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned from user function Get_Num_Geom.");
    return(ierr);
  }

  if(num_dims != 2 && num_dims != 3) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Incompatible space dimension for BSFC. Space dimension must be 2 or 3.");
    return(ZOLTAN_FATAL);
  }

  /*for heterogeneous systems where the size of an unsigned integer may be different,
    find the minimum size bytes used to store an unsigned integer*/
  i = sizeof( unsigned );
  ierr = MPI_Allreduce(&i, &size_of_unsigned, 1, MPI_INT, 
		       MPI_MIN, zz->Communicator);

  /* get application data (number of objects, ids, weights, and coords */
  num_local_objects = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned from user function Get_Num_Obj.");
    return(ierr);
  }

  if (num_local_objects > 0) {
    global_ids = ZOLTAN_MALLOC_GID_ARRAY(zz, num_local_objects);
    local_ids  = ZOLTAN_MALLOC_LID_ARRAY(zz, num_local_objects);

    if (!(global_ids) || (zz->Num_LID && !(local_ids))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
  }
  
  /*
   *  Get list of objects' IDs and weights.
   */
  if (num_local_objects > 0) {
    if (wgt_dim) {
      
      /* 
       *  Allocate space for object weights.
       */
      
      objs_wgt    = 
	(float *) ZOLTAN_MALLOC(wgt_dim*(num_local_objects)*sizeof(float));
      if (!objs_wgt) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return(ZOLTAN_MEMERR);
      }
      for (i = 0; i < wgt_dim*num_local_objects; i++) objs_wgt[i] = 0.;
    }  
  }
  Zoltan_Get_Obj_List(zz, global_ids, local_ids, wgt_dim, objs_wgt, &ierr);

  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
		   "Error returned from user function Zoltan_Get_Obj_List.");
    return(ierr);
  }
  /* if object weights are not defined, all object weights are 1.0 */
  if (wgt_dim < 1) {
    if (num_local_objects > 0) {
      objs_wgt    = (float *) ZOLTAN_MALLOC((num_local_objects)*sizeof(float));
      if (!objs_wgt) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return(ZOLTAN_MEMERR);
      }
      for (i = 0; i < num_local_objects; i++) objs_wgt[i] = 1.;
    }
    wgt_dim = 1;  /* set object weight dimension to 1 */
  }  
  
  sfc_vert_ptr = (BSFC_VERTEX_PTR) ZOLTAN_MALLOC(num_local_objects * sizeof(BSFC_VERTEX));
  if(num_local_objects != 0 && sfc_vert_ptr == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
  }
  coords = (double*) ZOLTAN_MALLOC(sizeof(double) * num_local_objects * num_dims);
  if(num_local_objects != 0 && coords == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  /* get the geometric coordinates of the objects */
  for(i=0;i<num_local_objects;i++) {
    zz->Get_Geom(zz->Get_Geom_Data, num_gid_entries, num_lid_entries,
		 &(global_ids[i*num_gid_entries]), &(local_ids[i*num_lid_entries]),
		 (coords+i*num_dims), &ierr);
    
    if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                     "Error returned from user defined Get_Geom function.");
      return(ierr);
    }  
  }

  /* go through and find bounding box for entire domain */  
  for(i=0;i<num_dims+num_dims;i++)
    bounding_box[i] = MAXDOUBLE;
  
  for(i=0;i<num_local_objects;i++) 
    for(j=0;j<num_dims;j++)  {
      if(coords[i*num_dims+j] < bounding_box[j])
	bounding_box[j] = coords[i*num_dims+j];
      if(coords[i*num_dims+j] > -bounding_box[j+num_dims])
	bounding_box[j+num_dims] = -coords[i*num_dims+j];
    }
  
  ierr = MPI_Allreduce(bounding_box, global_bounding_box, 2*num_dims,
		       MPI_DOUBLE, MPI_MIN, zz->Communicator);
  for(i=num_dims;i<num_dims+num_dims;i++)
    global_bounding_box[i] = - global_bounding_box[i];

  /* enlarge global_bounding_box slightly */
  for(i=0;i<num_dims;i++) {
    double delta = global_bounding_box[i+num_dims] - global_bounding_box[i];
    if(delta > 0 ) {
      global_bounding_box[i] = 
	global_bounding_box[i]*(1. - BSFC_BOUNDING_BOX_EPSILON);
      global_bounding_box[i+num_dims] =
	global_bounding_box[i+num_dims]*(1. + BSFC_BOUNDING_BOX_EPSILON);
    }
    else {
      global_bounding_box[i] = global_bounding_box[i] - BSFC_BOUNDING_BOX_EPSILON;
      global_bounding_box[i+num_dims] =
	global_bounding_box[i+num_dims] + BSFC_BOUNDING_BOX_EPSILON;
    }
  }
  /* done creating global bounding box */

  /* Normalize space coordinates and fill in sfc_vertex info */
  Zoltan_BSFC_create_info(zz, global_bounding_box, (global_bounding_box+num_dims), 
		  num_dims, num_local_objects, wgt_dim, sfc_vert_ptr,
		  coords);

  ZOLTAN_FREE(&coords);

  global_actual_work_allocated=(float*) ZOLTAN_MALLOC(sizeof(float)*wgt_dim* zz->Num_Proc);
  if(!global_actual_work_allocated) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  work_percent_array = (float*) ZOLTAN_MALLOC(sizeof(float) * zz->Num_Proc);
  if(!work_percent_array) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  total_weight_array = (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
  if(!total_weight_array) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  /*create bins, fill global weight vector and perform initial partition of the bins*/
  ierr = Zoltan_BSFC_create_bins(zz, num_local_objects, wgt_dim, sfc_vert_ptr, objs_wgt,
			 &amount_of_bits_used, size_of_unsigned, 
			 global_actual_work_allocated, work_percent_array, 
			 total_weight_array, &balanced_flag, &vert_in_cut_ptr,
			 &wgts_in_cut_ptr, &num_vert_in_cut, &number_of_cuts,
			 bins_per_proc, hashtable_divider, &plan,
			 &num_vert_sent, max_cuts_in_bin); 
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_BSFC_create_bins function.");
      return(ierr);
  }
 
  if(balanced_flag != BSFC_BALANCED) { 
    int* local_balanced_flag_array; /* used to indicate which partitions on this 
				       processor are already balanced - useful
				       for when more than one cut in a bin */
    int max_cuts_in_bin;
    int* ll_bins_head; /* used to indicate the beginning of the linklist */
    float* work_prev_allocated = NULL; /* stores the weights of all 
					  objects before a cut */

    if(num_vert_in_cut == 0 || zz->Proc == 0) 
      balanced_flag = BSFC_BALANCED;
    ierr = Zoltan_BSFC_create_refinement_info(zz, &number_of_cuts, 
				      global_actual_work_allocated, wgt_dim,
				      total_weight_array, work_percent_array,
				      num_vert_in_cut, vert_in_cut_ptr, 
				      wgts_in_cut_ptr, &work_prev_allocated);
    if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in create_refinement_info function.");
      return(ierr);
    }    

    ll_bins_head = (int*) ZOLTAN_MALLOC(sizeof(int) * (1+number_of_cuts));
    if(ll_bins_head == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory."); 
      return(ZOLTAN_MEMERR);
    }
    max_cuts_in_bin= number_of_cuts;
    if(number_of_cuts == 0)
      balanced_flag = BSFC_BALANCED;

    if(ll_bins_head != NULL)
      ll_bins_head[number_of_cuts] = 0;  /* the first link list starts off
					    at array location number_of_cuts-1 ! */
    for(i=0;i<number_of_cuts;i++)
      ll_bins_head[i] = -1;
   
    local_balanced_flag_array = (int*) ZOLTAN_MALLOC(sizeof(int) * (1+number_of_cuts));
    if(local_balanced_flag_array == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory."); 
      return(ZOLTAN_MEMERR);
    }
    for(i=0;i<number_of_cuts;i++)
      local_balanced_flag_array[i] = BSFC_BALANCED;
    local_balanced_flag_array[number_of_cuts] = balanced_flag;

    /* refine bins until a satisfactory partition tolerance is attained */
    while(balanced_flag != BSFC_BALANCED &&
	  refinement_level_counter < max_refinement_level) { 
      ierr = Zoltan_BSFC_refine_partition(zz, &balanced_flag, 
					  &amount_of_bits_used, num_vert_in_cut,
					  vert_in_cut_ptr, size_of_unsigned, 
					  wgt_dim, wgts_in_cut_ptr, 
					  work_percent_array, total_weight_array,
					  global_actual_work_allocated, 
					  number_of_cuts, &max_cuts_in_bin,
					  ll_bins_head, work_prev_allocated, 
					  subbins_per_bin, local_balanced_flag_array,
					  bin_refinement_method);
      if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_BSFC_refine_partition_level function.");
	return(ierr);
      }
      refinement_level_counter++;
    }
    ZOLTAN_FREE(&local_balanced_flag_array);
    ZOLTAN_FREE(&work_prev_allocated);
    ZOLTAN_FREE(&ll_bins_head);
  }
  
  /* if the objects were moved to different processors,
     we need to move them back now */
  if(plan != NULL) {
    BSFC_VERTEX_PTR recv_buf = NULL;
    int counter = 0;
    if (num_vert_sent > 0) {
      recv_buf = (BSFC_VERTEX_PTR) ZOLTAN_MALLOC(sizeof(BSFC_VERTEX)*num_vert_sent);
      if(recv_buf==NULL) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return(ZOLTAN_MEMERR);
      }
    }
    ierr = Zoltan_Comm_Do_Reverse(plan, comm_tag, (char*) vert_in_cut_ptr,
			      sizeof(BSFC_VERTEX), NULL, (char*) recv_buf);
    if(ierr == ZOLTAN_WARN) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Do_Reverse.");
    }
    else if(ierr == ZOLTAN_FATAL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Do_Reverse.");
      return(ierr);
    }      
    else if(ierr == ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Do_Reverse.");
      return(ierr);
    }
    /* put objects back in sfc_vert_ptr array the way they 
       were copied from this array */
    for(i=0;i<num_local_objects;i++)
      if(sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT)  {
	sfc_vert_ptr[i] = recv_buf[counter];
	counter++;
      }      

    ZOLTAN_FREE(&recv_buf);
    ierr = Zoltan_Comm_Destroy(&plan);
    if(ierr == ZOLTAN_WARN) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Destroy.");
    }
    else if(ierr == ZOLTAN_FATAL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Destroy.");
      return(ierr);
    }      
    else if(ierr == ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Destroy.");
      return(ierr);
    }
  }

  ZOLTAN_FREE(&vert_in_cut_ptr);
  ZOLTAN_FREE(&wgts_in_cut_ptr);
  ZOLTAN_FREE(&global_actual_work_allocated);
  ZOLTAN_FREE(&work_percent_array);  
  ZOLTAN_FREE(&total_weight_array);

  /* add up stuff to export */
  *num_export = 0;
  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].destination_proc != zz->Proc)
      (*num_export)++;

  ierr = Zoltan_Special_Malloc(zz, (void**) export_global_ids,
			   *num_export, ZOLTAN_SPECIAL_MALLOC_GID);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  ierr = Zoltan_Special_Malloc(zz, (void**) export_local_ids,
			   *num_export, ZOLTAN_SPECIAL_MALLOC_LID);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  ierr = Zoltan_Special_Malloc(zz, (void**) export_procs,
			   *num_export, ZOLTAN_SPECIAL_MALLOC_INT);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  
  /* fill in the export data */
  j = 0;
  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].destination_proc != zz->Proc) {
      *((*export_procs)+j) = sfc_vert_ptr[i].destination_proc;
      ZOLTAN_SET_GID(zz, (*export_global_ids+j*num_gid_entries), (global_ids+i*num_gid_entries));
      ZOLTAN_SET_LID(zz, (*export_local_ids+j*num_lid_entries), (local_ids+i*num_lid_entries));
      j++;
    }

  ZOLTAN_FREE(&objs_wgt);
  
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&sfc_vert_ptr);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}

/* create info before starting the multi-level refinement of the bins 
   NOTE:  this routine only works for objects with one weight!
*/

int Zoltan_BSFC_create_refinement_info(ZZ *zz, int* number_of_cuts, 
			       float* global_actual_work_allocated,
			       int wgt_dim, float* total_weight_array, 
			       float* work_percent_array,
			       int num_vert_in_cut, 
			       BSFC_VERTEX_PTR vert_in_cut_ptr,
			       float* wgts_in_cut_ptr, 
			       float** work_prev_allocated_ptr)
{
  char    yo[] = "create_refinement_info";
  float my_work, *work_array;
  int i = 0, j;

  /* find out how many cuts are in this bin. */
  my_work = global_actual_work_allocated[(zz->Proc)*wgt_dim];
  while(my_work > total_weight_array[0] * work_percent_array[zz->Proc-i])
    i++;
  *number_of_cuts = i;
  if(num_vert_in_cut == 0 || zz->Proc ==0) {
    *number_of_cuts = 0;
    return ZOLTAN_OK;
  }

  /* create link list for objects in the array.  link list
     is set up so that the objects in the array are
     traversed consecutively */
  for(i=0;i<(num_vert_in_cut-1);i++)
    vert_in_cut_ptr[i].next_sfc_vert_index = i+1;
  vert_in_cut_ptr[num_vert_in_cut-1].next_sfc_vert_index = -1;

  /* update work previously allocated to include work in all
     bins with higher keys than this bin */
  work_array = (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
  if(work_array == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  
  for(i=0;i<wgt_dim;i++)
    work_array[i] = 0;
  for(i=0;i<num_vert_in_cut;i++) 
    for(j=0;j<wgt_dim;j++)
      work_array[j] += wgts_in_cut_ptr[i*wgt_dim+j];
  
  *work_prev_allocated_ptr = 
    (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim * (*number_of_cuts+1));
  if(*work_prev_allocated_ptr == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
      
  /* update work previously allocated to include work in all bins
     with higher keys than this bin */
  for(i=0;i<wgt_dim;i++)
    *((*work_prev_allocated_ptr)+wgt_dim*(*number_of_cuts)+i) =
      global_actual_work_allocated[(zz->Proc)*wgt_dim+i] - work_array[i];

  ZOLTAN_FREE(&work_array);
  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
