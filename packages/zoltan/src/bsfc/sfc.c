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
#include "params_const.h"
#include "timer_const.h"
#include "comm_const.h"
#include "lb_util_const.h"
#include <values.h>
#include <limits.h>
#include "hilbert_const.h"
#include "sfc_const.h"
#include "sfc.h"
#include "all_allo_const.h"

/* okay values for this partitioning scheme, the user can set them in the as parameters also */
#define BINS_PER_PROC 10 
#define HASHTABLE_DIVIDER 20
#define MAX_CUTS_PER_BIN 10 /* maximum amount of cuts in a coarse level bin */

int sfc_create_bins(LB* lb, int num_local_objects, 
		    int wgt_dim, SFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[], int* amount_of_bits_used,
		    int sfc_keylength, int size_of_unsigned, unsigned imax, 
		    float* global_actual_work_allocated, float *work_percent_array, 
		    float* total_weight_array, int* balanced_flag,
		    SFC_VERTEX_PTR *vert_in_cut_ptr, float** wgts_in_cut_ptr, 
		    int* num_vert_in_cut, int* number_of_cuts, int bins_per_proc, 
		    int hashtable_divider, COMM_OBJ **plan, int* num_vert_sent);

static PARAM_VARS SFC_params[] = {
  { "SFC_BINS_PER_PROC", NULL, "INT" },
  { "SFC_HASHTABLE_DIVIDER", NULL, "INT" },
  { NULL, NULL, NULL } };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int LB_Set_SFC_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */

    status = LB_Check_Param(name, val, SFC_params, &result, &index);

    return(status);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/* Space filling curve (SFC) partioning routine */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int LB_sfc(
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
  char    yo[] = "LB_sfc";
  int wgt_dim = lb->Obj_Weight_Dim;              /* dimension of weights of each object */
  int num_dims;
  int ierr, i, j;
  double bounding_box[6];
  double global_bounding_box[6];
  int num_local_objects;
  double delta;
  SFC_VERTEX_PTR sfc_vert_ptr;
  int num_gid_entries = lb->Num_GID;
  int num_lid_entries = lb->Num_LID;
  int max_obj;
  LB_ID_PTR global_ids;
  LB_ID_PTR local_ids;
  float* objs_wgt;
  int size_of_unsigned;
  unsigned imax;
  float *global_actual_work_allocated = NULL, *total_weight_array = NULL;
  float *work_percent_array;
  int balanced_flag;
  SFC_VERTEX_PTR vert_in_cut_ptr = NULL;
  float* wgts_in_cut_ptr = NULL;
  int num_vert_in_cut;
  int number_of_cuts = 0; /* amount of cuts in the coarse bin that gets a cut */
  int bins_per_proc = 10;  /* bring this in as a parameter later on... */
  int amount_of_bits_used = 0;
  int bin_level = SFC_FIRST_LEVEL_FLAG;
  int hashtable_divider;
  COMM_OBJ *plan;
  int comm_tag = 8765; /* randomly chosen communication tag */
  int num_vert_sent;

  
  printf("in sfc partitioning\n");

  /* set up a couple of parameters */
  LB_Bind_Param(SFC_params, "SFC_BINS_PER_PROC", (void*) &bins_per_proc);
  LB_Bind_Param(SFC_params, "SFC_HASHTABLE_DIVIDER", (void*) &hashtable_divider); 
  bins_per_proc = BINS_PER_PROC;
  hashtable_divider = HASHTABLE_DIVIDER; 

  LB_Assign_Param_Vals(lb->Params, SFC_params, lb->Debug_Level, lb->Proc,
		       lb->Debug_Proc);

  /* Initializations in case of early exit. */
  *num_import = -1;  /* We don't compute the import map. */
  *num_export = -1;  

  /* get the dimension of the problem and make sure that it is either 2 or 3 */
  num_dims = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
  if(ierr != 0) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Geom.");
    return(ierr);
  }

  if(num_dims != 2 && num_dims != 3) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Incompatible space dimension for SFC. Space dimension must be 2 or 3.");
    ierr = LB_FATAL;
    return(ierr);
  }

  /*for heterogeneous systems where the size of an unsigned integer may be different,
    find the minimum size bytes used to store an unsigned integer*/
  i = sizeof( unsigned );
  ierr = MPI_Allreduce(&i, &size_of_unsigned, 1, MPI_INT, 
		       MPI_MIN, lb->Communicator);
  if(size_of_unsigned == sizeof( unsigned ))
    imax = IScale ;
  else 
    imax = pow(2, size_of_unsigned*8) - 1;
  
  /* get application data (number of objects, ids, weights, and coords */
  
  num_local_objects = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Obj.");
    return(ierr);
  }
  /*max_objects might need to be more than num_local_objects (max_objects taken from rcg/shared.c)*/
  max_obj = num_local_objects;
  global_ids = LB_MALLOC_GID_ARRAY(lb, (max_obj));
  local_ids  = LB_MALLOC_LID_ARRAY(lb, (max_obj));

  
  if (!(global_ids) || (lb->Num_LID && !(local_ids))) {
    LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    return(LB_MEMERR);
  }
  
  /*
   *  Get list of objects' IDs and weights.
   */
  if (num_local_objects > 0) {
    
    if (wgt_dim) {
      
      /* 
       *  Allocate space for object weights.
       */
      
      objs_wgt    = (float *) LB_MALLOC(wgt_dim*(num_local_objects)*sizeof(float));
      if (!objs_wgt) {
        LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        return(LB_MEMERR);
      }
      for (i = 0; i < wgt_dim*num_local_objects; i++) objs_wgt[i] = 0.;
    }  
  }
  LB_Get_Obj_List(lb, global_ids, local_ids, wgt_dim, objs_wgt, &ierr);

  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, 
		   "Error returned from user function LB_Get_Obj_List.");
    return(ierr);
  }
  /* if object weights are not defined, all object weights are 1.0 */
  if (num_local_objects > 0 && wgt_dim < 1) {
    objs_wgt    = (float *) LB_MALLOC((num_local_objects)*sizeof(float));
    if (!objs_wgt) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
    }
    for (i = 0; i < num_local_objects; i++) objs_wgt[i] = 1.;
    wgt_dim = 1;  /* set object weight dimension to 1 */
  }  
    
  
  sfc_vert_ptr = (SFC_VERTEX_PTR) LB_MALLOC(num_local_objects * sizeof(SFC_VERTEX));
  if(num_local_objects != 0 && sfc_vert_ptr == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
  }
  for(i=0;i<num_local_objects;i++) {
    lb->Get_Geom(lb->Get_Geom_Data, num_gid_entries, num_lid_entries,
		 &(global_ids[i*num_gid_entries]), &(local_ids[i*num_lid_entries]),
		 sfc_vert_ptr[i].coord, &ierr);
    
    if (ierr == LB_FATAL || ierr == LB_MEMERR) {
      LB_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user defined Get_Geom function.");
      return(ierr);
    }    

    
  }
  /* go through and find bounding box for entire domain */  
  for(i=0;i<num_dims+num_dims;i++)
    bounding_box[i] = MAXDOUBLE;
  
  for(i=0;i<num_local_objects;i++) 
    for(j=0;j<num_dims;j++)  {
      if(sfc_vert_ptr[i].coord[j] < bounding_box[j])
	bounding_box[j] = sfc_vert_ptr[i].coord[j];
      if(sfc_vert_ptr[i].coord[j] > -bounding_box[j+num_dims])
	bounding_box[j+num_dims] = -sfc_vert_ptr[i].coord[j];
    }
  
  ierr = MPI_Allreduce(bounding_box, global_bounding_box, 2*num_dims, MPI_DOUBLE, 
		       MPI_MIN, lb->Communicator);
  for(i=num_dims;i<num_dims+num_dims;i++)
    global_bounding_box[i] = - global_bounding_box[i];

  /* enlarge global_bounding_box slightly */
  for(i=0;i<num_dims;i++) {
    delta = global_bounding_box[i+num_dims] - global_bounding_box[i];
    if(delta > 0 ) {
      global_bounding_box[i] = global_bounding_box[i]*(1. - SFC_BOUNDING_BOX_EPSILON);
      global_bounding_box[i+num_dims] = global_bounding_box[i+num_dims]*(1. + SFC_BOUNDING_BOX_EPSILON);
    }
    else {
      global_bounding_box[i] = global_bounding_box[i] - SFC_BOUNDING_BOX_EPSILON;
      global_bounding_box[i+num_dims] = global_bounding_box[i+num_dims] + SFC_BOUNDING_BOX_EPSILON;
    }
  }
  /*done creating global bounding box*/
  /* Normalize space coordinates and fill in sfc_vertex info */
  sfc_create_info(lb, global_bounding_box, (global_bounding_box+num_dims), 
		  num_dims, num_local_objects, wgt_dim, sfc_vert_ptr, SFC_KEYLENGTH);

  global_actual_work_allocated = (float*) LB_MALLOC(sizeof(float) * wgt_dim * lb->Num_Proc);
  if(!global_actual_work_allocated) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
  }
  work_percent_array = (float*) LB_MALLOC(sizeof(float) * lb->Num_Proc);
  if(!work_percent_array) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
  }
  total_weight_array = (float*) LB_MALLOC(sizeof(float) * wgt_dim);
  if(!total_weight_array) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
  }
  /*create bins, fill global weight vector and perform initial partition of the bins*/
  ierr = sfc_create_bins(lb, num_local_objects, wgt_dim, sfc_vert_ptr, objs_wgt, &amount_of_bits_used, 
			 SFC_KEYLENGTH, size_of_unsigned, imax, global_actual_work_allocated, 
			 work_percent_array, total_weight_array, &balanced_flag, &vert_in_cut_ptr,
			 &wgts_in_cut_ptr, &num_vert_in_cut, &number_of_cuts, bins_per_proc,
			 hashtable_divider, &plan, &num_vert_sent); 
  if(ierr != LB_OK && ierr != LB_WARN) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error in sfc_create_bins function.");
      return(ierr);
  }

  /* for debugging */
/*  if(lb->Proc ==0) {
    printf("input balanced %d or unbalanced %d\n",SFC_BALANCED, SFC_NOT_BALANCED);
      scanf("%d",&balanced_flag);
      printf("\n");    
  }
  i = MPI_Bcast(&balanced_flag, 1, MPI_INT, 0, lb->Communicator);*/
  /* done debugging */
  printf("this routine will only do the coarse bin partitioning\n)");
  balanced_flag = SFC_BALANCED;
  /* refine bins until a satisfactory partition tolerance is attained */
  while(balanced_flag != SFC_BALANCED) {  
    int ll_bins_head[MAX_CUTS_PER_BIN]; /* used to indicate the beginning of the linklist
					   after refinement */
    printf("****doing refinement of the bins****\n");

    ierr = sfc_refine_partition_level(lb, &balanced_flag, &amount_of_bits_used,
				      &bin_level, num_vert_in_cut, vert_in_cut_ptr,
				      SFC_KEYLENGTH, size_of_unsigned, imax, wgt_dim,
				      wgts_in_cut_ptr, work_percent_array,
				      total_weight_array, global_actual_work_allocated, 
				      &number_of_cuts, ll_bins_head);
    if(ierr != LB_OK && ierr != LB_WARN) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error in sfc_refine_partition_level function.");
      return(ierr);
    }
  }
  /* if the objects were moved to different processors, we need to move them back now */
  if(plan != NULL) {
    SFC_VERTEX_PTR recv_buf;
    int counter = 0;
    recv_buf = (SFC_VERTEX_PTR) LB_MALLOC(sizeof(SFC_VERTEX) * num_vert_sent);
    if(recv_buf==NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
    }
    ierr = LB_Comm_Do_Reverse(plan, comm_tag, (char*) vert_in_cut_ptr,
			      sizeof(SFC_VERTEX), NULL, (char*) recv_buf);
    /* put objects back in sfc_vert_ptr array the way they were copied from this array */
    for(i=0;i<num_local_objects;i++)
      if(sfc_vert_ptr[i].cut_bin_flag == SFC_CUT)  {
	sfc_vert_ptr[i] = recv_buf[counter];
	counter++;
      }      

    LB_FREE(&recv_buf);
    ierr = LB_Comm_Destroy(&plan);
  }

  LB_FREE(&vert_in_cut_ptr);
  LB_FREE(&wgts_in_cut_ptr);
  LB_FREE(&global_actual_work_allocated);
  LB_FREE(&work_percent_array);  
  LB_FREE(&total_weight_array);
  LB_FREE(&objs_wgt);


  /* add up stuff to export */
  *num_export = 0;
  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].destination_proc != lb->Proc)
      (*num_export)++;
  ierr = LB_Special_Malloc(lb, (void**) export_global_ids, *num_export, LB_SPECIAL_MALLOC_GID);
  ierr = LB_Special_Malloc(lb, (void**) export_local_ids, *num_export, LB_SPECIAL_MALLOC_LID);
  ierr = LB_Special_Malloc(lb, (void**) export_procs, *num_export, LB_SPECIAL_MALLOC_INT);

  j = 0;
  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].destination_proc != lb->Proc) {
      *((*export_procs)+j) = sfc_vert_ptr[i].destination_proc;
      LB_SET_GID(lb, (*export_global_ids+j), (global_ids+i));
      LB_SET_LID(lb, (*export_local_ids+j), (local_ids+i));
      j++;
    }
  
  LB_FREE(&global_ids);
  LB_FREE(&local_ids);
  LB_FREE(&sfc_vert_ptr);

  MPI_Barrier(lb->Communicator);
  printf("proc %d is leaving sfc balancing routines after exporting %d objects\n",lb->Proc, *num_export);

  return LB_OK;
}



