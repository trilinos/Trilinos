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


#include "zz_const.h"
#include "zz_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines to set the partition sizes.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Set_Part_Sizes(ZZ *zz, int local_global,
    int len, int *wgt_idx, int *part_ids, float *part_sizes)
{
/*
 *  Function to set the desired partition sizes:
 *  Input:
 *    zz            --  The Zoltan structure to which this method
 *                      applies.
 *    local_global  --  Local or global partition numbers? (only global for now)
 *    len           --  Length of arrays wgt_idx, part_idx, part_sizes
 *    wgt_idx       --  Array of indices between 0 and Obj_Wgt_Dim-1
 *    part_ids      --  Array of partition ids (local or global)
 *    part_sizes    --  Array of floats that gives the desired partition 
 *                      size for each weight and each partition, i.e., 
 *                      part_sizes[i] corresponds to wgt_idx[i] and part_id[i]
 *
 *  Output:
 *    zz->LB.*      --  Appropriate fields set to designated values.
 */

  char *yo = "Zoltan_LB_Set_Part_Sizes";
  int i, j;
  float *temp;
  char msg[128];
  int error = ZOLTAN_OK;
  int part_dim = zz->Obj_Weight_Dim;        /* Current value of obj_wgt_dim. */

  if (part_dim==0) part_dim = 1;

  /* Check if we can use existing part_sizes array */
  if (zz->LB.Part_Sizes == NULL){

    /* Store max values, which give the dimensions of Part_Sizes. */
    zz->LB.Max_Part_Dim = part_dim;
    zz->LB.Max_Global_Parts = zz->LB.Num_Global_Parts;

    /* Allocate fresh Part_Sizes array . */
    zz->LB.Part_Sizes = (float *) ZOLTAN_MALLOC(
         part_dim * (zz->LB.Num_Global_Parts) * sizeof(float));
    if (zz->LB.Part_Sizes == NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      error = ZOLTAN_MEMERR;
      goto End;
    }
    for (i=0; i< part_dim * (zz->LB.Num_Global_Parts); i++){
      zz->LB.Part_Sizes[i] = -1.0;
    }
  }
  else if ((part_dim > zz->LB.Max_Part_Dim) ||
           (zz->LB.Num_Global_Parts > zz->LB.Max_Global_Parts)){

    /* Allocate larger LB.part_sizes array. */
    temp = zz->LB.Part_Sizes;
    zz->LB.Part_Sizes = (float *) ZOLTAN_MALLOC(
        part_dim * (zz->LB.Num_Global_Parts) * sizeof(float));
    if (zz->LB.Part_Sizes == NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      error = ZOLTAN_MEMERR;
      goto End;
    }
    /* Move data into new place. */
    for (i = 0; i < zz->LB.Num_Global_Parts; i++){
      for (j = 0; j < part_dim; j++){
        zz->LB.Part_Sizes[i*part_dim+j] = 
          temp[i*(zz->LB.Max_Part_Dim)+j];
      }
    }
    /* Update max values, which give the dimensions of Part_Sizes. */
    zz->LB.Max_Part_Dim = part_dim;
    zz->LB.Max_Global_Parts = zz->LB.Num_Global_Parts;
  }

  /* Insert new values into ps array. */
  for (i=0; i<len; i++){
    /* Error check. */
    if (part_ids[i] >= zz->LB.Num_Global_Parts){
      sprintf(msg, "partition number %d is too high", part_ids[i]);
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
      error = ZOLTAN_WARN;
    }
    else if (wgt_idx[i] >= part_dim){
      sprintf(msg, "partition weight index %d is too high", wgt_idx[i]);
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
      error = ZOLTAN_WARN;
    }
    else /* indices OK */
      zz->LB.Part_Sizes[part_ids[i]*part_dim+wgt_idx[i]] = part_sizes[i];
  }

End:
  return error;
}



int Zoltan_LB_Get_Part_Sizes(ZZ *zz, 
    int part_dim, int num_global_parts, float *part_sizes)
{
/*
 *  Function to get the scaled partition sizes.
 *
 *  Input:
 *    zz            --  The Zoltan structure to which this method
 *                      applies.
 *    part_dim      --  The number of object weights per partition.
 *                      (This usually equals lb->Obj_Wgt_Dim.)
 *    num_global_parts -- Number of global partitions.
 *                      (This usually equals lb->Num_Global_Parts)
 *
 *  Output:
 *    part_sizes    --  Array of floats that gives the set partition 
 *                      sizes, scaled such that they sum to one.
 */
  int i, j;
  float *temp_part_sizes=NULL, *sum=NULL;
  int error = ZOLTAN_OK;
  char msg[128];
  static char *yo = "Zoltan_LB_Get_Part_Sizes";

  if (part_dim > zz->LB.Max_Part_Dim){
    sprintf(msg, "Input part_dim=%d is too large", part_dim); 
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    error = ZOLTAN_FATAL;
    goto End;
  }
  if (num_global_parts > zz->LB.Max_Global_Parts){
    sprintf(msg, "Input num_global_parts=%d is too large", num_global_parts); 
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    error = ZOLTAN_FATAL;
    goto End;
  }
  if (part_sizes == NULL){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input argument part_sizes is NULL.");
    error = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->LB.Part_Sizes == NULL){
    /* Uniform partition sizes is the default. */
    /* EBEB What if LB.Part_Sizes is NULL on some procs but not all? */
    zz->LB.Uniform_Parts = 1;
    for (i = 0; i < num_global_parts*part_dim; i++)
      part_sizes[i] = 1.0 / (float)num_global_parts;
  }
  else {
   /* Get the partition sizes set by the user (application).
    * LB.Part_Sizes is stored distributed, so first we need to
    * gather all the partition sizes using collective communication.
    * Only get the data we need (part_dim * num_global_parts).
    */
    zz->LB.Uniform_Parts = 0;
    sum = (float *)ZOLTAN_MALLOC(part_dim*sizeof(float));

    /* Pack LB.Part_Sizes into temp array if leading dimensions differ */
    if (part_dim != zz->LB.Max_Part_Dim){
      temp_part_sizes = (float *)ZOLTAN_MALLOC(num_global_parts*part_dim
        *sizeof(float));
      if (temp_part_sizes == NULL){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
        error = ZOLTAN_MEMERR;
        goto End;
      }
      for (i = 0; i < num_global_parts; i++){
        for (j = 0; j < part_dim; j++){
          temp_part_sizes[i*part_dim+j] = 
            zz->LB.Part_Sizes[i*(zz->LB.Max_Part_Dim)+j];
        }
      }
    }

    /* Reduce over all procs */
    MPI_Allreduce((void*) (temp_part_sizes ? temp_part_sizes : 
      zz->LB.Part_Sizes), (void*) part_sizes, num_global_parts*part_dim, 
      MPI_FLOAT, MPI_MAX, zz->Communicator);
  
    /* Check for errors. Scale the sizes so they sum to one for each weight. */
    for (j = 0; j < part_dim; j++) 
      sum[j] = 0.0;

    for (i = 0; i < num_global_parts; i++){
      for (j = 0; j < part_dim; j++){
        if (part_sizes[i*part_dim+j]<0){
          sprintf(msg, "Partition weight (%1d,%1d) is invalid or has not been set.", i, j); 
	  ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
          error = ZOLTAN_FATAL;
          goto End;
        }
        sum[j] += part_sizes[i*part_dim+j];
      }
    }

    /* Check for sum[j] == 0 (error). */
    for (j = 0; j < part_dim; j++) {
      if (sum[j] == 0.0) {
        sprintf(msg, "Sum of weights (component %1d) is zero.", j);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
        error = ZOLTAN_FATAL;
        goto End;
      }
    }

    /* Normalize partition sizes */
    for (i = 0; i < num_global_parts; i++)
      for (j = 0; j < part_dim; j++)
        part_sizes[i*part_dim+j] /= sum[j];

  }
 
End:
  if (temp_part_sizes) ZOLTAN_FREE(&temp_part_sizes);
  if (sum)             ZOLTAN_FREE(&sum);

  return error;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
