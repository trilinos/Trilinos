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


int Zoltan_Divide_Machine(
   ZZ *zz,             /* The Zoltan structure (not used now, will be
                          used for pointer to machine details */
   int obj_wgt_dim,    /* Number of different weights (loads). */
   float *part_sizes,  /* Array of partition sizes, containing percentage of 
                          work per partition. (length= obj_wgt_dim*num_parts) */
   int proc,           /* my processor number in global sense */
   MPI_Comm comm,      /* communicator for part of machine to be divided */
   int *set,           /* set that proc is in after divide (lowest global
                          numbered processor in set 0) */
   int *proclower,     /* lowest numbered processor in first set */
   int *procmid,       /* lowest numbered processor in second set */
   int *num_procs,     /* on input, # of procs to be divided
                          on exit, # of procs in the set that proc is in */
   int *partlower,     /* lowest numbered partition in first set */
   int *partmid,       /* lowest numbered partition in second set */
   int *num_parts,     /* on input, # of partitions to be divided
                          on exit, # of parts in the set that proc is in */
   double *fractionlo  /* actual division of machine: % of work to be assigned
                          to first set (length obj_wgt_dim) */
)
{
int i, j, k;
int np = 0;     /* Number of partitions on procmid */
int fpartmid;   /* First partition on procmid */
int totalparts; /* Total number of partitions in input set. */
int totalprocs; /* Total number of processors in input set. */
int dim = obj_wgt_dim;
double *sum = NULL;

/* This routine divides the current machine (defined by the communicator)
 * into two pieces.
 * For now, it simply divides the machine in half.  In the future, it will
 * be a more complicated routine taking into account the architecture of
 * the machine and communication network. 
 * The two resulting sets contain contiguously numbered processors 
 * and partitions.
 */

  if (dim<1) dim = 1;   /* In case obj_wgt_dim==0. */

  /* The following statement assumes that proclower is being set correctly in
     the calling routine if Tflops_Special flag is set */
  if (!zz->Tflops_Special)
     MPI_Allreduce(&proc, proclower, 1, MPI_INT, MPI_MIN, comm);

  totalparts = *partlower + *num_parts;
  totalprocs = *proclower + *num_procs;

  /* Compute procmid as roughly half the number of processors. */
  /* Then partmid is the lowest-numbered partition on procmid. */

  *procmid = *proclower + (*num_procs - 1)/2 + 1;
  if (*procmid < totalprocs)
    Zoltan_LB_Proc_To_Part(zz, *procmid, &np, &fpartmid);
  if (np > 0)
    *partmid = fpartmid;
  else {
    /* No partitions on procmid; find next part number in procs > procmid */
    i = *procmid;
    while (np == 0 && (++i) < totalprocs) {
      Zoltan_LB_Proc_To_Part(zz, i, &np, &fpartmid);
    }
    if (np) 
      *partmid = fpartmid;
    else
      *partmid = totalparts;
  }

  /* Check special cases */

  if (!zz->LB.Single_Proc_Per_Part && *partmid != totalparts) {
    i = Zoltan_LB_Part_To_Proc(zz, *partmid, NULL);
    if (i != *procmid) {

      /* Partition is spread across several processors. 
         Don't allow mid to fall within a partition; reset procmid so that it
         falls at a partition boundary.  */

      if (i != *proclower) {
        /* set procmid to lowest processor containing partmid */
        *procmid = i;
      }
      else { /* i == *proclower */
        /* Move mid to next partition so that procmid != proclower */
        (*partmid)++;
        *procmid = Zoltan_LB_Part_To_Proc(zz, *partmid, NULL);
      }
    }
  }

  /* Sum up desired partition sizes. */
  sum = (double *)ZOLTAN_MALLOC(dim*sizeof(double));

  for (k=0; k<dim; k++){
    sum[k] = 0.0;
    fractionlo[k] = 0.0;
  }
  for (i = 0; i < *num_parts; i++) {
    j = *partlower + i;
    for (k=0; k<dim; k++){
      if (j < *partmid)
        fractionlo[k] += (double) part_sizes[j*dim+k];
      sum[k] += (double) part_sizes[j*dim+k];
    }
  }
  for (k=0; k<dim; k++)
    if (sum[k] != 0.0) fractionlo[k] /= sum[k];

  if (proc < *procmid) {
    *set = 0;
    *num_parts = *partmid - *partlower;
    *num_procs = *procmid - *proclower;
  } 
  else {
    *set = 1;
    *num_parts = totalparts - *partmid;
    *num_procs = totalprocs - *procmid;
  }

  ZOLTAN_FREE(&sum);
  return ZOLTAN_OK;
}

int Zoltan_Divide_Parts(
   ZZ *zz,             /* The Zoltan structure (not used now, will be
                          used for pointer to machine details */
   int obj_wgt_dim,    /* Number of different weights (loads). */
   float *part_sizes,  /* Array of partition sizes, containing percentage of 
                          work per partition. (length= obj_wgt_dim*num_parts) */
   int num_parts,      /* Input: # of partitions to be divided */
   int *partlower,     /* lowest numbered partition in first set */
   int *partmid,       /* lowest numbered partition in second set */
   double *fractionlo  /* actual division of machine: % of work to be assigned
                          to first set (array if obj_wgt_dim>1) */
)
{
int i, j, k;
int dim = obj_wgt_dim;
double *sum = NULL;

/* This SERIAL routine divides the current group of partitions
 * into two pieces with roughly equal numbers of partitions per piece. 
 * It is designed to be used within a single processor to divide its
 * partitions into two sets (e.g., in serial_rcb).
 */

  if (obj_wgt_dim<1) dim = 1; /* In case obj_wgt_dim==0. */

  /* Compute procmid as roughly half the number of processors. */
  /* Then partmid is the lowest-numbered partition on procmid. */

  *partmid = *partlower + (num_parts - 1)/2 + 1;

  sum = (double *)ZOLTAN_MALLOC(dim*sizeof(double));

  for (k=0; k<dim; k++){
    sum[k] = 0.0;
    fractionlo[k] = 0.0;
  }
  for (i = 0; i < num_parts; i++) {
    j = *partlower + i;
    for (k=0; k<dim; k++){
      if (j < *partmid)
        fractionlo[k] += (double) part_sizes[j*dim+k];
      sum[k] += (double) part_sizes[j*dim+k];
    }
  }
  for (k=0; k<dim; k++)
    if (sum[k] != 0.0) fractionlo[k] /= sum[k];

  ZOLTAN_FREE(&sum);
  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
