// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "order_const.h"

#include <stdio.h>

/* MPI tags */
#undef TAG1
#undef TAG2
#define TAG1 32111
#define TAG2 32112


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for manipulating permutations in Zoltan.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 * Compute an array that contains the cumulative sum of objects
 * on each processor.
 *
 * Memory for the vtxdist array is allocated here,
 * but must be freed by the calling routine.
 *
 */
int Zoltan_Get_Distribution(ZZ *zz, int **vtxdist)
{
  int ierr = ZOLTAN_OK, num_obj;
  char *yo = "Zoltan_Get_Distribution";

  num_obj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
    /* Return error code */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Get_Num_Obj.");
    return (ierr);
  }
  
  *vtxdist = (int *) ZOLTAN_MALLOC((zz->Num_Proc+1)*sizeof(int));
  if (num_obj>0){
    if (!(*vtxdist)){
      /* Not enough memory */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory.");
      return ZOLTAN_MEMERR;
    }
  }
  
  /* Construct *vtxdist[i] = the number of objects on all procs < i. */
  /* Scan to compute partial sums of the number of objs */
  MPI_Scan (&num_obj, *vtxdist, 1, MPI_INT, MPI_SUM, zz->Communicator);
  /* Gather data from all procs */
  MPI_Allgather (&((*vtxdist)[0]), 1, MPI_INT, 
                 &((*vtxdist)[1]), 1, MPI_INT, zz->Communicator);
  (*vtxdist)[0] = 0;
  
  return ZOLTAN_OK;
}


int Zoltan_Inverse_Perm(
  ZZ  *zz,		/* Input: Zoltan struct */
  int *perm, 		/* Input: Permutation to invert. */
  int *inv_perm, 	/* Output: Inverse permutation of perm. */
  int *vtxdist, 	/* Input: Distribution of the vectors. */
  char *order_type, 	/* Input: Local or global ordering? */
  int start_index	/* Input: Do permutations start with 0 or 1? */
  )
{
  int i, ierr, num_obj, nrecv, offset;
  int *proclist, *sendlist, *recvlist;
  ZOLTAN_COMM_OBJ *comm_plan;
  char msg[256];
  char *yo = "Zoltan_Inverse_Perm";

  ierr = ZOLTAN_OK;
  proclist = sendlist = recvlist = NULL;
  comm_plan = NULL;

  /* num_obj = local number of objects (elements in the perm vectors) */
  num_obj = vtxdist[(zz->Proc)+1] - vtxdist[zz->Proc];

  /* Verify that input permutation is really a permutation */
  /* Also check that start_index is correct. */

  /* Convert permutation vector to 0-base if necessary */
  if (start_index>0){
    for (i=0; i<num_obj; i++)
      perm[i] -= start_index;
  }

  if (strcmp(order_type, "LOCAL")==0){
    /* Local inverse */
    for (i=0; i<num_obj; i++)
      inv_perm[perm[i]] = i;
  }
  else if (strcmp(order_type, "GLOBAL")==0){
    /* Global inverse; use Zoltan Comm package */
    proclist = (int *) ZOLTAN_MALLOC (5*num_obj*sizeof(int));
    sendlist = &proclist[num_obj];
    recvlist = &proclist[3*num_obj];
    /* Set up comm plan. We know where to send. */
    /* Send pairs of (i, perm[i]) to other procs */
    offset = vtxdist[zz->Proc];
    for (i=0; i<num_obj; i++){
      sendlist[2*i] = offset+i;
      sendlist[2*i+1] = perm[i];
      proclist[i] =  Zoltan_Get_Processor_Graph(vtxdist, zz->Num_Proc, perm[i]);
    }
    ierr = Zoltan_Comm_Create(&comm_plan, num_obj, proclist, 
             zz->Communicator, TAG1, &nrecv);
    if (ierr != ZOLTAN_OK){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_Comm_Create");
      goto error;
    }
    if (nrecv != num_obj){
      /* This should never happen. */
      sprintf(msg, "Internal error: nrecv (%3d) != num_obj (%3d). Invalid permutation.\n", nrecv, num_obj); 
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      ierr = ZOLTAN_FATAL;
      goto error;
    }
    /* Do the communication. */
    ierr = Zoltan_Comm_Do(comm_plan, TAG2, (char *)sendlist, 2*sizeof(int), (char *) recvlist);
    if (ierr != ZOLTAN_OK){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_Comm_Do");
      goto error;
    }

    /* Permute data locally. */
    for (i=0; i<num_obj; i++){
      /* inv_perm[perm[i]] = i; */
      inv_perm[recvlist[2*i+1]-offset] = recvlist[2*i];
    }
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unknown order_type.");
    ierr = ZOLTAN_FATAL;
    goto error;
  }

  /* Convert permutation vectors back to their right start_index */
  if (start_index>0){
    for (i=0; i<num_obj; i++){
      perm[i] += start_index;
      inv_perm[i] += start_index;
    }
  }

error:
  /* Free the comm_plan, proclist, sendlist, recvlist. */
  if (comm_plan) Zoltan_Comm_Destroy( &comm_plan);
  if (proclist ) ZOLTAN_FREE(&proclist);

  return (ierr);
}


/* Find out which proc owns a certain index by binary search */
int Zoltan_Get_Processor_Graph(int *vtxdist, int p, int index)
{
  int lo, hi, mid;

  lo = 0;
  hi = p;

  /* Check for values out of range */
  if (index<vtxdist[0] || index>=vtxdist[p])
    return (-1);

  /* Binary search */
  while(hi > lo+1){
    mid = (lo+hi)/2;
    if (index > vtxdist[mid])
      lo = mid;
    else if (index < vtxdist[mid])
      hi = mid;
    else /* index == vtxdist[mid] */
      return mid;
  }
  /* printf("DEBUG: proc %d owns object %d\n", lo, index); */

  return lo;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
