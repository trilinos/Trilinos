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

#include "zz_const.h"
#include "order_const.h"

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
  int i, num_obj;
  char *yo = "Zoltan_Get_Distribution";

  num_obj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  if (ierr){
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
  int order_type, 	/* Input: Local or global ordering? */
  int start_index	/* Input: Do permutations start with 0 or 1? */
  )
{
  int i, num_obj, offset;
  int *proclist, *sendlist, *recvlist;
  ZOLTAN_COMM_OBJ *comm_plan;
  char *yo = "Zoltan_Inverse_Perm";
  static int owner();

  /* num_obj = local number of objects (elements in the perm vectors) */
  num_obj = vtxdist[(zz->Proc)+1] - vtxdist[zz->Proc];

  /* Verify that input permutation is really a permutation */
  /* Also check that start_index is correct. */

  /* Convert permutation vector to 0-base if necessary */
  if (start_index>0){
    for (i=0; i<num_obj; i++)
      perm[i] -= start_index;
  }

  if (order_type == ZOLTAN_LOCAL){
    /* Local inverse */
    for (i=0; i<num_obj; i++)
      inv_perm[perm[i]] = i;
  }
  else if (order_type == ZOLTAN_GLOBAL){
    /* Global inverse; use Zoltan Comm package */
    proclist = (int *) ZOLTAN_MALLOC (5*num_obj*sizeof(int));
    sendlist = &proclist[num_obj];
    recvlist = &proclist[3*num_obj];
    /* Set up comm plan. We know where to send. */
    /* Send pairs of (i, perm[i]) to other procs */
    offset = vtxdist[zz->Proc];
    for (i=0; i<num_obj; i++){
      sendlist[i] = perm[i];
      proclist[2*i] = offset+i;
      proclist[2*i+1] = owner(vtxdist, zz->Num_Proc, perm[i]);
    }
    ierr = Zoltan_Comm_Create(&comm_plan, num_obj, proclist, comm, TAG1, &nrecv);
    if (nrecv != num_obj){
      /* This should never happen. */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Internal error: nrecv != num_obj. Permutation may be invalid.");
    }
    /* Do the communication. */
    ierr = Zoltan_Comm_Do(comm_plan, TAG2, (char *)sendlist, 2*sizeof(int), (char *) recvlist);

    /* Permute data locally. */
    for (i=0; i<num_obj; i++){
      /* inv_perm[perm[i]] = i; */
      inv_perm[recvlist[2*i+1]-offset] = recvlist[2*i];
    }
    /* Free the comm_plan, proclist, sendlist, recvlist. */
    Zoltan_Comm_Destroy( &comm_plan);
    ZOLTAN_FREE(proclist);
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unknown order_type.");
  }

  /* Convert permutation vectors back to their right start_index */
  if (start_index>0){
    for (i=0; i<num_obj; i++){
      perm[i] += start_index;
      inv_perm[i] += start_index;
    }
  }

  return ZOLTAN_OK;
}


/* Find out which proc owns a certain index by binary search */
static int owner(int *vtxdist, int p, int index)
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
  return lo;
}

