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

#ifndef __ZOLTAN_MATRIX_INPUT
#define __ZOLTAN_MATRIX_INPUT


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"

static PARAM_VARS MM_params[] = {
  /* Add parameters here. MM_PARAMETER is a dummy. */
  {"MM_PARAMETER",     NULL,  "INT",     0},
  {NULL,               NULL,  NULL,     0}
};
    
/***********************************************************************
** Allow user of parallel hypergraph methods to input a sparse matrix, 
** and have Zoltan create the hypergraph that represents communication
** required for matrix/vector multiplication (y = Ax).  Then call
** Zoltan's PHG code to do the work.
**
** Zoltan then can return a partitioning of non-zeros (of A), a
** partitioning of rows (for y) and a partitioning of columns (for x).
** (Details to be worked out yet.)
**
** See Erik Boman's message to Isorropia-developers at
** http://software.sandia.gov/pipermail/isorropia-developers/2006-August/000065.html
**
** This avoids the necessity of having the user figure out how our
** hypergraph partitioning codes work.
**
** We may want to change the name of Zoltan_Matrix_Multiply, since
** it sounds like we are performing a MM instead of partitioning data
** to balance one.
*************************************************************************/

static int process_matrix_input(ZZ *zz, ZOLTAN_MM_DATA *mmd);
static int make_hypergraph(ZZ *zz, ZOLTAN_MM_DATA *mmd);
static ZOLTAN_MM_DATA *MM_Initialize_Structure();


void Zoltan_MM_Free_Structure(ZZ *zz)
{
  ZOLTAN_MM_DATA *mmd = (ZOLTAN_MM_DATA *)zz->LB.Data_Structure;

  if (mmd != NULL){

    /* TODO release of any structures allocated in mmd */

    Zoltan_Destroy(&(mmd->zzLib));  /* we created this PHG problem */

    ZOLTAN_FREE(&zz->LB.Data_Structure);
  }
}

int Zoltan_MM_Initialize_Params(ZZ *zz, int *mmval)
{
  int ierr = ZOLTAN_OK;
  int mmparam;

  Zoltan_Bind_Param(MM_params, "MM_PARAMETER", &mmparam);  
  mmparam = 0;  /* default */
  ierr = Zoltan_Assign_Param_Vals(zz->Params, MM_params, zz->Debug_Level,
          zz->Proc, zz->Debug_Proc);

  if (ierr == ZOLTAN_OK){
    *mmval = mmparam;
  }
  return ierr;
}

#define ROW_TYPE 1
#define COL_TYPE 2

int Zoltan_Matrix_Multiply(
ZZ *zz,                    /* The Zoltan structure  */
float *part_sizes,         /* Input:  Array of size zz->Num_Global_Parts
                                containing the percentage of work assigned
                                to each partition. */
ZOLTAN_MM_RESULT **results) /* Output: undefined as of yet */
{
  char *yo = "Zoltan_Matrix_Multiply";
  int ierr = ZOLTAN_OK;
  ZZ *zzLib = NULL;
  int npins;

  ZOLTAN_MM_DATA *mmd = MM_Initialize_Structure();

  if (mmd == NULL){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* Process parameters.  We don't have any yet.
  ierr = Zoltan_MM_Initialize_Params(ZZ *zz, pointers to parameters);
   */

  /* Call application defined query functions to get the non-zeros.
   * User has a choice of giving us compressed rows or compressed columns.
   * We need global row and column IDs.  Any process can give us any of the
   * non-zeros, but each non-zero must be supplied by only one process.
   *
   * I'm assuming row and column IDs are long ints, not more general objects.
   * We know we are working with a matrix, not something abstracted as
   * a matrix.  I'm assuming the IDs are contiguous as well.
   */

  ierr = ZOLTAN_FATAL;

  if ((zz->Get_CSC_Size != NULL) && (zz->Get_CSC != NULL)){
    ierr = ZOLTAN_OK;
    mmd->input_type = COL_TYPE;
  }
  else if ((zz->Get_CSR_Size != NULL) && (zz->Get_CSR != NULL)){
    mmd->input_type = ROW_TYPE;
    ierr = ZOLTAN_OK;
  }

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Required query functions are not defined.\n");
    goto End;
  }

  if (mmd->input_type == COL_TYPE){
    zz->Get_CSC_Size(zz->Get_CSC_Size_Data, &mmd->nRC, &npins, &ierr);
  }
  else {
    zz->Get_CSR_Size(zz->Get_CSR_Size_Data, &mmd->nRC, &npins, &ierr);
  }

  mmd->rcGID = (long int *)ZOLTAN_MALLOC(mmd->nRC * sizeof(long int));
  mmd->pinIndex = (long int *)ZOLTAN_MALLOC((mmd->nRC+1) * sizeof(long int));
  mmd->pinGID = (long int *)ZOLTAN_MALLOC(npins * sizeof(long int));
  if ((mmd->nRC && !mmd->rcGID) || !mmd->pinIndex || (npins && !mmd->pinGID)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  if (mmd->input_type == COL_TYPE){
    zz->Get_CSC(zz->Get_CSC_Data, mmd->nRC, npins, 
                mmd->rcGID, mmd->pinIndex, mmd->pinGID, &ierr);
  }
  else {
    zz->Get_CSR(zz->Get_CSR_Data, mmd->nRC, npins, 
                mmd->rcGID, mmd->pinIndex, mmd->pinGID, &ierr);
  }
 
  mmd->pinIndex[mmd->nRC] = npins;

  /*
   * Determine globally how many rows and columns were returned, and if
   * their IDs are 0 based or 1 based.  (Assume they are contiguous.)
   */

  ierr = process_matrix_input(zz, mmd);

  /*
   * Convert the matrix to a hypergraph.  Set the query functions we will use
   * with PHG.  We will call one routine to do
   * this, but eventually we could have more than one way to do this.
   */

  ierr = make_hypergraph(zz, mmd);

  /* Update our LB.Data_Structure with everything we need to respond to queries
   * about our vertices and hyperedges.
   * 
   * Update any PHG parameters that need to be set, including PHG queries.
   */

  zzLib = Zoltan_Create(zz->Communicator);
  Zoltan_Set_Param(zzLib, "LB_METHOD", "HYPERGRAPH"); 
  Zoltan_Set_Param(zzLib, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param(zzLib, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zzLib, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zzLib, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zzLib, "RETURN_LISTS", "ALL");    /* not sure we need "all" */
  Zoltan_Set_Param(zzLib, "OBJ_WEIGHT_DIM", "1"); 
  Zoltan_Set_Param(zzLib, "EDGE_WEIGHT_DIM", "0");
  Zoltan_Set_Param(zzLib, "ADD_OBJ_WEIGHT", "pins");

  mmd->zzLib = zzLib;

  /*
   * Call Zoltan_PHG to partition the rows and columns (a.k.a. vertices)

  ierr = Zoltan_PHG(zzLib, part_sizes, etc.);
   */

  /*
   * "Partition" the non-zeros that were returned by user in CSR or CSC in some
   * way that acheives balance and minimizes row and column cuts.  (Most
   * likely assigning to the process that either their row or their column
   * went to.)  Save info required so that when user calls 
   * Zoltan_Get_Sparse_Matrix_NonZero_Partition_Assigned(i, j) for their
   * non-zero, we can answer it.  (Maybe it should have a different name...)
   */

  /* 
   * Write to a results structure either the global row/column partitionings
   * or a way to get them (function pointer).  Write the local non-zero
   * partitionings or a way to get them (function pointer).  Maybe row/column
   * global information will be too big, and just save import/export lists.
   */

End:
  return ierr;
}

static int make_hypergraph(ZZ *zz, ZOLTAN_MM_DATA *mmd)
{
  int ierr = ZOLTAN_OK;

  return ierr;
}

static int process_matrix_input(ZZ *zz, ZOLTAN_MM_DATA *mmd)
{
  int ierr = ZOLTAN_OK;
  int i;
  long int minID, maxID, minPinID, maxPinID;
  long int npins=0;
  long int vals[2], gvals[2];

  /* get global number of rows, columns and pins, range of row and column IDs */

  if (mmd->nRC > 0)
    npins = mmd->pinIndex[mmd->nRC];

  MPI_Allreduce(&npins, &mmd->nNonZeroes, 1, MPI_LONG, MPI_SUM, zz->Communicator);

  maxID = maxPinID = -1;

  if (mmd->nRC > 0){
    minID = maxID = mmd->rcGID[0];

    for (i=1; i<mmd->nRC; i++){
      if (mmd->rcGID[i] < minID) minID = mmd->rcGID[i];
      else if (mmd->rcGID[i] > maxID) maxID = mmd->rcGID[i];
    }
    if (npins > 0){
      minPinID = maxPinID = mmd->pinGID[0];
      for (i=1; i<npins; i++){
        if (mmd->pinGID[i] < minPinID) minPinID = mmd->pinGID[i];
        else if (mmd->pinGID[i] > maxPinID) maxPinID = mmd->pinGID[i];
      }
    }
  }
  vals[0] = maxID;
  vals[1] = maxPinID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MAX, zz->Communicator);

  maxID = gvals[0];
  maxPinID = gvals[1];

  if (npins == 0){
    minPinID = maxPinID;
    if (mmd->nRC == 0){
      minID = maxID;
    }
  }
  vals[0] = minID;
  vals[1] = minPinID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MIN, zz->Communicator);

  minID = gvals[0];
  minPinID = gvals[1];

  if (mmd->input_type == ROW_TYPE){
    mmd->rowBaseID = minID;
    mmd->colBaseID = minPinID;
    mmd->nRows = maxID - minID + 1;
    mmd->nCols = maxPinID - minPinID + 1;
  }
  else{
    mmd->rowBaseID = minPinID;
    mmd->colBaseID = minID;
    mmd->nRows = maxPinID - minPinID + 1;
    mmd->nCols = maxID - minID + 1;
  }

  /*
   * Here we could check the input.  Make sure that row and column IDs
   * are contigous values starting at the baseID.
   */

  return ierr;
}
static ZOLTAN_MM_DATA *MM_Initialize_Structure()
{
  /* 
   * This is the structure we save at zz->LB.Data_Structure
   */
  ZOLTAN_MM_DATA *mmd = (ZOLTAN_MM_DATA *)ZOLTAN_MALLOC(sizeof(ZOLTAN_MM_DATA));

  if (mmd == NULL){
    return NULL;
  }

  mmd->input_type = 0;
  mmd->nRC = 0;
  mmd->rcGID = NULL;
  mmd->pinIndex = NULL;
  mmd->pinGID = NULL;

  mmd->rowBaseID=0;
  mmd->colBaseID=0;
  mmd->nRows=0;
  mmd->nCols=0;
  mmd->nNonZeroes=0;

  mmd->zzLib = NULL;

  return mmd;
}

int Zoltan_Matrix_Multiply_Free(ZZ *zz)
{
  /* Free our structures */
  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_MATRIX_INPUT */
