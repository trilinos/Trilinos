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

static PARAM_VARS MATRIX_MULTIPLY_params[] = {
  /* Add parameters here. */
  {NULL,                              NULL,  NULL,     0}
};
    
/***********************************************************************
** Allow user of parallel hypergraph methods to input a sparse matrix, 
** and have Zoltan create the hypergraph that represents communication
** required for matrix/vector multiplication (y = Ax).  Then call
** Zoltan's PHG code to do the work.
**
** Zoltan then can return a partitioning of non-zeroes (of A), a
** partitioning of rows (for y) and a partitioning of columns (for x).
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

/************* these will go in a header file ***************************/
/*
** The data we will save at the user's ZZ->LB.Data_Structure
*/
struct Zoltan_MM_Data_Struct{
  ZZ *zzLib;    /* PHG problem created by Zoltan from user's sparse matrix */
}

typedef struct Zoltan_MM_Data_Struct ZOLTAN_MM_DATA;

/*
** The structure we will return to user, no clue what's in it - maybe
** just function pointers so they can query partitioning.
*/

struct Zoltan_MM_Results_Struct{
  void (*you_tell_me)();
}

typedef struct Zoltan_MM_Results)Struct ZOLTAN_MM_RESULT;

/************************************************************************/

static ZOLTAN_MM_DATA *MM_Initialize_Structure()
{
  ZOLTAN_MM_DATA *mmd = (ZOLTAN_MM_DATA *)ZOLTAN_MALLOC(sizeof(ZOLTAN_MM_DATA));

  if (mmd == NULL){
    return NULL;
  }

  mmd->zzLib = NULL;
}

void Zoltan_MM_Free_Structure(ZZ *zz)
{
  ZOLTAN_MM_DATA *mmd = (ZOLTAN_MM_DATA *)zz->LB.Data_Structure;

  if (mmd != NULL){

    Zoltan_Destroy(&(zz->zzLib));  /* we created this PHG problem */

    ZOLTAN_FREE(&zz->LB.Data_Structure);
  }
}

int Zoltan_Matrix_Multiply(
ZZ *zz,                    /* The Zoltan structure  */
float *part_sizes,         /* Input:  Array of size zz->Num_Global_Parts
                                containing the percentage of work assigned
                                to each partition. */
ZOLTAN_MM_RESULT **results) /* Output: functions to call to get partitionings */
                                that are computed here. To be determined. */
{
  char *yo = "Zoltan_Matrix_Multiply";
  int ierr = ZOLTAN_OK;
  ZZ *zzLib = NULL;

  ZOLTAN_MM_DATA *mmd = MM_Initialize_Structure();

  if (mmd == NULL){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* Call application defined query functions to get the non-zeroes.
   * User has a choice of giving us compressed rows or compressed columns.
   * We need global row and column IDs.  Any process can give us any of the
   * non-zeroes, but each non-zero must be supplied by only one process.
   *
   * I'm assuming row and column IDs are long ints, not more general objects.
   *
   * If the ccs_size_fn is defined then the user is giving us compressed
   * columns.  If the ccr_size_fn is defined the user is giving us compressed
   * rows.
   *
   * typedef void ZOLTAN_CCS_SIZE_FN(void *data, int *num_columns,
   *                                 int *num_non_zeroes, int *ierr)
   *
   * typedef void ZOLTAN_CRS_SIZE_FN(void *data, int *num_rows,
   *                                 int *num_non_zeroes, int *ierr)
   *
   * typedef void ZOLTAN_CCS_FN(void *data, int num_columns, int num_non_zeroes,
   *         long int *column_gids, long int *row_gid_index, long int *row_gids,
   *         int *ierr)
   *
   * typedef void ZOLTAN_CRS_FN(void *data, int num_rows, int num_non_zeroes,
   *         long int *row_gids, long int *col_gid_index, long int *col_gids,
   *         int *ierr)
   */

  /* Hash rows to processes, columns to processes (do we need to worry about
   * balance, or just divide them up?)
   */

  /* Create hypergraph and query functions.  This is a simple way to do
   * it and may change.
   *  
   * Each row and column becomes a "vertex".  (So M rows, N columns yields
   * M+N vertices).  We need to remember which row or column each vertex
   * originally represented.  We may want vertex weight to be pin weight.
   * (So non-zeroes are balanced later on if they go to owner of their row
   * or column.)
   *
   * Create hyperedges for your rows and your columns, don't know if we should
   * have edge weights.
   */

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
   */

  ierr = Zoltan_PHG(zzLib, part_sizes, etc.);

  /*
   * "Partition" the non-zeroes that were returned by user in CRS or CCS in some
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

  /*
   * End:
   */
}

int Zoltan_Matrix_Multiply_Free(ZZ *zz)
{
  /* Free our structures */
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_MATRIX_INPUT */
