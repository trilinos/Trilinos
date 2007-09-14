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

#ifndef __ZOLTAN_MATRIX_PARTITION
#define __ZOLTAN_MATRIX_PARTITION

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*
 * Subroutines to support sparse matrix partitioning.
 * Currently, only 1d (row or column) methods are supported,
 * but the design allows general 2d methods.
 */
#include "matrix_partition.h"
/*******************/

#include "phg.h"   /* for Zoltan_Convert_To_CSR */
#include "zz_util_const.h"   /* for Zoltan_Hash */

/* Query functions set up by Zoltan_Matrix_Partition for the
 * Zoltan library
 */

/* Parameters for Zoltan method SPARSE_MATRIX
 */
static PARAM_VARS MP_params[] = {
  {"MATRIX_APPROACH",         NULL,  "STRING",     0},
  {"MATRIX_METHOD",         NULL,  "STRING",     0},
  {NULL,                  NULL,  NULL,     0}
};

int Zoltan_MP_Set_Param(
  char *name,                     /* name of variable */
  char *val)                      /* value of variable */
{
  /* associates value to named variable for MP partitioning parameters */
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */
  int status;
  int i;

  char *valid_approach[] = {
        "ROWS", "ROW", 
        "COLS", "COLUMNS", "COLUMN",
        "GENERAL", "NONZEROS", "NONZERO",
         NULL };

  status = Zoltan_Check_Param(name, val, MP_params, &result, &index);

  if (status == 0){
    if (strcasecmp(name, "MATRIX_APPROACH") == 0){
      status = 2;
      for (i=0; valid_approach[i] != NULL; i++){
        if (strcasecmp(val, valid_approach[i]) == 0){
          status = 0;
          break;
        }
      }
    }
  }
  return(status);
}


    
/***********************************************************************
** Allow user of parallel hypergraph methods to input a sparse matrix, 
** and have Zoltan create the hypergraph that best solves their problem.
**
** See Erik Boman's message to Isorropia-developers at
** http://software.sandia.gov/pipermail/isorropia-developers/2006-August/000065.html
**
** This avoids the necessity of having the user figure out how our
** hypergraph partitioning codes work.
**
** An IJTYPE is the type of a row or column ID.
**
** MATRIX_APPROACH:
**    MP_ROWS    partition the rows (objects are rows, 
**                for hypergraph model hyperedges are columns)
**                
**    MP_COLS    partition the columns (objects are columns
**                for hypergraph model hyperedges are rows)
**
**    MP_NONZEROS    partition the nonzeros 
**
** If you want to add a new approach to Zoltan_Matrix_Partition, you
** need to write functions analogous to the 5 functions written
** for the MP_ROWS and MP_COLS approaches: mp_1d_setup, mp_1d_results,
** mp_1d_get_nzs, mp_1d_get_rows, and mp_1d_get_columns.
**
*************************************************************************/

/*
 * General matrix partition functions
 */
static int process_matrix_input(ZZ *zz, ZOLTAN_MP_DATA *mpd);
static int make_mirror(ZOLTAN_MP_DATA *mpd);
static int allocate_copy(IJTYPE **to, IJTYPE *from, IJTYPE len);

static ZOLTAN_MP_DATA *MP_Initialize_Structure(ZZ *zz);
static int MP_Initialize_Params(ZZ *zz, ZOLTAN_MP_DATA *mpd);

/*
 * Functions specifically to process approaches MP_ROWS and MP_COLS
 */
static int mp_1d_setup(ZOLTAN_MP_DATA *mpd);
static int mp_1d_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  int *export_procs, int *export_to_part);
static int mp_1d_get_nzs(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nNzs, IJTYPE *I, IJTYPE *J, int *ijProcs, int *ijParts);
static int mp_1d_get_rows(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nRows, IJTYPE *rowIDs, int *rowProcs, int *rowParts);
static int mp_1d_get_columns(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nCols, IJTYPE *colIDs, int *colProcs, int *colParts);
/*
 * Functions specifically to process approach MP_NZ. 
 */
static int mp_2d_setup(ZOLTAN_MP_DATA *mpd);
static int mp_2d_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  int *export_procs, int *export_to_part);
static int mp_2d_get_nzs(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nNzs, IJTYPE *I, IJTYPE *J, int *ijProcs, int *ijParts);

/*
 * Functions to create a search structure to locate rows, columns
 * or non-zeros in the ZOLTAN_MP_DATA structure.  To support
 * sparse matrix queries that can be called when Zoltan_Matrix_Partition
 * has completed.  The lookup function Zoltan_Lookup_Obj is global.
 */
static void free_obj_lookup_table(obj_lookup **lu);
static obj_lookup *create_obj_lookup_table(IJTYPE numObjs, IJTYPE *objGIDs);
static obj_lookup *create_obj_lookup_table2(IJTYPE sizeI, IJTYPE sizeJ,
                   IJTYPE *listI, IJTYPE *indexJ, IJTYPE *listJ);

/****************************************************************/
/****************************************************************/
/* API:
 *
 *  Zoltan_Matrix_Partition()     must be called by all processes
 *
 *  Queries (not all queries are defined for all MATRIX_APPROACH):
 *
 *  Zoltan_MP_Get_Row_Assignment()
 *  Zoltan_MP_Get_Column_Assignment()
 *  Zoltan_MP_Get_NonZero_Assignment()
 *
 *  To evaluate partitioning:
 *
 *  Zoltan_Matrix_Partition_Eval()   must be called by all procs
 */

int Zoltan_MP_Get_Row_Assignment(ZZ *zz, int nRows, IJTYPE *rowIDs,
        int *rowProcs, int *rowParts)
{
  char *yo = "Zoltan_MP_Get_Row_Assignment";
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;
  int ierr = ZOLTAN_FATAL;

  /*
   * Row assignment depends on the MATRIX_APPROACH.  We can only
   * return row assignments for rows that this process gave us
   * non-zeroes for.
   */

  if (mpd){
    switch (mpd->approach)
    {
      case MP_ROW_TYPE:
      case MP_COLUMN_TYPE:
        ierr = mp_1d_get_rows(mpd, nRows, rowIDs, rowProcs, rowParts);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
    }
  }

  return ierr;
}
int Zoltan_MP_Get_Column_Assignment(ZZ *zz, int nCols, IJTYPE *colIDs,
        int *colProcs, int *colParts)
{
  char *yo = "Zoltan_MP_Get_Column_Assignment";
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;
  int ierr = ZOLTAN_FATAL;

  /*
   * Column assignment depends on the MATRIX_APPROACH.  We can only
   * return column assignments for columns that this process gave us
   * non-zeroes for.
   */

  if (mpd){
    switch (mpd->approach)
    {
      case MP_ROW_TYPE:
      case MP_COLUMN_TYPE:
        ierr = mp_1d_get_columns(mpd, nCols, colIDs, colProcs, colParts);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
    }
  }

  return ierr;
}
int Zoltan_MP_Get_NonZero_Assignment(ZZ *zz, int nNZ, 
        IJTYPE *rowIDs, IJTYPE *colIDs, int *nzProcs, int *nzParts)
{
  char *yo = "Zoltan_MP_Get_NonZero_Assignment";
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;
  int ierr = ZOLTAN_FATAL;

  /*
   * Nonzero assignment depends on the MATRIX_APPROACH.  We can only
   * return nonzero assignments for nonzeroes that were returned
   * by this process in the query functions.
   */

  if (mpd){
    switch (mpd->approach)
    {
      case MP_ROW_TYPE:
      case MP_COLUMN_TYPE:
        ierr = mp_1d_get_nzs(mpd, nNZ, rowIDs, colIDs, nzProcs, nzParts);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
    }
  }

  return ierr;
}
int Zoltan_Matrix_Partition(ZZ *zz)
{
  char *yo = "Zoltan_Matrix_Partition";
  int ierr = ZOLTAN_OK;
  ZZ *zzLib = NULL;
  IJTYPE nnz;

  /* The persistent structure we will save at zz->LB.Data_Structure */

  ZOLTAN_MP_DATA *mpd = MP_Initialize_Structure(zz);

  if (mpd == NULL){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  if (zz->LB.Data_Structure){
    zz->LB.Free_Structure(zz);
  }
  zz->LB.Data_Structure = mpd;

  /* Process parameters, write to mpd.  */

  ierr = MP_Initialize_Params(zz, mpd);

  if (ierr != ZOLTAN_OK){
    goto End;
  }

  /*
   * Create the Zoltan problem that we will define from the input sparse matrix.
   */

  zzLib = Zoltan_Create(zz->Communicator);
  mpd->zzLib = zzLib;

  /* Call application defined query functions to get the non-zeros.
   * User has a choice of giving us compressed rows or compressed columns.
   * We need global row and column IDs.  Any process can give us any of the
   * non-zeros, but each non-zero must be supplied by only one process.
   */

  ierr = ZOLTAN_FATAL;

  if ((zz->Get_CSC_Size != NULL) && (zz->Get_CSC != NULL)){
    ierr = ZOLTAN_OK;
    mpd->input_type = COL_TYPE;
  }
  else if ((zz->Get_CSR_Size != NULL) && (zz->Get_CSR != NULL)){
    mpd->input_type = ROW_TYPE;
    ierr = ZOLTAN_OK;
  }

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Required query functions are not defined.\n");
    goto End;
  }

  if (mpd->input_type == COL_TYPE){
    zz->Get_CSC_Size(zz->Get_CSC_Size_Data, &mpd->numRC, &nnz, &ierr);
  }
  else {
    zz->Get_CSR_Size(zz->Get_CSR_Size_Data, &mpd->numRC, &nnz, &ierr);
  }

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "query function failed\n");
    goto End;
  }

  mpd->rcGID = (IJTYPE *)ZOLTAN_MALLOC(mpd->numRC * sizeof(IJTYPE));
  mpd->nzIndex = (IJTYPE *)ZOLTAN_MALLOC((mpd->numRC+1) * sizeof(IJTYPE));
  mpd->nzGID = (IJTYPE *)ZOLTAN_MALLOC(nnz * sizeof(IJTYPE));
  if ((mpd->numRC && !mpd->rcGID) || !mpd->nzIndex || (nnz && !mpd->nzGID)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  if (mpd->input_type == COL_TYPE){
    zz->Get_CSC(zz->Get_CSC_Data, mpd->numRC, nnz, 
                mpd->rcGID, mpd->nzIndex, mpd->nzGID, &ierr);
  }
  else {
    zz->Get_CSR(zz->Get_CSR_Data, mpd->numRC, nnz, 
                mpd->rcGID, mpd->nzIndex, mpd->nzGID, &ierr);
  }

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "query function 2 failed\n");
    goto End;
  }
 
  mpd->nzIndex[mpd->numRC] = nnz;

  /*
   * Determine globally how many rows and columns were returned, and if
   * their IDs are 0 based or 1 based.  (Assume they are contiguous.)
   * Some row/column numbers may not appear in input if they have no nzs.
   */

  ierr = process_matrix_input(zz, mpd);

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "process_matrix_input\n");
    goto End;
  }

  /*
   * Set the query functions and parameters required for this problem.
   * Precompute hypergraph if necessary.  This is the code that needs
   * to be written if you are implementing a new MATRIX_APPROACH.
   */
  switch (mpd->approach)
    {
      case MP_ROW_TYPE:
      case MP_COLUMN_TYPE:
        ierr = mp_1d_setup(mpd);
        if (ierr != ZOLTAN_OK){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "setup for partitioning\n");
        }
        break;
      case MP_NZ_TYPE:
        ierr = mp_2d_setup(mpd);
        if (ierr != ZOLTAN_OK){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "setup for partitioning\n");
        }
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
        ierr = ZOLTAN_FATAL;
    }

  if (ierr != ZOLTAN_OK){
    goto End;
  }

  /*
   * Call Zoltan_LB_Partition to partition the objects
   */

  int changes, num_gid_entries, num_lid_entries, num_import, num_export;
  unsigned int *import_global_ids, *import_local_ids;
  unsigned int *export_global_ids, *export_local_ids;
  int *import_to_part, *export_to_part;
  int *import_procs, *export_procs;

  ierr = Zoltan_LB_Partition(zzLib,
      &changes, &num_gid_entries, &num_lid_entries, 
      &num_import,
      &import_global_ids, &import_local_ids, &import_procs, &import_to_part,
      &num_export,
      &export_global_ids, &export_local_ids, &export_procs, &export_to_part);

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "partitioning\n");
    goto End;
  }

  /*
   * Save the data required to respond to the queries that are
   * supported by each MATRIX_APPROACH.
   */

  switch (mpd->approach)
    {
      case MP_ROW_TYPE:
      case MP_COLUMN_TYPE:
        /*
         * Save row or column assignment for all rows or columns
         * in my non-zeroes.  Each of my non-zeroes is assigned to
         * the row or column that it's row or column is assigned to.
         *
         * MP_ROW_TYPE - we only give row and nonzero assignments
         * MP_COLUMN_TYPE - we only give column and nonzero assignments
         */
        ierr = mp_1d_result(mpd, num_export,
            export_global_ids, export_local_ids, export_procs, export_to_part);
        if (ierr != ZOLTAN_OK){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo,"processing results of partitioning\n");
        }
        break;
      case MP_NZ_TYPE:
        /*
         * Save nonzero assignment for all local nonzeros.
         */
        ierr = mp_2d_result(mpd, num_export,
            export_global_ids, export_local_ids, export_procs, export_to_part);
        if (ierr != ZOLTAN_OK){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo,"processing results of partitioning\n");
        }
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
        ierr = ZOLTAN_FATAL;
    }

  if (ierr != ZOLTAN_OK){
    goto End;
  }


End:
  Zoltan_Multifree(__FILE__, __LINE__, 8,
    &import_global_ids, &import_local_ids, &import_procs, &import_to_part,
    &export_global_ids, &export_local_ids, &export_procs, &export_to_part);

  return ierr;
}
/****************************************************************/
/****************************************************************/
/*
 * Called by Zoltan_Destroy()
 */
void Zoltan_MP_Free_Structure(ZZ *zz)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;

  if (mpd != NULL){

    ZOLTAN_FREE(&mpd->rcGID);
    ZOLTAN_FREE(&mpd->nzIndex);
    ZOLTAN_FREE(&mpd->nzGID);

    ZOLTAN_FREE(&mpd->crGID);
    ZOLTAN_FREE(&mpd->mirrorNzIndex);
    ZOLTAN_FREE(&mpd->mirrorNzGID);

    ZOLTAN_FREE(&mpd->vtxGID);

    ZOLTAN_FREE(&mpd->rowproc);
    ZOLTAN_FREE(&mpd->rowpart);
    ZOLTAN_FREE(&mpd->colproc);
    ZOLTAN_FREE(&mpd->colpart);
    ZOLTAN_FREE(&mpd->nzproc);
    ZOLTAN_FREE(&mpd->nzpart);

    free_obj_lookup_table(&mpd->row_lookup);
    free_obj_lookup_table(&mpd->col_lookup);
    free_obj_lookup_table(&mpd->nz_lookup);

    Zoltan_Destroy(&(mpd->zzLib));  /* we created this PHG problem */

    ZOLTAN_FREE(&zz->LB.Data_Structure);
  }
}

/****************************************************************/
/****************************************************************/
/*
 * Functions used by Zoltan_Matrix_Partition()
 */

static int MP_Initialize_Params(ZZ *zz, ZOLTAN_MP_DATA *mpd)
{
  char *yo="MP_Initialize_Params";
  int ierr = ZOLTAN_OK;
  char approach[MAX_PARAM_STRING_LEN];

  Zoltan_Bind_Param(MP_params, "MATRIX_APPROACH", approach);  
  Zoltan_Bind_Param(MP_params, "MATRIX_METHOD", mpd->method);  

  strncpy(approach, "ROWS", MAX_PARAM_STRING_LEN);
  strncpy(mpd->method, "HYPERGRAPH", MAX_PARAM_STRING_LEN);

  ierr = Zoltan_Assign_Param_Vals(zz->Params, MP_params, zz->Debug_Level,
          zz->Proc, zz->Debug_Proc);

  if (ierr == ZOLTAN_OK){
    /* Update this if you add a new MATRIX_APPROACH */
    if ((!strcasecmp(approach, "rows")) ||
        (!strcasecmp(approach, "row")) )
      mpd->approach = MP_ROW_TYPE;
    else if ((!strcasecmp(approach, "columns")) ||
             (!strcasecmp(approach, "cols"))    ||
             (!strcasecmp(approach, "col")) )
      mpd->approach = MP_COLUMN_TYPE;
    else if ((!strcasecmp(approach, "general")) ||
             (!strcasecmp(approach, "nonzeros")) ||
             (!strcasecmp(approach, "nonzero")) )
      mpd->approach = MP_NZ_TYPE;
    else {
    }

    if (zz->Num_GID > sizeof(IJTYPE)){
      /*
       * Maybe they want to use long ints, we only handle ints.  Change
       * IJTYPE to make code use long ints.
       */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input row and column IDs are too large\n");
      ierr = ZOLTAN_FATAL;
    }
  }
  else{
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_Assign_Param_Vals\n");
  }
  return ierr;
}


static int process_matrix_input(ZZ *zz, ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;
  int i;
  long int lnzs, gnzs;
  IJTYPE minID, maxID, minNzID, maxNzID;
  IJTYPE nnz=0;
  long int vals[2], gvals[2];

  /* get global number of rows, columns and nzs, range of row and column IDs */

  if (mpd->numRC > 0)
    nnz = mpd->nzIndex[mpd->numRC];

  lnzs = (long int)nnz;

  MPI_Allreduce(&lnzs, &gnzs, 1, MPI_LONG, MPI_SUM, zz->Communicator);

  mpd->nNonZeros = (IJTYPE)gnzs;

  maxID = maxNzID = -1;

  if (mpd->numRC > 0){
    minID = maxID = mpd->rcGID[0];

    for (i=1; i<mpd->numRC; i++){
      if (mpd->rcGID[i] < minID) minID = mpd->rcGID[i];
      else if (mpd->rcGID[i] > maxID) maxID = mpd->rcGID[i];
    }
    if (nnz > 0){
      minNzID = maxNzID = mpd->nzGID[0];
      for (i=1; i<nnz; i++){
        if (mpd->nzGID[i] < minNzID) minNzID = mpd->nzGID[i];
        else if (mpd->nzGID[i] > maxNzID) maxNzID = mpd->nzGID[i];
      }
    }
  }
  vals[0] = (long int)maxID;
  vals[1] = (long int)maxNzID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MAX, zz->Communicator);

  maxID = (IJTYPE)gvals[0];
  maxNzID = (IJTYPE)gvals[1];

  if (nnz == 0){
    minNzID = maxNzID;
    if (mpd->numRC == 0){
      minID = maxID;
    }
  }
  vals[0] = (long int)minID;
  vals[1] = (long int)minNzID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MIN, zz->Communicator);

  minID = (IJTYPE)gvals[0];
  minNzID = (IJTYPE)gvals[1];

  if (mpd->input_type == ROW_TYPE){
    mpd->rowBaseID = minID;
    mpd->colBaseID = minNzID;
    mpd->nRows = maxID - minID + 1;
    mpd->nCols = maxNzID - minNzID + 1;
  }
  else{
    mpd->rowBaseID = minNzID;
    mpd->colBaseID = minID;
    mpd->nRows = maxNzID - minNzID + 1;
    mpd->nCols = maxID - minID + 1;
  }

  return ierr;
}
/* TODO - need to copy structure too.
 */
static ZOLTAN_MP_DATA *MP_Initialize_Structure(ZZ *zz)
{
  ZOLTAN_MP_DATA *mpd = 
    (ZOLTAN_MP_DATA *)ZOLTAN_CALLOC(sizeof(ZOLTAN_MP_DATA), 1);

  if (mpd == NULL){
    ZOLTAN_PRINT_ERROR(zz->Proc, "MP_Initialize_Structure", "Memory error\n");
    return NULL;
  }

  return mpd;
}
static int make_mirror(ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;
  int nnz, nIDs;

  /*
   * If sparse matrix was supplied in CSC, create another in CSR format,
   * and vice versa.
   */
  ZOLTAN_FREE(&mpd->crGID);
  ZOLTAN_FREE(&mpd->mirrorNzIndex);
  ZOLTAN_FREE(&mpd->mirrorNzGID);

  nIDs = mpd->numRC;
  nnz = mpd->nzIndex[nIDs];

  mpd->numCR = mpd->numRC;
  allocate_copy(&mpd->crGID, mpd->nzGID, nnz);
  allocate_copy(&mpd->mirrorNzIndex, mpd->nzIndex, nIDs+1);
  allocate_copy(&mpd->mirrorNzGID, mpd->rcGID, nIDs);

  /* TODO  Convert_to_CSR thinks some of these fields are ints
   *         when IJTYPEs might not be ints.  FIX, but OK for now.
   *
   */

  ierr = Zoltan_Convert_To_CSR(mpd->zzLib, nnz, (int *)mpd->nzIndex,
             /* The following get overwritten with mirror */
             (int *)&(mpd->numCR),
             (ZOLTAN_ID_PTR *)&(mpd->mirrorNzGID),
             (int **)&(mpd->mirrorNzIndex),
             (ZOLTAN_ID_PTR *)&(mpd->crGID));

  return ierr;
}
static int allocate_copy(IJTYPE **to, IJTYPE *from, IJTYPE len)
{
  *to = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * len);
  memcpy(*to, from, sizeof(IJTYPE) * len);
  return ZOLTAN_OK;
}
/****************************************************************/
/****************************************************************/
/* Functions to support a hash table mapping global IDs to local
 * IDs, which are indices into an array.  The global ID can be
 * a single IJTYPE or an IJTYPE pair.
 */
static void free_obj_lookup_table(obj_lookup **lu)
{
  obj_lookup *l = *lu;

  if (l == NULL) return;

  ZOLTAN_FREE(&l->htTop);
  ZOLTAN_FREE(&l->ht);
  ZOLTAN_FREE(lu);
}
/*
 * Look up array index for global ID or for (i,j) pair.
 * This is a global function, in case other parts of Zoltan
 * library want to use the lookup tables in the ZOLTAN_MP_DATA
 * object that is saved at zz->LB.Data_Structure.
 */
int Zoltan_Lookup_Obj(obj_lookup *lu, IJTYPE I, IJTYPE J)
{
  struct obj_node *hn;
  unsigned int hashVal;
  ZOLTAN_ID_TYPE ids[2];

  ids[0] = I;
  if (lu->key_size > 1) ids[1] = J;

  hashVal = Zoltan_Hash(ids, lu->key_size, lu->table_size);

  for (hn=lu->ht[hashVal]; hn != NULL; hn = hn->next){
    if (hn->i == I){
      if ((lu->key_size == 1) || (hn->j == J)){
        return hn->objLID;
      }
    }
  }
  return -1;
}
/*
** Create a hash table mapping the list of GIDs in objGIDs to their
** array index in the objGIDs list.  GIDs in objGIDs are unique.
*/
static obj_lookup *create_obj_lookup_table(IJTYPE numObjs, IJTYPE *objGIDs)
{
  IJTYPE i;
  unsigned int hashVal;
  struct obj_node *hn;
  obj_lookup *lu = NULL;
  int size_id;

  lu = (obj_lookup *)ZOLTAN_MALLOC(sizeof(obj_lookup));
  if (!lu){
    return NULL;
  }

  lu->key_size = 1;
  lu->table_size = numObjs / 4;
  if (lu->table_size < 4) lu->table_size = 4;

  lu->ht = 
  (struct obj_node **)ZOLTAN_CALLOC(sizeof(struct obj_node*) , lu->table_size);

  hn = lu->htTop = 
  (struct obj_node *)ZOLTAN_MALLOC(sizeof(struct obj_node) * numObjs);

  if (lu->table_size && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu);
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    return NULL;
  }

  size_id = sizeof(IJTYPE) / sizeof(ZOLTAN_ID_TYPE);
  if (size_id < 1) size_id = 1;

  for (i=0; i<numObjs; i++){

    hn->i = objGIDs[i];
    hn->j = 0;
    hn->objLID = i;

    hashVal = Zoltan_Hash((ZOLTAN_ID_PTR)&objGIDs[i], size_id, lu->table_size);

    hn->next = lu->ht[hashVal];
    lu->ht[hashVal] = hn;

    hn++;
  }

  return lu;
}
/*
** Create a hash table mapping the (I,J) pairs to their location in
** listJ.  IDs in listI are unique.  Number of pairs is sizeJ.
*/
static obj_lookup *create_obj_lookup_table2(IJTYPE sizeI, IJTYPE sizeJ,
                   IJTYPE *listI, IJTYPE *indexJ, IJTYPE *listJ)
{
  IJTYPE i, j;
  unsigned int hashVal;
  struct obj_node *hn;
  obj_lookup *lu = NULL;
  ZOLTAN_ID_TYPE ids[2];

  if (sizeof(IJTYPE) > sizeof(ZOLTAN_ID_TYPE)){
    ZOLTAN_PRINT_ERROR(0, "create_obj_lookup_table2", 
                       "code rewrite required.\n");
    exit(1);
  }

  lu = (obj_lookup *)ZOLTAN_MALLOC(sizeof(obj_lookup));
  if (!lu){
    return NULL;
  }

  lu->key_size = 2;
  lu->table_size = sizeJ / 4;
  if (lu->table_size < 4) lu->table_size = 4;

  lu->ht = 
  (struct obj_node **)ZOLTAN_CALLOC(sizeof(struct obj_node*) , lu->table_size);

  hn = lu->htTop = 
  (struct obj_node *)ZOLTAN_MALLOC(sizeof(struct obj_node) * sizeJ);

  if (lu->table_size && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu);
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    return NULL;
  }

  for (i=0; i<sizeI; i++){
    for (j = indexJ[i]; j < indexJ[i+1]; j++){
      hn->i = listI[i];
      hn->j = listJ[j];
      hn->objLID = j;

      ids[0] = listI[i];
      ids[1] = listJ[j];

      hashVal = Zoltan_Hash(ids, 2, lu->table_size);

      hn->next = lu->ht[hashVal];
      lu->ht[hashVal] = hn;

      hn++;
    }
  }
  return lu;
}

/****************************************************************/
/****************************************************************/
/* Functions that support MATRIX_APPROACH = MP_ROWS or MP_COLS
 *
 * The objects to be partitioned are the matrix rows or columns,
 * and we are using Zoltan's PHG to do it.
 *
 * MP_ROWS: the objects are the sparse matrix rows, and the
 *  hyperedges are the columns of the sparse matrix.
 *
 * MP_COLS: the objects are the sparse matrix columns, and
 *  the hyperedges are the rows.
 */

ZOLTAN_NUM_OBJ_FN mp_1d_get_num_obj;
ZOLTAN_OBJ_LIST_FN mp_1d_get_obj_list;
ZOLTAN_HG_SIZE_CS_FN mp_1d_get_size_matrix;
ZOLTAN_HG_CS_FN mp_1d_get_matrix;
static int get_proc_part(ZOLTAN_MP_DATA *mpd, int num_export,
  unsigned int *gids, int *procs, int *parts,
  IJTYPE nobj, IJTYPE *objList, int *objProcs, int *objParts);
static int mp_1d_obj_to_proc(ZOLTAN_MP_DATA *mpd, IJTYPE objID);
static int mp_1d_my_objects(ZOLTAN_MP_DATA *mpd, IJTYPE *nobj, IJTYPE **objIDs);
/*
 *
 */
static int mp_1d_get_nzs(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nNzs, IJTYPE *I, IJTYPE *J, int *ijProcs, int *ijParts)
{
  if (mpd->approach == MP_ROW_TYPE){
    /* The nz belongs to the partition its row belongs to */
    return mp_1d_get_rows(mpd, nNzs, I, ijProcs, ijParts);
  }
  else{
    /* The nz belongs to the partition its column belongs to */
    return mp_1d_get_columns(mpd, nNzs, J, ijProcs, ijParts);
  }
}
static int mp_1d_get_columns(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nCols, IJTYPE *colIDs, int *colProcs, int *colParts)
{
  int idx, i;
  char *yo = "mp_1d_get_columns";
  int nTotalCols = ((mpd->input_type == ROW_TYPE) ? mpd->numCR : mpd->numRC);

  if (mpd->approach != MP_COLUMN_TYPE){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, 
     "chosen MATRIX_APPROACH does not support obtaining column partitions\n");
    return ZOLTAN_FATAL; 
  }

  for (i=0; i<nCols; i++){ 
    idx = Zoltan_Lookup_Obj(mpd->col_lookup, colIDs[i], 0);
    if ((idx < 0) || (idx >= nTotalCols)){
      ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Unable to determine column partition\n");
     return ZOLTAN_FATAL; 
    }
    if (colProcs) colProcs[i] = mpd->colproc[idx];
    if (colParts) colParts[i] = mpd->colpart[idx];
  }
  return ZOLTAN_OK;
}
static int mp_1d_get_rows(ZOLTAN_MP_DATA *mpd, 
        IJTYPE nRows, IJTYPE *rowIDs, int *rowProcs, int *rowParts)
{
  int idx, i;
  char *yo = "mp_1d_get_rows";
  int nTotalRows = ((mpd->input_type == ROW_TYPE) ? mpd->numRC : mpd->numCR);

  if (mpd->approach != MP_ROW_TYPE){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, 
     "chosen MATRIX_APPROACH does not support obtaining row partitions\n");
    return ZOLTAN_FATAL; 
  }

  for (i=0; i<nRows; i++){ 
    idx = Zoltan_Lookup_Obj(mpd->row_lookup, rowIDs[i], 0);
    if ((idx < 0) || (idx >= nTotalRows)){
      ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Unable to determine row partition\n");
     return ZOLTAN_FATAL; 
    }
    if (rowProcs) rowProcs[i] = mpd->rowproc[idx];
    if (rowParts) rowParts[i] = mpd->rowpart[idx];
  }
  return ZOLTAN_OK;
}
static int mp_1d_obj_to_proc(ZOLTAN_MP_DATA *mpd, IJTYPE objID)
{
  int nProcs = mpd->zzLib->Num_Proc;

  return objID % nProcs;
}
static int mp_1d_my_objects(ZOLTAN_MP_DATA *mpd, IJTYPE *nobj, IJTYPE **objIDs)
{
  IJTYPE i, nmyids=0, *obj=NULL;
  IJTYPE baseid, lastid;
  int ierr = ZOLTAN_OK;
  int me = mpd->zzLib->Proc;

  if (mpd->approach == MP_ROW_TYPE){
    baseid = mpd->rowBaseID;
    lastid = mpd->rowBaseID +  mpd->nRows - 1;
  }
  else{
    baseid = mpd->colBaseID;
    lastid = mpd->colBaseID +  mpd->nCols - 1;
  }

  for (i=baseid; i<=lastid; i++){
    if (mp_1d_obj_to_proc(mpd, i) == me){
      nmyids++;
    }
  }

  obj = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * nmyids);
  if (nmyids && !obj){
    ierr = ZOLTAN_MEMERR;
  }
  else{
    *objIDs = obj;
    *nobj = nmyids;

    for (i=baseid; i<=lastid; i++){
      if (mp_1d_obj_to_proc(mpd, i) == me){
        *obj++ = i;
      }
    }
  }

  return ierr;
}

static int mp_1d_setup(ZOLTAN_MP_DATA *mpd)
{
  char *yo = "mp_1d_setup";
  int ierr = ZOLTAN_OK;

  /* TODO Copy parameters from parent ZZ instance. */
  /* Set specific parmeters for Zoltan subproblem. */
  Zoltan_Set_Param(mpd->zzLib, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(mpd->zzLib, "NUM_LID_ENTRIES", "1");

  Zoltan_Set_Param(mpd->zzLib, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(mpd->zzLib, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(mpd->zzLib, "ADD_OBJ_WEIGHT", "PINS");

  /* Request export list that has proc/part for every object */
  Zoltan_Set_Param(mpd->zzLib, "RETURN_LISTS", "PARTITION");

  /* Set LB_METHOD from MP_METHOD, default is HYPERGRAPH. */
  /* Zoltan_Set_Param(mpd->zzLib, "LB_METHOD", "HYPERGRAPH"); */
  Zoltan_Set_Param(mpd->zzLib, "LB_METHOD", mpd->method);

  Zoltan_Set_Num_Obj_Fn(mpd->zzLib, mp_1d_get_num_obj, mpd);
  Zoltan_Set_Obj_List_Fn(mpd->zzLib, mp_1d_get_obj_list, mpd);
  Zoltan_Set_HG_Size_CS_Fn(mpd->zzLib, mp_1d_get_size_matrix, mpd);
  Zoltan_Set_HG_CS_Fn(mpd->zzLib, mp_1d_get_matrix, mpd);

  ierr = make_mirror(mpd);

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "make_mirror error\n");
  }
  else{
    ierr = mp_1d_my_objects(mpd, &mpd->nMyVtx, &mpd->vtxGID);
  }

  return ierr;
}

int mp_1d_get_num_obj(void *data, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;

  return mpd->nMyVtx;
}

void mp_1d_get_obj_list(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts, int *ierr)
{
  IJTYPE i;
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;
 
  for (i=0; i<mpd->nMyVtx; i++){
    gids[i] = (ZOLTAN_ID_TYPE)mpd->vtxGID[i];
    lids[i] = (ZOLTAN_ID_TYPE)i;
  }
}

void mp_1d_get_size_matrix(void *data, int *num_lists, int *num_nzs,
  int *format, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;

  if (((mpd->approach == MP_ROW_TYPE) && (mpd->input_type == ROW_TYPE)) ||
      ((mpd->approach == MP_COLUMN_TYPE) && (mpd->input_type == COL_TYPE)) ){
    *num_lists = mpd->numCR;
    *num_nzs  = mpd->mirrorNzIndex[mpd->numCR];
  }
  else{
    *num_lists = mpd->numRC;
    *num_nzs  = mpd->nzIndex[mpd->numRC];
  }
  *format = ZOLTAN_COMPRESSED_EDGE;
}

void mp_1d_get_matrix(void *data, int num_gid_entries, int num_vtx_edge,
  int num_nzs, int format, 
  ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, ZOLTAN_ID_PTR nz_GID, int *ierr)
{
  int i;
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;
  IJTYPE *edgeGID = (IJTYPE *)vtxedge_GID;
  IJTYPE *nzGID = (IJTYPE *)nz_GID;

  if (((mpd->approach == MP_ROW_TYPE) && (mpd->input_type == ROW_TYPE)) ||
      ((mpd->approach == MP_COLUMN_TYPE) && (mpd->input_type == COL_TYPE)) ){
    for (i=0; i<mpd->numCR; i++){
      edgeGID[i] = mpd->crGID[i];
      vtxedge_ptr[i] = mpd->mirrorNzIndex[i];
    }
    for (i=0; i< mpd->mirrorNzIndex[mpd->numCR]; i++){
      nzGID[i] = mpd->mirrorNzGID[i];
    }
  }
  else{
    for (i=0; i<mpd->numRC; i++){
      edgeGID[i] = mpd->rcGID[i];
      vtxedge_ptr[i] = mpd->nzIndex[i];
    }
    for (i=0; i< mpd->nzIndex[mpd->numRC]; i++){
      nzGID[i] = mpd->nzGID[i];
    }
  }
}

static int mp_1d_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  int *export_procs, int *export_to_part)
{
  char *yo = "mp_1d_result";
  int ierr = ZOLTAN_OK;
  int *proclist=NULL, *partlist=NULL;
  obj_lookup *lu=NULL;
  IJTYPE nobj=0;
  IJTYPE *objList=NULL;

  if (mpd->approach == MP_ROW_TYPE){
    if (mpd->input_type == ROW_TYPE){
      nobj = mpd->numRC;
      objList = mpd->rcGID;
    }
    else{
      nobj = mpd->numCR;
      objList = mpd->crGID;
    }
    lu = mpd->row_lookup = create_obj_lookup_table(nobj, objList);
    proclist = mpd->rowproc = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
    partlist = mpd->rowpart = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
  }
  else{
    if (mpd->input_type == ROW_TYPE){
      nobj = mpd->numCR;
      objList = mpd->crGID;
    }
    else{
      nobj = mpd->numRC;
      objList = mpd->rcGID;
    }
    lu = mpd->col_lookup = create_obj_lookup_table(nobj, objList);
    proclist = mpd->colproc = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
    partlist = mpd->colpart = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
  }

  if (nobj && (!lu || !proclist || !partlist)){
    free_obj_lookup_table(&lu);
    ZOLTAN_FREE(&proclist);
    ZOLTAN_FREE(&partlist);
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Out of memory.\n");
    ierr = ZOLTAN_MEMERR;
  }
  else{
    ierr = get_proc_part(mpd, num_export, 
         export_global_ids, export_procs, export_to_part, 
         nobj, objList, proclist, partlist);
  }

  return ierr;
}

static int get_proc_part(ZOLTAN_MP_DATA *mpd, int num_export,
  unsigned int *gids, int *procs, int *parts,
  IJTYPE nobj, IJTYPE *objList, int *objProcs, int *objParts)
{
  char *yo = "get_proc_part";
  obj_lookup *myIDs=NULL;
  ZOLTAN_COMM_OBJ *plan=NULL;
  int *toProc=NULL, *outprocs=NULL, *outparts=NULL;
  IJTYPE *idReqs=NULL;
  int tag = 35000, nrecv=0, idx, i;
  MPI_Comm comm = mpd->zzLib->Communicator;
  int ierr = ZOLTAN_OK;

  /*
   * Create lookup table for my gids so I can reply to other procs.
   */
  myIDs = create_obj_lookup_table(num_export, gids);
  if (myIDs == NULL){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "create_obj_lookup_table\n");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /*
   * Create a communication plan, requesting new proc/part for
   * all IDs in my objList.
   */
  toProc = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
  if (nobj && !toProc){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "memory\n");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  for (i=0; i<nobj; i++){
    toProc[i] = mp_1d_obj_to_proc(mpd, objList[i]);
  }

  ierr = Zoltan_Comm_Create(&plan, nobj, toProc, comm, tag, &nrecv);
  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "comm create\n");
    goto End;
  }
  /*
   * Send list of obj IDs I need proc/part for, get list from others.
   */
  idReqs = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * nrecv);
  outprocs = (int *)ZOLTAN_MALLOC(sizeof(int) * nrecv);
  outparts = (int *)ZOLTAN_MALLOC(sizeof(int) * nrecv);
  if (nrecv && (!idReqs || !outprocs || !outparts)){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "memory\n");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  ierr = Zoltan_Comm_Do(plan, --tag, (char *)objList, sizeof(IJTYPE),
                        (char *)idReqs);
  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "comm do\n");
    goto End;
  }

  /*
   * Create list of proc/part for all requests I received.
   */
  for (i=0; i<nrecv; i++){
    idx = Zoltan_Lookup_Obj(myIDs, idReqs[i], 0);
    if (idx == -1){
      ierr = ZOLTAN_FATAL;
      ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Zoltan_Lookup_Obj\n");
      goto End;
    }
    outprocs[i] = procs[idx];
    outparts[i] = parts[idx];
  }

  /*
   * Send proc/parts to others and also get my requests.
   */
  ierr = Zoltan_Comm_Do_Reverse(plan, --tag, (char *)outprocs,
             sizeof(unsigned int), NULL, (char *)objProcs);
  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "comm do reverse\n");
    goto End;
  }

  ierr = Zoltan_Comm_Do_Reverse(plan, --tag, (char *)outparts,
             sizeof(unsigned int), NULL, (char *)objParts);
  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "comm do reverse 2\n");
    goto End;
  }

End:
  Zoltan_Comm_Destroy(&plan); 
  free_obj_lookup_table(&myIDs);
  ZOLTAN_FREE(&toProc);
  ZOLTAN_FREE(&idReqs);
  ZOLTAN_FREE(&outprocs);
  ZOLTAN_FREE(&outparts);

  return ierr;
}

/****************************************************************/
/****************************************************************/
/* Functions that support MATRIX_APPROACH = MP_NONZEROS
 *
 * The objects to be partitioned are the matrix nonzeros.
 *
 */

ZOLTAN_NUM_OBJ_FN mp_2d_get_num_obj;
ZOLTAN_OBJ_LIST_FN mp_2d_get_obj_list;
ZOLTAN_NUM_GEOM_FN mp_2d_get_num_geom;
ZOLTAN_GEOM_MULTI_FN mp_2d_get_geom;
/* ZOLTAN_HG_SIZE_CS_FN mp_2d_get_size_matrix; */
/* ZOLTAN_HG_CS_FN mp_2d_get_matrix; */

static int mp_2d_setup(ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;

  Zoltan_Set_Param(mpd->zzLib, "NUM_GID_ENTRIES", "2");
  Zoltan_Set_Param(mpd->zzLib, "NUM_LID_ENTRIES", "2");

  Zoltan_Set_Param(mpd->zzLib, "LB_APPROACH", "PARTITION");

  /* Request export list that has proc/part for every object */
  Zoltan_Set_Param(mpd->zzLib, "RETURN_LISTS", "PARTITION");

  /* Set LB_METHOD from MP_METHOD, default is HYPERGRAPH. */
  if (strcmp(mpd->method, "HYPERGRAPH")==0) /* TODO */
    strcpy(mpd->method, "BLOCK"); 
  Zoltan_Set_Param(mpd->zzLib, "LB_METHOD", mpd->method);

  Zoltan_Set_Num_Obj_Fn(mpd->zzLib, mp_2d_get_num_obj, mpd);
  Zoltan_Set_Obj_List_Fn(mpd->zzLib, mp_2d_get_obj_list, mpd);
  /* Zoltan_Set_HG_Size_CS_Fn(mpd->zzLib, mp_2d_get_size_matrix, mpd); */
  /* Zoltan_Set_HG_CS_Fn(mpd->zzLib, mp_2d_get_matrix, mpd); */

/*
  ierr = make_mirror(mpd);

  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "make_mirror error\n");
  }
  else{
    ierr = mp_1d_my_objects(mpd, &mpd->nMyVtx, &mpd->vtxGID);
  }
*/

  return ierr;
}

/* Basic 2d query functions. Each nonzero is a Zoltan object. */
int mp_2d_get_num_obj(void *data, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;

  return mpd->nzIndex[mpd->numRC];
}

void mp_2d_get_obj_list(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts, int *ierr)
{
  IJTYPE i, j;
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  char *yo = "mp_2d_get_obj_list";
 
  if (mpd) 
    *ierr = ZOLTAN_OK;
  else
    *ierr = ZOLTAN_FATAL;

  if (!(num_gid_entries==2)){
    ZOLTAN_PRINT_ERROR(0, yo, "Assumed num_gid_entries==2");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  /* Loop over nonzeros. Each ID is an (i,j) pair. */
  /* Note: Row major, i.e. outer loop = rows (i), inner loop = columns (j). */
  for (i=0; i<mpd->numRC; i++){
    lids[2*i] = (ZOLTAN_ID_TYPE)i;
    gids[2*i] = (ZOLTAN_ID_TYPE)mpd->rcGID[i];
    for (j=mpd->nzIndex[i]; j<mpd->nzIndex[i+1]; j++){
      lids[2*i+1] = (ZOLTAN_ID_TYPE)j;
      gids[2*i+1] = (ZOLTAN_ID_TYPE)mpd->nzGID[j];
    }
  }
}

/* Geometric 2d query functions. Each nonzero has coordinates (i,j) */
int mp_2d_get_num_geom(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;

  return 2;
}

void mp_2d_get_geom(
       void *data, int num_gid_entries, int num_lid_entries, 
       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, 
       int num_dim, double *geom_vec, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  int i;
  char *yo = "mp_2d_get_geom";

  if (mpd) 
    *ierr = ZOLTAN_OK;
  else
    *ierr = ZOLTAN_FATAL;

  if (!(num_gid_entries==2)){
    ZOLTAN_PRINT_ERROR(0, yo, "Assumed num_gid_entries==2");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  /* Map matrix so A(1,1) is at (0,0) and is top left. */
  for (i=0; i<num_obj; i++){
    geom_vec[2*i]   = (double) global_ids[2*i+1] - 1.;
    geom_vec[2*i+1] = (double) -(global_ids[2*i]-1.);
  }

}


/* TODO: Hypergraph query functions. These should build fine-grain hgraph. */

/* Store partition results in mpd struct. */
static int mp_2d_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  int *export_procs, int *export_to_part)
{
  char *yo = "mp_2d_result";
  int ierr = ZOLTAN_OK;
  int *proclist=NULL, *partlist=NULL;
  obj_lookup *lu=NULL;
  IJTYPE nobj=0;
  IJTYPE i, j;

  nobj = mpd->nMyVtx;
  // TODO What are i and j and arrays for create_obj_lookup_table2?
  lu = mpd->nz_lookup = create_obj_lookup_table2(i, j, NULL, NULL, NULL);
  proclist = mpd->nzproc = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);
  partlist = mpd->nzpart = (int *)ZOLTAN_MALLOC(sizeof(int) * nobj);

  if (nobj && (!lu || !proclist || !partlist)){
    free_obj_lookup_table(&lu);
    ZOLTAN_FREE(&proclist);
    ZOLTAN_FREE(&partlist);
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Out of memory.\n");
    ierr = ZOLTAN_MEMERR;
  }
  else{
    /* Copy partition results to proclist and partlist. */
    memcpy(proclist, export_procs, nobj*sizeof(int));
    memcpy(partlist, export_to_part, nobj*sizeof(int));
  }

  return ierr;
}

/****************************************************************/
/****************************************************************/

/* For small sparse matrices used in testing, process 0 prints
 * out the matrix with each nz value equal to the nz's partition.
 * Matrix should be small enough that we can print it out.
 */
void Zoltan_MP_Debug_Partitioning(ZZ *zz)
{
  char *yo = "Zoltan_MP_Debug_Partitioning";
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)zz->LB.Data_Structure;
  IJTYPE *colIDs, *rowIDs, *nzIDs, *nzIdx, *IDs, *counts=NULL, *recvIDs2=NULL;
  IJTYPE numRows, numCols, numIDs, row, col, numNzs, numIJs, totNumNzs, baseID, idx;
  IJTYPE totNumIJs=0, *recvIDs=NULL, *nzs=NULL, *Ival=NULL;
  int rc, me, nprocs, i, j, nextIdx, oops, numParts, part2D;
  int totColCuts, totRowCuts, width, totCount;
  MPI_Comm comm = zz->Communicator;
  int *parts=NULL, *recvParts=NULL, *recvDisp=NULL, *recvCounts=NULL;
  int **A=NULL;
  int *colCuts=NULL, *rowCuts=NULL, *countParts=NULL;

  if ((mpd->nRows > 100) || (mpd->nCols > 100)){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Matrix is too large.\n");
    exit(0);
  }

  if (sizeof(IJTYPE) != sizeof(unsigned int)){
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, 
      "REWRITE required due to change in IJTYPE\n");
    exit(0);
  }

  me = mpd->zzLib->Proc;
  nprocs = mpd->zzLib->Num_Proc;

  /* Find out where the non-zeroes are */

  if (me == 0){
    A = (int **)ZOLTAN_MALLOC(sizeof(int *) * mpd->nRows);
    for (i=0; i<mpd->nRows; i++){
      A[i] = (int *)ZOLTAN_MALLOC(sizeof(int) * mpd->nCols);
      for (j=0; j < mpd->nCols; j++){
        A[i][j] = -2;       /* not a nz */
      }
    }
  }

  if (mpd->input_type == ROW_TYPE){
    rowIDs = mpd->rcGID;
    numRows = mpd->numRC;
    colIDs = mpd->crGID;
    numCols = mpd->numCR;
    nzIDs = mpd->nzGID;            /* GID of nz column */
    nzIdx = mpd->nzIndex;
  }
  else{
    rowIDs = mpd->crGID;
    numRows = mpd->numCR;
    colIDs = mpd->rcGID;
    numCols = mpd->numRC;
    nzIDs = mpd->mirrorNzGID;      /* GID of nz column */
    nzIdx = mpd->mirrorNzIndex;
  }

  numNzs = nzIdx[numRows];
  numIJs = numNzs * 2;

  nzs = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * numIJs);
  for (i=0,nextIdx=0; i<numRows; i++){
    row = rowIDs[i] - mpd->rowBaseID;
    for (j=nzIdx[i]; j<nzIdx[i+1]; j++){
      col = nzIDs[j] - mpd->colBaseID;
      nzs[nextIdx++] = row;
      nzs[nextIdx++] = col;
    }
  }

  if (me == 0){
    recvCounts = (int *)ZOLTAN_MALLOC(sizeof(int) * nprocs);
    recvDisp = (int *)ZOLTAN_MALLOC(sizeof(int) * nprocs);
  }

  MPI_Gather(&numIJs, 1, MPI_UNSIGNED, recvCounts, 1, MPI_INT, 0, comm);

  if (me == 0){
    totNumIJs = 0;
    for (i=0; i<nprocs; i++){
      recvDisp[i] = totNumIJs;
      totNumIJs += recvCounts[i]; 
    }
    recvIDs = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * totNumIJs);
  }

  MPI_Gatherv(nzs, numIJs, MPI_UNSIGNED,
              recvIDs, recvCounts, recvDisp, MPI_UNSIGNED, 0, comm);

  if (me == 0){
    for (i=0,nextIdx=0; i<nprocs; i++){
      for (j=0; j<recvCounts[i]; j+=2){
        row = recvIDs[recvDisp[i]+j];
        col = recvIDs[recvDisp[i]+j+1];
        if ( (row < 0) || (row > mpd->nRows) ||
             (col < 0) || (col > mpd->nCols)){
          ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Bad ID\n");
          exit(0);
        }
        A[row][col] = -1;     /* it's a nz */
      }   
    }
    totNumNzs = totNumIJs/2;
  }

  ZOLTAN_FREE(&recvIDs);
  ZOLTAN_FREE(&recvCounts);
  ZOLTAN_FREE(&recvDisp);
  ZOLTAN_FREE(&nzs);

  /* Get the partition assignments of the rows, columns or non-zeroes. */

  if ((mpd->row_lookup) || (mpd->col_lookup)){
    part2D = 0;                 /* 1 dimensional partitioning */
  }
  else if (mpd->nz_lookup){
    part2D = 1;                 /* 2 dimensional partitioning */
  }
  else{
    ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "Partitioning not saved\n");
    return;
  }

  if (part2D){
    parts = (int *)ZOLTAN_MALLOC(sizeof(int) * numNzs);
    Ival = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * numNzs);
    for (i=0; i<numRows; i++){
      for (j=0; j<(nzIdx[i+1]-nzIdx[i]); j++){
        Ival[nzIdx[i]+j] = rowIDs[i];
      }
    }

    rc = Zoltan_MP_Get_NonZero_Assignment(zz, numNzs,
                   Ival, nzIDs,    /* I, J coords of my non-zeroes */
                   NULL, parts);

    IDs = Ival;       /* will also send nzIDs */
    numIDs = numNzs;
  }
  else{
    if (mpd->row_lookup){
      parts = (int *)ZOLTAN_MALLOC(sizeof(int) * numRows);
      rc = Zoltan_MP_Get_Row_Assignment(zz, numRows, rowIDs, NULL, parts);
      IDs = rowIDs;
      numIDs = numRows;
    }
    else if (mpd->col_lookup){
      parts = (int *)ZOLTAN_MALLOC(sizeof(int) * numCols);
      rc = Zoltan_MP_Get_Column_Assignment(zz, numCols, colIDs, NULL, parts);
      IDs = colIDs;
      numIDs = numCols;
    }
  }

  if (me == 0){
    recvCounts = (int *)ZOLTAN_MALLOC(sizeof(int) * nprocs);
    recvDisp = (int *)ZOLTAN_MALLOC(sizeof(int) * nprocs);
  }

  MPI_Gather(&numIDs, 1, MPI_UNSIGNED, recvCounts, 1, MPI_INT, 0, comm);

  if (me == 0){
    totCount = 0;
    for (i=0; i<nprocs; i++){
      recvDisp[i] = totCount;
      totCount += recvCounts[i];
    }
    recvIDs = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * totCount);
    recvParts = (int *)ZOLTAN_MALLOC(sizeof(int) * totCount);

    if (part2D){
      recvIDs2 = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * totCount);
    }
  }

  MPI_Gatherv(IDs, numIDs, MPI_UNSIGNED, 
              recvIDs, recvCounts, recvDisp, MPI_UNSIGNED, 0, comm);

  if (part2D){
    MPI_Gatherv(nzIDs, numIDs, MPI_UNSIGNED, 
              recvIDs2, recvCounts, recvDisp, MPI_UNSIGNED, 0, comm);
  }

  MPI_Gatherv(parts, numIDs, MPI_INT, 
              recvParts, recvCounts, recvDisp, MPI_INT, 0, comm);

  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&Ival);
  ZOLTAN_FREE(&counts);

  if (me == 0){
    numParts = 0;

    if (part2D){
      for (nextIdx=0; nextIdx<totCount; nextIdx++){
        i = recvIDs[nextIdx] - mpd->rowBaseID; 
        j = recvIDs2[nextIdx] - mpd->colBaseID; 
        A[i][j] = recvParts[nextIdx];
        if (recvParts[nextIdx] > numParts) numParts = recvParts[nextIdx];
      }
    }
    else{
      if (mpd->row_lookup){
        numIDs = mpd->nRows;
        baseID = mpd->rowBaseID;
      }
      else{
        numIDs = mpd->nCols;
        baseID = mpd->colBaseID;
      }
    
      parts = (int *)ZOLTAN_MALLOC(sizeof(int) * numIDs);
    
      for (i=0; i<numIDs; i++){
        parts[i] = -1;
      }
    
      for (i=0; i<nprocs; i++){
        for (j=recvDisp[i]; j < recvDisp[i] + recvCounts[i]; j++){
          idx = recvIDs[j] - baseID;
    
          if ((idx < 0) || (idx >= numIDs)){
            ZOLTAN_PRINT_ERROR(mpd->zzLib->Proc, yo, "bad id\n");
            exit(0);
          }
          if (parts[idx] < 0){
            parts[idx] = recvParts[j];
            if (parts[idx] > numParts) numParts = parts[idx];
            if (mpd->row_lookup){
              for (nextIdx=0; nextIdx<mpd->nCols; nextIdx++){
                if (A[idx][nextIdx] == -1) A[idx][nextIdx] = parts[idx];
              }
            }
            else{
              for (nextIdx=0; nextIdx<mpd->nRows; nextIdx++){
                if (A[nextIdx][idx] == -1) A[nextIdx][idx] = parts[idx];
              }
            }
          }
        }
      }
      ZOLTAN_FREE(&parts);
    }
    numParts++;
    ZOLTAN_FREE(&recvIDs);
    ZOLTAN_FREE(&recvIDs2);
    ZOLTAN_FREE(&recvParts);
  }

  /* Process zero can print out the matrix */

  fflush(stdout);
  MPI_Barrier(comm);
  if (me == 0){
    /* count the cuts in rows and columns */
    oops = 0;
    colCuts = (int *)ZOLTAN_CALLOC(sizeof(int) , mpd->nCols);
    rowCuts = (int *)ZOLTAN_CALLOC(sizeof(int) , mpd->nRows);
    countParts = (int *)ZOLTAN_CALLOC(sizeof(int), numParts);
    totColCuts = totRowCuts = 0;
    for (i=0; i<mpd->nRows; i++){
      for (j=0; j<mpd->nCols; j++){
        if (A[i][j] >= 0) countParts[A[i][j]] = i+1;
      }
      rowCuts[i] = -1;
      for (j=0; j<numParts; j++){
        if (countParts[j] > i) rowCuts[i]++;
      }
      if (rowCuts[i] > 0){
        totRowCuts += rowCuts[i];
      }
      else{
        rowCuts[i] = 0;
      }
    }
    memset(countParts, 0, sizeof(int) * numParts);
    for (j=0; j<mpd->nCols; j++){
      for (i=0; i<mpd->nRows; i++){
        if (A[i][j] >= 0) countParts[A[i][j]] = j+1;
      }
      colCuts[j] = -1;
      for (i=0; i<numParts; i++){
        if (countParts[i] > j) colCuts[j]++;
      }
      if (colCuts[j] > 0){
        totColCuts += colCuts[j];
      }
      else{
        colCuts[j] = 0;
      }
    }
    ZOLTAN_FREE(&countParts);

    /* print out the sparse matrix */

    if (numParts > 10) width = 3;
    else               width = 2;

    printf("\n     ");
    for (j=0; j<mpd->nCols; j++){
      if (width==3)
        printf("%3d",(j + mpd->colBaseID) % 100);
      else
        printf("%2d",(j + mpd->colBaseID) % 10);
    }
    printf("| cuts\n");
    j = 5 + (mpd->nCols * ((width==3) ? 3 : 2));
    for (i=0; i<j; i++) printf("=");
    printf("|\n");

    for (i=0; i<mpd->nRows; i++){
      printf("%3d |",i + mpd->rowBaseID);
      if (width==3){
        for (j=0; j<mpd->nCols; j++){
          if (A[i][j] == -2) printf(" - "); 
          else if (A[i][j] == -1){
            printf(" e ");
            oops = 1;
          }
          else{
            printf("%3d",A[i][j]);
          }
        }
        printf("|%3d\n",rowCuts[i]);
      }
      else{
        for (j=0; j<mpd->nCols; j++){
          if (A[i][j] == -2) printf(" -"); 
          else if (A[i][j] == -1){
            printf(" e");
            oops = 1;
          }
          else{
            printf("%2d",A[i][j]);
          }
        }
        printf("|%2d\n",rowCuts[i]);
      }
    }
    j = 5 + (mpd->nCols * ((width==3) ? 3 : 2));
    for (i=0; i<j; i++) printf("=");
    printf("|\n cuts");
    for (j=0; j<mpd->nCols; j++){
      if (width == 3)
        printf("%3d",colCuts[j]);
      else
        printf("%2d",colCuts[j]);
    }
    printf("\n\n"); 

    ZOLTAN_FREE(&rowCuts);
    ZOLTAN_FREE(&colCuts);

    printf("Row cuts: total %d average %f\n",
           totRowCuts,(double)totRowCuts/mpd->nRows);
    printf("Col cuts: total %d average %f\n",
           totColCuts,(double)totColCuts/mpd->nCols);
    printf("Number of non-zeroes: %d\n",totNumNzs);
  
    if (oops){
      printf("An \"e\" means an error - we don't know the partition\n");
    }
    printf("\n");
  }
  fflush(stdout);
  MPI_Barrier(comm);
  fflush(stdout);
  MPI_Barrier(comm);
  
/*End:*/
  ZOLTAN_FREE(&nzs);
  ZOLTAN_FREE(&recvCounts);
  ZOLTAN_FREE(&recvDisp);
  ZOLTAN_FREE(&recvIDs);
  ZOLTAN_FREE(&recvParts);
  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&counts);

  if (A != NULL){
    for (i=0; i<mpd->nRows; i++) {ZOLTAN_FREE(&A[i]);}
    ZOLTAN_FREE(&A);
  }
}

/****************************************************************/
/****************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_MATRIX_PARTITION */
