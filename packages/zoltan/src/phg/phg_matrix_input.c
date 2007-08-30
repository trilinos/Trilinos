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

/* TODO: check all the time whether mpd->numRC is zero before proceeding */

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*
 * Assumption is we are creating a hypergraph and using Zoltan_PHG
 * to partition objects, but eventually this code should come out
 * of the phg directory and be more general, using other partitioning
 * methods on the sparse matrix.
 */
#include "phg.h"
/*******************/
#include "zz_util_const.h"

ZOLTAN_NUM_OBJ_FN Zoltan_MP_Get_Num_Obj;
ZOLTAN_OBJ_LIST_FN Zoltan_MP_Get_Obj_List;
ZOLTAN_HG_SIZE_CS_FN Zoltan_MP_Get_Size_Matrix;
ZOLTAN_HG_CS_FN Zoltan_MP_Get_Matrix;

static PARAM_VARS MP_params[] = {
  {"LB_APPROACH",         NULL,  "STRING",     0},
  {"NUM_GID_ENTRIES",     NULL,  "INT",     0},
  {NULL,                  NULL,  NULL,     0}
};
    
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
** This source file uses the terms "pins" and "non-zeros" interchangeably.
**
** An IJTYPE is the type of a row or column ID.
**
** LB_METHOD = SPARSE_MATRIX
** LB_APPROACH:
**    PHG_ROWS    use Zoltan's PHG to partition the rows (objects are
**                rows, hyperedges are columns)
**                
**    PHG_COLS    use Zoltan's PHG to partition the columns (objects
**                are columns, hyperedges are rows)
**
** To find the code that must be added if you add a new approach,
** search for "LB_APPROACH" in this source file.
*************************************************************************/
struct vtx_node{
  IJTYPE vtxGID;
  int    vtxLID;
  struct vtx_node *next;
};
typedef struct _vtx_lookup{
  struct vtx_node *htTop;
  struct vtx_node **ht;
  IJTYPE table_size;
  int    numVtx;
  int    idLength;
}vtx_lookup;

#if 0
static vtx_lookup *create_vtx_lookup_table(ZOLTAN_MP_DATA *mpd,
        IJTYPE *ids1, IJTYPE *ids2, IJTYPE len1, IJTYPE len2);
static void free_vtx_lookup_table(vtx_lookup **vl);
static int lookup_vtx(vtx_lookup *vl, IJTYPE vtx_id);

static int make_hypergraph(ZZ *zz, ZOLTAN_MP_DATA *mpd);
#endif

static int process_matrix_input(ZZ *zz, ZOLTAN_MP_DATA *mpd);
static int make_mirror(ZOLTAN_MP_DATA *mpd);
static int allocate_copy(IJTYPE **to, IJTYPE *from, IJTYPE len);

static ZOLTAN_MP_DATA *MP_Initialize_Structure();
static int MP_Initialize_Params(ZZ *zz, ZOLTAN_MP_DATA *mpd);

static int phg_rows_or_columns(ZOLTAN_MP_DATA *mpd);
static int phg_rows_or_columns_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  unsigned int *export_procs, unsigned int *export_to_part);

/****************************************************************/
/****************************************************************/
/* API:
 *
 *  Zoltan_Matrix_Partition()     must be called by all processes
 *
 *  Queries (not all queries are defined for all LB_APPROACH):
 *
 *  Zoltan_MP_Get_Local_Rows()   assignment of rows of my non-zeroes
 *  Zoltan_MP_Get_Local_Cols()  assignment of columns of my non-zeroes
 *  Zoltan_MP_Get_Local_NZ()      assignment of my non-zeroes
 *
 *  To evaluate partitioning:
 *
 *  Zoltan_Matrix_Partition_Eval()   must be called by all procs
 */
int Zoltan_Matrix_Partition(ZZ *zz)
{
  char *yo = "Zoltan_Matrix_Partition";
  int ierr = ZOLTAN_OK;
  ZZ *zzLib = NULL;
  IJTYPE npins;

  /* The persistent structure we will save at zz->LB.Data_Structure */

  ZOLTAN_MP_DATA *mpd = MP_Initialize_Structure();

  if (mpd == NULL){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

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
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Required query functions are not defined.\n");
    goto End;
  }

  if (mpd->input_type == COL_TYPE){
    zz->Get_CSC_Size(zz->Get_CSC_Size_Data, &mpd->numRC, &npins, &ierr);
  }
  else {
    zz->Get_CSR_Size(zz->Get_CSR_Size_Data, &mpd->numRC, &npins, &ierr);
  }

  mpd->rcGID = (IJTYPE *)ZOLTAN_MALLOC(mpd->numRC * sizeof(IJTYPE));
  mpd->pinIndex = (IJTYPE *)ZOLTAN_MALLOC((mpd->numRC+1) * sizeof(IJTYPE));
  mpd->pinGID = (IJTYPE *)ZOLTAN_MALLOC(npins * sizeof(IJTYPE));
  if ((mpd->numRC && !mpd->rcGID) || !mpd->pinIndex || (npins && !mpd->pinGID)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  if (mpd->input_type == COL_TYPE){
    zz->Get_CSC(zz->Get_CSC_Data, mpd->numRC, npins, 
                mpd->rcGID, mpd->pinIndex, mpd->pinGID, &ierr);
  }
  else {
    zz->Get_CSR(zz->Get_CSR_Data, mpd->numRC, npins, 
                mpd->rcGID, mpd->pinIndex, mpd->pinGID, &ierr);
  }
 
  mpd->pinIndex[mpd->numRC] = npins;

  /*
   * Determine globally how many rows and columns were returned, and if
   * their IDs are 0 based or 1 based.  (Assume they are contiguous.)
   * Some row/column numbers may not appear in input if they have no pins.
   */

  ierr = process_matrix_input(zz, mpd);

  /*
   * Set the query functions and parameters required for this problem.
   * Precompute hypergraph if necessary.  This is the code that needs
   * to be written if you are implementing a new LB_APPROACH.
   */
  switch (mpd->approach)
    {
      case PHG_ROWS:
      case PHG_COLUMNS:
        phg_rows_or_columns(mpd);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
        ierr = ZOLTAN_FATAL;
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

  /*
   * Save the data required to respond to the queries that are
   * supported by each LB_APPROACH.
   */

  switch (mpd->approach)
    {
      case PHG_ROWS:
      case PHG_COLUMNS:
        /*
         * Save row or column assignment for all rows or columns
         * in my non-zeroes.  Each of my non-zeroes is assigned to
         * the row or column that it's row or column is assigned to.
         *
         * PHG_ROWS - we don't repond to Local_Cols query.
         * PHG_COLUMNS - we don't repond to Local_Rows query.
         */
        phg_rows_or_columns_result(mpd, num_export,
            export_global_ids, export_local_ids, export_procs, export_to_part);
        break;
      default:
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid approach.\n");
        ierr = ZOLTAN_FATAL;
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
    ZOLTAN_FREE(&mpd->pinIndex);
    ZOLTAN_FREE(&mpd->pinGID);

    ZOLTAN_FREE(&mpd->crGID);
    ZOLTAN_FREE(&mpd->mirrorPinIndex);
    ZOLTAN_FREE(&mpd->mirrorPinGID);

    ZOLTAN_FREE(&mpd->vtxGID);
    ZOLTAN_FREE(&mpd->vtxWgt);
    ZOLTAN_FREE(&mpd->hindex);
    ZOLTAN_FREE(&mpd->hvertex);

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

  Zoltan_Bind_Param(MP_params, "LB_APPROACH", approach);  
  Zoltan_Bind_Param(MP_params, "NUM_GID_ENTRIES", &mpd->gidLen);  

  strncpy(approach, "PHG_ROWS", MAX_PARAM_STRING_LEN);
  mpd->gidLen = 1;

  ierr = Zoltan_Assign_Param_Vals(zz->Params, MP_params, zz->Debug_Level,
          zz->Proc, zz->Debug_Proc);

  if (ierr == ZOLTAN_OK){
    /* Update this if you add a new LB_APPROACH */
    if (!strcasecmp(approach, "phg_rows"))
      mpd->approach = PHG_ROWS;
    else if (!strcasecmp(approach, "phg_columns"))
      mpd->approach = PHG_COLUMNS;

    if (mpd->gidLen > sizeof(IJTYPE)){
      /*
       * Maybe they want to use long ints, we only handle ints.  Change
       * IJTYPE to make code use long ints.
       */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input row and column IDs are too large");
      ierr = ZOLTAN_FATAL;
    }
  }
  return ierr;
}


static int process_matrix_input(ZZ *zz, ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;
  int i;
  long int lpins, gpins;
  IJTYPE minID, maxID, minPinID, maxPinID;
  IJTYPE npins=0;
  long int vals[2], gvals[2];

  /* get global number of rows, columns and pins, range of row and column IDs */

  if (mpd->numRC > 0)
    npins = mpd->pinIndex[mpd->numRC];

  lpins = (long int)npins;

  MPI_Allreduce(&lpins, &gpins, 1, MPI_LONG, MPI_SUM, zz->Communicator);

  mpd->nNonZeros = (IJTYPE)gpins;

  maxID = maxPinID = -1;

  if (mpd->numRC > 0){
    minID = maxID = mpd->rcGID[0];

    for (i=1; i<mpd->numRC; i++){
      if (mpd->rcGID[i] < minID) minID = mpd->rcGID[i];
      else if (mpd->rcGID[i] > maxID) maxID = mpd->rcGID[i];
    }
    if (npins > 0){
      minPinID = maxPinID = mpd->pinGID[0];
      for (i=1; i<npins; i++){
        if (mpd->pinGID[i] < minPinID) minPinID = mpd->pinGID[i];
        else if (mpd->pinGID[i] > maxPinID) maxPinID = mpd->pinGID[i];
      }
    }
  }
  vals[0] = (long int)maxID;
  vals[1] = (long int)maxPinID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MAX, zz->Communicator);

  maxID = (IJTYPE)gvals[0];
  maxPinID = (IJTYPE)gvals[1];

  if (npins == 0){
    minPinID = maxPinID;
    if (mpd->numRC == 0){
      minID = maxID;
    }
  }
  vals[0] = (long int)minID;
  vals[1] = (long int)minPinID;

  MPI_Allreduce(vals, gvals, 2, MPI_LONG, MPI_MIN, zz->Communicator);

  minID = (IJTYPE)gvals[0];
  minPinID = (IJTYPE)gvals[1];

  if (mpd->input_type == ROW_TYPE){
    mpd->rowBaseID = minID;
    mpd->colBaseID = minPinID;
    mpd->nRows = maxID - minID + 1;
    mpd->nCols = maxPinID - minPinID + 1;
  }
  else{
    mpd->rowBaseID = minPinID;
    mpd->colBaseID = minID;
    mpd->nRows = maxPinID - minPinID + 1;
    mpd->nCols = maxID - minID + 1;
  }

  return ierr;
}
/* TODO - need to copy structure too.
 */
static ZOLTAN_MP_DATA *MP_Initialize_Structure()
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)ZOLTAN_CALLOC(sizeof(ZOLTAN_MP_DATA), 1);

  if (mpd == NULL){
    return NULL;
  }

  return mpd;
}
static int make_mirror(ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;
  int npins, nIDs;

  /*
   * If sparse matrix was supplied in CSC, create another in CSR format,
   * and vice versa.
   */
  ZOLTAN_FREE(&mpd->crGID);
  ZOLTAN_FREE(&mpd->mirrorPinIndex);
  ZOLTAN_FREE(&mpd->mirrorPinGID);

  nIDs = mpd->numRC;
  npins = mpd->pinIndex[nIDs];

  mpd->numCR = mpd->numRC;
  allocate_copy(&mpd->crGID, mpd->pinGID, npins);
  allocate_copy(&mpd->mirrorPinIndex, mpd->pinIndex, nIDs+1);
  allocate_copy(&mpd->mirrorPinGID, mpd->rcGID, nIDs);

  /* TODO  Convert_to_CSR thinks some of these fields are ints
   *         when IJTYPEs might not be ints.  FIX, but OK for now.
   *
   */

  ierr = Zoltan_Convert_To_CSR(mpd->zzLib, npins, (int *)mpd->pinIndex,
             /* The following get overwritten with mirror */
             (int *)&(mpd->numCR),
             (ZOLTAN_ID_PTR *)&(mpd->mirrorPinGID),
             (int **)&(mpd->mirrorPinIndex),
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
/* A search structure and routines for finding a vertex ID      
 */
#if 0
static void free_vtx_lookup_table(vtx_lookup **lu)
{
  vtx_lookup *l = *lu;

  if (l == NULL) return;

  ZOLTAN_FREE(&l->htTop);
  ZOLTAN_FREE(&l->ht);
  ZOLTAN_FREE(lu);
}

/*
 * Look up array index (local id) for vertex global ID
 */
static int lookup_vtx(vtx_lookup *lu, IJTYPE vtxGID)
{
  char *yo="lookup_vtx";
  struct vtx_node *hn;
  int i, num_id_entries;

  if (lu->table_size < 1) return -1;
  if (lu->numVtx < 1) return -1;

  if (num_id_entries < 1){
    ZOLTAN_PRINT_ERROR(0, yo, "code rewrite required.\n");
    exit(0);  /* this should not happen, but just in case... */
  }

  i = Zoltan_Hash((ZOLTAN_ID_PTR)vtxGID, lu->idLength, (unsigned int) lu->table_size);

  for (hn=lu->ht[i]; hn != NULL; hn = hn->next){
    if (hn->vtxGID == vtxGID){
      return hn->vtxLID;
    }
  }
  return -1;
}

/*
 * Given two lists of vertex ids, map them to local ids (array indices).  The
 * ids in a list may be repeated, the ids in each list are disjoint sets.
 * Create a lookup table.
 */
static vtx_lookup *create_vtx_lookup_table(ZOLTAN_MP_DATA *mpd,
          IJTYPE *ids1, IJTYPE *ids2, IJTYPE len1, IJTYPE len2)
{
  int i, ii, j;
  struct vtx_node *hn;
  vtx_lookup *lu = NULL;
  IJTYPE tsize, numids, len, found;
  IJTYPE *ids;

  lu = (vtx_lookup *)ZOLTAN_MALLOC(sizeof(vtx_lookup));
  if (!lu){
    return NULL;
  }

  lu->idLength = mpd->gidLen;
  numids = len1 + len2;
  tsize = numids / 4;
  if (tsize < 10) tsize = 10;

  lu->ht = (struct vtx_node **)ZOLTAN_CALLOC(sizeof(struct vtx_node*) , tsize);
  hn = lu->htTop = (struct vtx_node *)ZOLTAN_MALLOC(sizeof(struct vtx_node) * numids);

  if (tsize && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu);
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    return NULL;
  }

  lu->table_size = tsize;
  lu->numVtx = 0;

  for (ii=0; ii<2; ii++){
    ids = ((ii) ? ids2 : ids1);
    len = ((ii) ? len2 : len1);
    for (i=0; i<len; i++){
  
      found = lookup_vtx(lu, ids[i]);
  
      if (found >= 0) continue;
  
      hn->vtxGID = ids[i];
      hn->vtxLID = lu->numVtx;

      j = Zoltan_Hash((ZOLTAN_ID_PTR)(ids+i), mpd->gidLen, tsize);
  
      hn->next = lu->ht[j];
      lu->ht[j] = hn;
  
      hn++;
      lu->numVtx++;
    }
  }

  if (lu->numVtx < numids){
    lu->htTop = (struct vtx_node *)ZOLTAN_REALLOC(lu->htTop, sizeof(struct vtx_node) * lu->numVtx);
  }

  return lu;
}
#endif
/****************************************************************/
/****************************************************************/
/* Functions that support LB_APPROACH = PHG_ROWS or PHG_COLS
 *
 * The objects to be partitioned are the matrix rows or columns,
 * and we are using Zoltan's PHG to do it.
 *
 * PHG_ROWS: the objects are the sparse matrix rows, and the
 *  hyperedges are the columns of the sparse matrix.
 *
 * PHG_COLS: the objects are the sparse matrix columns, and
 *  the hyperedges are the rows.
 */

static int my_objects(ZOLTAN_MP_DATA *mpd, IJTYPE *nobj, IJTYPE **objIDs)
{
  int nProcs, me;
  IJTYPE i, n, *obj;
  int ierr = ZOLTAN_OK;

  nProcs = mpd->zzLib->Num_Proc;
  me = mpd->zzLib->Proc;
  n = 0;
  obj=NULL;

  if (mpd->approach == PHG_ROWS){
    for (i=mpd->rowBaseID; i <mpd->rowBaseID + mpd->nRows; i++){
      if (i % nProcs == me){
        n++;
      }
    }
  }
  else{
    for (i=mpd->colBaseID; i <mpd->colBaseID + mpd->nCols; i++){
      if (i % nProcs == me){
        n++;
      }
    }
  }

  obj = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * n);
  if (n && !obj){
    ierr = ZOLTAN_MEMERR;
  }
  else{
    *objIDs = obj;
    *nobj = n;

    if (mpd->approach == PHG_ROWS){
      for (i=mpd->rowBaseID; i <mpd->rowBaseID + mpd->nRows; i++){
        if (i % nProcs == me){
          *obj++ = i;
        }
      }
    }
    else{
      for (i=mpd->colBaseID; i <mpd->colBaseID + mpd->nCols; i++){
        if (i % nProcs == me){
          *obj++ = i;
        }
      }
    }
  }

  return ierr;
}

static int phg_rows_or_columns(ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;

  Zoltan_Set_Param(mpd->zzLib, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(mpd->zzLib, "NUM_LID_ENTRIES", "1");

  Zoltan_Set_Param(mpd->zzLib, "LB_METHOD", "HYPERGRAPH");
  Zoltan_Set_Param(mpd->zzLib, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param(mpd->zzLib, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(mpd->zzLib, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(mpd->zzLib, "ADD_OBJ_WEIGHT", "PINS");
  Zoltan_Set_Param(mpd->zzLib, "EDGE_WEIGHT_DIM", "0");

  Zoltan_Set_Num_Obj_Fn(mpd->zzLib, Zoltan_MP_Get_Num_Obj, mpd);
  Zoltan_Set_Obj_List_Fn(mpd->zzLib, Zoltan_MP_Get_Obj_List, mpd);
  Zoltan_Set_HG_Size_CS_Fn(mpd->zzLib, Zoltan_MP_Get_Size_Matrix, mpd);
  Zoltan_Set_HG_CS_Fn(mpd->zzLib, Zoltan_MP_Get_Matrix, mpd);

  ierr = make_mirror(mpd);

  if (ierr == ZOLTAN_OK){
    ierr = my_objects(mpd, &mpd->nMyVtx, &mpd->vtxGID);
  }

  return ierr;
}

int Zoltan_MP_Get_Num_Obj(void *data, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;

  return mpd->nMyVtx;
}

void Zoltan_MP_Get_Obj_List(void *data, int num_gid_entries, int num_lid_entries,
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

void Zoltan_MP_Get_Size_Matrix(void *data, int *num_lists, int *num_pins,
  int *format, int *ierr)
{
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;

  if (((mpd->approach == PHG_ROWS) && (mpd->input_type == ROW_TYPE)) ||
      ((mpd->approach == PHG_COLUMNS) && (mpd->input_type == COL_TYPE)) ){
    *num_lists = mpd->numCR;
    *num_pins  = mpd->mirrorPinIndex[mpd->numCR];
  }
  else{
    *num_lists = mpd->numRC;
    *num_pins  = mpd->pinIndex[mpd->numRC];
  }
  *format = ZOLTAN_COMPRESSED_EDGE;
}

void Zoltan_MP_Get_Matrix(void *data, int num_gid_entries, int num_vtx_edge,
  int num_pins, int format, 
  ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, ZOLTAN_ID_PTR pin_GID, int *ierr)
{
  int i;
  ZOLTAN_MP_DATA *mpd = (ZOLTAN_MP_DATA *)data;
  *ierr = ZOLTAN_OK;
  IJTYPE *edgeGID = (IJTYPE *)vtxedge_GID;
  IJTYPE *pinGID = (IJTYPE *)pin_GID;

  if (((mpd->approach == PHG_ROWS) && (mpd->input_type == ROW_TYPE)) ||
      ((mpd->approach == PHG_COLUMNS) && (mpd->input_type == COL_TYPE)) ){
    for (i=0; i<mpd->numCR; i++){
      edgeGID[i] = mpd->crGID[i];
      vtxedge_ptr[i] = mpd->mirrorPinIndex[i];
    }
    for (i=0; i< mpd->mirrorPinIndex[mpd->numCR]; i++){
      pinGID[i] = mpd->mirrorPinGID[i];
    }
  }
  else{
    for (i=0; i<mpd->numRC; i++){
      edgeGID[i] = mpd->rcGID[i];
      vtxedge_ptr[i] = mpd->pinIndex[i];
    }
    for (i=0; i< mpd->pinIndex[mpd->numRC]; i++){
      pinGID[i] = mpd->pinGID[i];
    }
  }
}

static int phg_rows_or_columns_result(ZOLTAN_MP_DATA *mpd,
  int num_export,
  unsigned int *export_global_ids, unsigned int *export_local_ids, 
  unsigned int *export_procs, unsigned int *export_to_part)
{
  char *yo = "phg_rows_or_columns_result";
  int ierr = ZOLTAN_OK;

  return ierr;
}


/****************************************************************/
/****************************************************************/
/* Code to make a hypergraph out of sparse matrix by creating an
   object from each row and column, and a hyperedge from each
   row and column.  This code is not used currently.  The idea
   is to partition both rows and columns where the objects are
   the rows and columns themselves. 

   This code is not tested.
*/

#if 0

static IJTYPE vertex_GID(IJTYPE id, int rc, ZOLTAN_MP_DATA *mpd);
static int object_owner(ZZ *zz, IJTYPE gid);
static int send_recv_rows_columns(int rc, ZZ *zz, ZOLTAN_MP_DATA *mpd, int tag,
  int numRC, IJTYPE *rcGID, IJTYPE *pinIndex, IJTYPE *pinGID,
  int *nrecv, IJTYPE **recvIDs, IJTYPE **recvLen, IJTYPE **recvPins); 

static int make_hypergraph(ZZ *zz, ZOLTAN_MP_DATA *mpd)
{
  int ierr = ZOLTAN_OK;
  int i, j, k, lid, size_hvertex;
  IJTYPE numIDs;
  int nrecvA, nrecvB;
  IJTYPE *recvIDsA=NULL, *recvIDsB=NULL;
  IJTYPE *recvLenA=NULL, *recvLenB=NULL;
  IJTYPE *recvPinsA=NULL, *recvPinsB=NULL;
  IJTYPE *hesize=NULL;
  IJTYPE *IDs, *pins, *nPins, *vptr;
  vtx_lookup *lu_table=NULL;

  /*
   * The user has input a sparse input in CSC or CSR format.  Create the
   * mirror CSR or CSC format.  
   */
  ierr = make_mirror(mpd);

  if (ierr != ZOLTAN_OK){
    return ierr;
  }

  /* Send rows and columns to owners, and convert row and column IDs
   * to hypergraph vertex IDs.
   */

  ierr = send_recv_rows_columns(mpd->input_type, zz, mpd, 40000,
      mpd->numRC, mpd->rcGID, mpd->pinIndex, mpd->pinGID,
      &nrecvA, &recvIDsA, &recvLenA, &recvPinsA);

  ierr = send_recv_rows_columns(
    (mpd->input_type == ROW_TYPE ? COL_TYPE : ROW_TYPE), 
    zz, mpd, 39000,
    mpd->numCR, mpd->crGID, mpd->mirrorPinIndex, mpd->mirrorPinGID,
    &nrecvB, &recvIDsB, &recvLenB, &recvPinsB);


 /* 
  * Create local portion of hypergraph from all data just received. 
  */ 

  lu_table = create_vtx_lookup_table(mpd, recvIDsA, recvIDsB, nrecvA, nrecvB);

  mpd->nHedges = mpd->nMyVtx = lu_table->numVtx;

  mpd->vtxGID = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * mpd->nMyVtx);
  mpd->hindex = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * (mpd->nMyVtx + 1));

  /* count size of each hyperedge */

  hesize = (IJTYPE *)ZOLTAN_CALLOC(sizeof(IJTYPE) , mpd->nHedges);

  for (i=0; i<2; i++){
    numIDs = (i ? nrecvB : nrecvA);
    IDs = (i ? recvIDsB : recvIDsA);
    nPins = (i ? recvLenB: recvLenA);
    pins = (i ? recvPinsB : recvPinsA);

    for (j=0; j < numIDs; j++){
      lid = lookup_vtx(lu_table, IDs[j]);
      hesize[lid] += nPins[j];  /* # pins I recv'd for this row or column */
    }
  }

  for (i=0; i<mpd->nHedges; i++){
    /* "+ 1" because vtx is also part of assoc. hyperedge */
    mpd->hindex[i+1] = mpd->hindex[i] + hesize[i] + 1; 
  }
  size_hvertex = mpd->hindex[mpd->nHedges];

  /* write out vertices in each hyperedge */

  mpd->hvertex = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * size_hvertex);
  memset(hesize, 0, sizeof(int) * mpd->nHedges);

  for (i=0; i<2; i++){
    numIDs = (i ? nrecvB : nrecvA);
    IDs = (i ? recvIDsB : recvIDsA);
    nPins = (i ? recvLenB: recvLenA);
    pins = (i ? recvPinsB : recvPinsA);

    for (j=0; j < numIDs; j++){
      lid = lookup_vtx(lu_table, IDs[j]);
      vptr = mpd->hvertex + mpd->hindex[lid] + hesize[lid];
      if (hesize[lid] == 0){
        *vptr++ = IDs[j];      /* vertex is part of its assoc. hyperedge */
        hesize[lid] += 1;
      }
      for (k=0; k < nPins[j]; k++){
        *vptr++ = *pins++;
      }
      hesize[lid] += nPins[j];
    }
  }

  /* Set vertex weight 
   *    We'll use ADD_OBJ_WEIGHT for this, let zoltan calculate weight
   *

  mpd->vtxWgt = (double *)ZOLTAN_MALLOC(sizeof(double) * mpd->nMyVtx);
  for (i=0; i<mpd->nMyVtx; i++){
    mpd->vtxWgt[i] = (double)(mpd->hindex[i+1] - mpd->hindex[i]);
  }
   */

  ZOLTAN_FREE(&hesize);
  ZOLTAN_FREE(&recvIDsA);
  ZOLTAN_FREE(&recvIDsB);
  ZOLTAN_FREE(&recvLenA);
  ZOLTAN_FREE(&recvLenB);
  ZOLTAN_FREE(&recvPinsA);
  ZOLTAN_FREE(&recvPinsB);

 /*
  * TODO: can we free the input sparse matrix and/or mirror now, or
  * do we need them to reply to queries.
  */

  /*
   * Now set parameters and query functions
   */

  Zoltan_Set_Param(mpd->zzLib, "LB_METHOD", "HYPERGRAPH"); 
  Zoltan_Set_Param(mpd->zzLib, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param(mpd->zzLib, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(mpd->zzLib, "RETURN_LISTS", "ALL");    /* not sure we need "all" */
  Zoltan_Set_Param(mpd->zzLib, "OBJ_WEIGHT_DIM", "0"); 
  Zoltan_Set_Param(mpd->zzLib, "EDGE_WEIGHT_DIM", "0");
  Zoltan_Set_Param(mpd->zzLib, "ADD_OBJ_WEIGHT", "pins");

  return ierr;
}

static IJTYPE vertex_GID(IJTYPE id, int rc, ZOLTAN_MP_DATA *mpd)
{
  int vtxGID;
  /*
   * Each row and column of the sparse matrix will yield a vertex and
   * hyperedge of the hypergraph.
   */
  if (rc == ROW_TYPE){       /* "id" is a row GID */
    if ( (id >= mpd->rowBaseID) && (id < (mpd->rowBaseID + mpd->nRows))){
      vtxGID = id - mpd->rowBaseID;
    }
  }
  else if (rc == COL_TYPE){  /* "id" is a column GID */
    if ( (id >= mpd->colBaseID) && (id < (mpd->colBaseID + mpd->nCols))){
      vtxGID = mpd->nRows + (id - mpd->colBaseID);
    }
  }
  return vtxGID;
}

static int object_owner(ZZ *zz, IJTYPE gid)
{
  /* If there are many rows/cols with no pins, we may get
   * some imbalance here.  If this is a problem, we need
   * to do more work to determine what all the GIDs are
   * and map them to consecutive global numbers.
   */
  int owner = gid % zz->Num_Proc;
  return owner;
}

static int send_recv_rows_columns(int rc, ZZ *zz, ZOLTAN_MP_DATA *mpd, int tag,
  int numRC, IJTYPE *rcGID, IJTYPE *pinIndex, IJTYPE *pinGID,           /* input */
  int *nrecv, IJTYPE **recvIDs, IJTYPE **recvLen, IJTYPE **recvPins) /* output */
{
  int i;
  IJTYPE vtxgid, len;
  IJTYPE *ids;
  int ierr = ZOLTAN_OK;
  int *owner = NULL;
  IJTYPE *numNonZeros = NULL;
  ZOLTAN_COMM_OBJ *plan;

  int cr = ((rc == ROW_TYPE) ? COL_TYPE : ROW_TYPE);

  /*
   * Each row and column of the input sparse matrix is a "vertex"
   * in the hypergraph.  And each row and column is a hyperedge,
   * it's vertices being the non-zeros in the row or column.  We
   * send rows and columns to their owners.
   */

  owner = (int *)ZOLTAN_MALLOC(sizeof(int) * numRC);
  numNonZeros = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * numRC);

  for (i=0; i<numRC; i++){
    vtxgid = vertex_GID(rcGID[i], rc, mpd);
    owner[i] = object_owner(zz, vtxgid);

    numNonZeros[i] = pinIndex[i+1] - pinIndex[i];
  }

  ierr = Zoltan_Comm_Create(&plan, numRC, owner, zz->Communicator, tag, (int *)&len);

  /* Send row or column numbers to owners */

  *nrecv = len;
  *recvIDs = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * len);

  ids = *recvIDs;

  tag--;

  ierr = Zoltan_Comm_Do(plan, tag, (char *)rcGID, sizeof(IJTYPE), (char *)ids);

  /* Convert input matrix row/column IDs to hg vertex IDs */
  
  for (i=0; i<len; i++){
    ids[i] = vertex_GID(ids[i], rc, mpd);
  }

  /* Send number of non-zeros in row or column */

  *recvLen = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * len);

  tag--;

  ierr = Zoltan_Comm_Do(plan, tag, (char *)numNonZeros, sizeof(IJTYPE), (char *)*recvLen);

  /* Send pins in row or column */

  tag--;

  ierr = Zoltan_Comm_Resize(plan, (int *)numNonZeros, tag, (int *)&len);

  *recvPins = (IJTYPE *)ZOLTAN_MALLOC(sizeof(IJTYPE) * len);
  ids = *recvPins;

  tag--;

  ierr = Zoltan_Comm_Do(plan, tag, (char *)pinGID, sizeof(IJTYPE), (char *)ids);

  ierr = Zoltan_Comm_Destroy(&plan);

  /* Convert input matrix row/column IDs to hg vertex IDs */
  
  for (i=0; i<len; i++){
    ids[i] = vertex_GID(ids[i], cr, mpd);
  }

  return ierr;
}

#endif
/****************************************************************/
/****************************************************************/


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_MATRIX_INPUT */
