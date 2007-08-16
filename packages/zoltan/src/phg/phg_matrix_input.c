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

/* TODO: check all the time whether mmd->numRC is zero before proceeding */

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"
#include "zz_util_const.h"

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
**
** This source file uses "pins" and "non-zeros" interchangeably.
*************************************************************************/
struct vtx_node{
  int vtxGID;
  int vtxLID;
  struct vtx_node *next;
};
typedef struct _vtx_lookup{
  struct vtx_node *htTop;
  struct vtx_node **ht;
  int table_size;
  int numVtx;
}vtx_lookup;

static vtx_lookup *create_vtx_lookup_table(int *ids1, int *ids2, int len1, int len2);
static void free_vtx_lookup_table(vtx_lookup **vl);
static int lookup_vtx(vtx_lookup *vl, int vtx_id);

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

  /*
   * Create the Zoltan problem that we will be solving with the
   * input sparse matrix.
   */

  zzLib = Zoltan_Create(zz->Communicator);
  Zoltan_Set_Param(zzLib, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zzLib, "NUM_LID_ENTRIES", "1");
  mmd->zzLib = zzLib;

  /* Call application defined query functions to get the non-zeros.
   * User has a choice of giving us compressed rows or compressed columns.
   * We need global row and column IDs.  Any process can give us any of the
   * non-zeros, but each non-zero must be supplied by only one process.
   *
   * I'm assuming row and column IDs are ints, not more general objects.
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
    zz->Get_CSC_Size(zz->Get_CSC_Size_Data, &mmd->numRC, &npins, &ierr);
  }
  else {
    zz->Get_CSR_Size(zz->Get_CSR_Size_Data, &mmd->numRC, &npins, &ierr);
  }

  mmd->rcGID = (int *)ZOLTAN_MALLOC(mmd->numRC * sizeof(int));
  mmd->pinIndex = (int *)ZOLTAN_MALLOC((mmd->numRC+1) * sizeof(int));
  mmd->pinGID = (int *)ZOLTAN_MALLOC(npins * sizeof(int));
  if ((mmd->numRC && !mmd->rcGID) || !mmd->pinIndex || (npins && !mmd->pinGID)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  if (mmd->input_type == COL_TYPE){
    zz->Get_CSC(zz->Get_CSC_Data, mmd->numRC, npins, 
                mmd->rcGID, mmd->pinIndex, mmd->pinGID, &ierr);
  }
  else {
    zz->Get_CSR(zz->Get_CSR_Data, mmd->numRC, npins, 
                mmd->rcGID, mmd->pinIndex, mmd->pinGID, &ierr);
  }
 
  mmd->pinIndex[mmd->numRC] = npins;

  /*
   * Determine globally how many rows and columns were returned, and if
   * their IDs are 0 based or 1 based.  (Assume they are contiguous.)
   * Some row/column numbers may not appear in input if they have no pins.
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

  Zoltan_Set_Param(zzLib, "LB_METHOD", "HYPERGRAPH"); 
  Zoltan_Set_Param(zzLib, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param(zzLib, "LB_APPROACH", "PARTITION");
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

static int vertex_GID(int id, int rc, ZOLTAN_MM_DATA *mmd)
{
  int vtxGID=-1;
  /*
   * Each row and column of the sparse matrix will yield a vertex and
   * hyperedge of the hypergraph.
   */
  if (rc == ROW_TYPE){       /* "id" is a row GID */
    if ( (id >= mmd->rowBaseID) && (id < (mmd->rowBaseID + mmd->nRows))){
      vtxGID = id - mmd->rowBaseID;
    }
  }
  else if (rc == COL_TYPE){  /* "id" is a column GID */
    if ( (id >= mmd->colBaseID) && (id < (mmd->colBaseID + mmd->nCols))){
      vtxGID = mmd->nRows + (id - mmd->colBaseID);
    }
  }
  return vtxGID;
}

static int object_owner(ZZ *zz, int gid)
{
  /* If there are many rows/cols with no pins, we may get
   * some imbalance here.  If this is a problem, we need
   * to do more work to determine what all the GIDs are
   * and map them to consecutive global numbers.
   */
  int owner = gid % zz->Num_Proc;
  return owner;
}

static int allocate_copy(int **to, int *from, int len)
{
  *to = (int *)ZOLTAN_MALLOC(sizeof(int) * len);
  memcpy(*to, from, sizeof(int) * len);
  return ZOLTAN_OK;
}

static int send_recv_rows_columns(int rc, ZZ *zz, ZOLTAN_MM_DATA *mmd, int tag,
  int numRC, int *rcGID, int *pinIndex, int *pinGID,           /* input */
  int *nrecv, int **recvIDs, int **recvLen, int **recvPins) /* output */
{
  int i, vtxgid, len;
  int ierr = ZOLTAN_OK;
  int *owner = NULL;
  int *numNonZeros = NULL;
  int *ids;
  ZOLTAN_COMM_OBJ *plan;

  int cr = ((rc == ROW_TYPE) ? COL_TYPE : ROW_TYPE);

  /*
   * Each row and column of the input sparse matrix is a "vertex"
   * in the hypergraph.  And each row and column is a hyperedge,
   * it's vertices being the non-zeros in the row or column.  We
   * send rows and columns to their owners.
   */

  owner = (int *)ZOLTAN_MALLOC(sizeof(int) * numRC);
  numNonZeros = (int *)ZOLTAN_MALLOC(sizeof(int) * numRC);

  for (i=0; i<numRC; i++){
    vtxgid = vertex_GID(rcGID[i], rc, mmd);
    owner[i] = object_owner(zz, vtxgid);

    numNonZeros[i] = pinIndex[i+1] - pinIndex[i];
  }

  ierr = Zoltan_Comm_Create(&plan, numRC, owner, zz->Communicator, tag, &len);

  /* Send row or column numbers to owners */

  *nrecv = len;
  *recvIDs = (int *)ZOLTAN_MALLOC(sizeof(int) * len);

  ids = *recvIDs;

  tag--;

  ierr = Zoltan_Comm_Do(plan, tag, (char *)rcGID, sizeof(int), (char *)ids);

  /* Convert input matrix row/column IDs to hg vertex IDs */
  
  for (i=0; i<len; i++){
    ids[i] = vertex_GID(ids[i], rc, mmd);
  }

  /* Send number of non-zeros in row or column */

  *recvLen = (int *)ZOLTAN_MALLOC(sizeof(int) * len);

  tag--;

  ierr = Zoltan_Comm_Do(plan, tag, (char *)numNonZeros, sizeof(int), (char *)*recvLen);

  /* Send pins in row or column */

  tag--;

  ierr = Zoltan_Comm_Resize(plan, numNonZeros, tag, &len);

  *recvPins = (int *)ZOLTAN_MALLOC(sizeof(int) * len);
  ids = *recvPins;

  tag--;

  ierr = Zoltan_Comm_Do(plan, tag, (char *)pinGID, sizeof(int), (char *)ids);

  ierr = Zoltan_Comm_Destroy(&plan);

  /* Convert input matrix row/column IDs to hg vertex IDs */
  
  for (i=0; i<len; i++){
    ids[i] = vertex_GID(ids[i], cr, mmd);
  }

  return ierr;
}

static int make_hypergraph(ZZ *zz, ZOLTAN_MM_DATA *mmd)
{
  int ierr = ZOLTAN_OK;
  int i, j, k, numIDs, lid, size_hvertex;
  int nrecvA, nrecvB;
  int *recvIDsA=NULL, *recvIDsB=NULL;
  int *recvLenA=NULL, *recvLenB=NULL;
  int *recvPinsA=NULL, *recvPinsB=NULL;
  int *hesize=NULL;
  int *IDs, *nPins, *pins, *vptr;
  vtx_lookup *lu_table=NULL;
  ZOLTAN_ID_PTR objPtr, pinPtr;

  /*
   * The user has input a sparse input in CSC or CSR format.  Create the
   * mirror CSR or CSC format.  Convert_To_CSR program converts in place.
   */

  mmd->numCR = mmd->numRC;
  allocate_copy(&mmd->crGID, mmd->rcGID, mmd->numRC);
  allocate_copy(&mmd->mirrorPinIndex, mmd->pinIndex, mmd->numRC+1);
  allocate_copy(&mmd->mirrorPinGID, mmd->pinGID, mmd->pinIndex[mmd->numRC]);

  objPtr = (ZOLTAN_ID_PTR)mmd->crGID;
  pinPtr = (ZOLTAN_ID_PTR)mmd->mirrorPinGID;

  ierr = Zoltan_Convert_To_CSR(mmd->zzLib,    /* TODO test this */
             mmd->pinIndex[mmd->numRC],  
             mmd->pinIndex,
             &(mmd->numCR),
             &objPtr,
             &(mmd->mirrorPinIndex),
             &pinPtr);

  mmd->crGID = (int *)objPtr;
  mmd->mirrorPinIndex = (int *)pinPtr;

  /* Send rows and columns to owners, and convert row and column IDs
   * to hypergraph vertex IDs.
   */

  ierr = send_recv_rows_columns(mmd->input_type, zz, mmd, 40000,
      mmd->numRC, mmd->rcGID, mmd->pinIndex, mmd->pinGID,
      &nrecvA, &recvIDsA, &recvLenA, &recvPinsA);

  ierr = send_recv_rows_columns(
    (mmd->input_type == ROW_TYPE ? COL_TYPE : ROW_TYPE), 
    zz, mmd, 39000,
    mmd->numCR, mmd->crGID, mmd->mirrorPinIndex, mmd->mirrorPinGID,
    &nrecvB, &recvIDsB, &recvLenB, &recvPinsB);


 /* 
  * Create local portion of hypergraph from all data just received. 
  */ 

  lu_table = create_vtx_lookup_table(recvIDsA, recvIDsB, nrecvA, nrecvB);

  mmd->nHedges = mmd->nMyVtx = lu_table->numVtx;

  mmd->vtxGID = (int *)ZOLTAN_MALLOC(sizeof(int) * mmd->nMyVtx);
  mmd->vtxWgt = (double *)ZOLTAN_MALLOC(sizeof(double) * mmd->nMyVtx);
  mmd->hindex = (int *)ZOLTAN_MALLOC(sizeof(int) * (mmd->nMyVtx + 1));

  /* count size of each hyperedge */

  hesize = (int *)ZOLTAN_CALLOC(sizeof(int) , mmd->nHedges);

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

  for (i=0; i<mmd->nHedges; i++){
    /* "+ 1" because vtx is also part of assoc. hyperedge */
    mmd->hindex[i+1] = mmd->hindex[i] + hesize[i] + 1; 
  }
  size_hvertex = mmd->hindex[mmd->nHedges];

  /* write out vertices in each hyperedge */

  mmd->hvertex = (int *)ZOLTAN_MALLOC(sizeof(int) * size_hvertex);
  memset(hesize, 0, sizeof(int) * mmd->nHedges);

  for (i=0; i<2; i++){
    numIDs = (i ? nrecvB : nrecvA);
    IDs = (i ? recvIDsB : recvIDsA);
    nPins = (i ? recvLenB: recvLenA);
    pins = (i ? recvPinsB : recvPinsA);

    for (j=0; j < numIDs; j++){
      lid = lookup_vtx(lu_table, IDs[j]);
      vptr = mmd->hvertex + mmd->hindex[lid] + hesize[lid];
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

  /* Set vertex weight to size of associated hyperedge? */

  for (i=0; i<mmd->nMyVtx; i++){
    mmd->vtxWgt[i] = (double)(mmd->hindex[i+1] - mmd->hindex[i]);
  }

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

  return ierr;
}

static int process_matrix_input(ZZ *zz, ZOLTAN_MM_DATA *mmd)
{
  int ierr = ZOLTAN_OK;
  int i;
  int minID, maxID, minPinID, maxPinID;
  int npins=0;
  int vals[2], gvals[2];

  /* get global number of rows, columns and pins, range of row and column IDs */

  if (mmd->numRC > 0)
    npins = mmd->pinIndex[mmd->numRC];

  MPI_Allreduce(&npins, &mmd->nNonZeros, 1, MPI_LONG, MPI_SUM, zz->Communicator);

  maxID = maxPinID = -1;

  if (mmd->numRC > 0){
    minID = maxID = mmd->rcGID[0];

    for (i=1; i<mmd->numRC; i++){
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
    if (mmd->numRC == 0){
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
  mmd->numRC = 0;
  mmd->rcGID = NULL;
  mmd->pinIndex = NULL;
  mmd->pinGID = NULL;

  mmd->rowBaseID=0;
  mmd->colBaseID=0;
  mmd->nRows=0;
  mmd->nCols=0;
  mmd->nNonZeros=0;

  mmd->zzLib = NULL;

  return mmd;
}

/**********************************************************
 * A search structure and routines for finding a vertex ID
 */

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
static int lookup_vtx(vtx_lookup *lu, int vtxGID)
{
  struct vtx_node *hn;
  int i;
  unsigned int dummy;

  if (lu->table_size < 1) return -1;
  if (lu->numVtx < 1) return -1;

  dummy = (unsigned int)vtxGID;

  i = Zoltan_Hash(&dummy, 1, (unsigned int) lu->table_size);

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
static vtx_lookup *create_vtx_lookup_table(int *ids1, int *ids2, int len1, int len2)
{
  int i, ii, j, tsize, numids, found, len;
  unsigned int lookup_val;
  struct vtx_node *hn;
  vtx_lookup *lu = NULL;
  int *ids;

  lu = (vtx_lookup *)ZOLTAN_MALLOC(sizeof(vtx_lookup));
  if (!lu){
    return NULL;
  }

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

      lookup_val = (unsigned int)ids[i];
  
      j = Zoltan_Hash(&lookup_val, 1, tsize);
  
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

/**********************************************************/
int Zoltan_Matrix_Multiply_Free(ZZ *zz)
{
  /* Free our structures */
  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_MATRIX_INPUT */
