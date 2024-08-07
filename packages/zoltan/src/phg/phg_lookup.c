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

#include "phg_lookup.h"
#include "zz_util_const.h"

/*****************************************************************************/
void phg_free_objects(zoltan_objects *zo)
{
  if (zo == NULL) return;
  ZOLTAN_FREE(&(zo->vtxHash));
}
/*****************************************************************************/
void phg_free_pins(zoltan_pins *zp)
{
  if (zp == NULL) return;
  ZOLTAN_FREE(&(zp->edgeGID));
  ZOLTAN_FREE(&(zp->esizes));
  ZOLTAN_FREE(&(zp->pinGID));
  ZOLTAN_FREE(&(zp->edgeHash));
}
/*****************************************************************************/
void phg_free_ews(zoltan_ews *zew)
{
  if (zew == NULL) return;
  ZOLTAN_FREE(&(zew->edgeGID));
  ZOLTAN_FREE(&(zew->edgeHash));
  ZOLTAN_FREE(&(zew->wgt));
}
/*****************************************************************************/
void phg_free_temp_edges(zoltan_temp_edges *zte)
{
  if (zte == NULL) return;
  ZOLTAN_FREE(&(zte->edgeGID));
  ZOLTAN_FREE(&(zte->pinGID));
  ZOLTAN_FREE(&(zte->pinHash));
}
/*****************************************************************************/
void phg_free_temp_vertices(zoltan_temp_vertices *ztv)
{
  if (ztv == NULL) return;
  ZOLTAN_FREE(&(ztv->vtxGID));
  ZOLTAN_FREE(&(ztv->vtxOwner));
  ZOLTAN_FREE(&(ztv->vtxGNO));
}

/*****************************************************************************/
int phg_map_GIDs_to_processes(ZZ *zz, ZOLTAN_ID_PTR eid, int size, 
                             int lenGID, int **hashedProc, int nprocs)
{
int i, j;
int *procList;
static char *yo = "map_GIDs_to_processes";

  *hashedProc = NULL;

  if (size < 1){
    return ZOLTAN_OK;
  }

  procList = (int *)ZOLTAN_MALLOC(sizeof(int) * size);

  if (!procList){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); 
    return ZOLTAN_MEMERR;
  }

  for (i=0; i<size; i++){
    j = Zoltan_Hash(eid, lenGID, nprocs);
    procList[i] = j;
    eid += lenGID;
  }

  *hashedProc = procList;

  return ZOLTAN_OK;
}

/****************************************************************************/
/* 
 * Create, access and delete a hash table mapping a GID its
 * position in a list.
 */

void phg_free_GID_lookup_table(phg_GID_lookup **lu)
{
  phg_GID_lookup *l = *lu;

  if (l == NULL) return;

  ZOLTAN_FREE(&l->htTop);
  ZOLTAN_FREE(&l->ht);
  ZOLTAN_FREE(lu);
}
/****************************************************************************/
phg_GID_lookup *phg_create_GID_lookup_table(ZOLTAN_ID_PTR gids, int size, int lenGID)
{
  int i, j, tsize;
  phg_GID_lookup *lu = NULL;

  lu = (phg_GID_lookup *)ZOLTAN_MALLOC(sizeof(phg_GID_lookup));
  if (!lu){
    return NULL;
  }

  tsize = size * 1.25;

  lu->htTop = (struct Hash_Node *)ZOLTAN_MALLOC(sizeof(struct Hash_Node)*size);
  lu->ht = (struct Hash_Node **)ZOLTAN_CALLOC(sizeof(struct Hash_Node*), tsize);

  if (tsize && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    ZOLTAN_FREE(&lu);
    return NULL;
  }

  lu->table_size = tsize;
  lu->numGIDs = size;
  lu->lenGID = lenGID;

  for (i=0; i<size; i++){
    lu->htTop[i].gid = gids + (i * lenGID);
    lu->htTop[i].gno = i;

    j = Zoltan_Hash(lu->htTop[i].gid, lenGID, tsize);
    
    lu->htTop[i].next = lu->ht[j];
    lu->ht[j] = lu->htTop + i;
  }
 
  return lu;
}

/****************************************************************************/
/* gid list is not unique.  rewrite it as a list of unique gids */

phg_GID_lookup *phg_create_GID_lookup_table2(ZOLTAN_ID_PTR gids, int ngids, int lenGID)
{
  int i, j, k, tsize, found;
  struct Hash_Node *hn;
  ZOLTAN_ID_PTR nextGID, nextUniqueGID;
  phg_GID_lookup *lu = NULL;

  tsize = ngids;    /* actually may be larger than number of unique ids */
  
  nextGID = nextUniqueGID = gids;

  lu = (phg_GID_lookup *)ZOLTAN_MALLOC(sizeof(phg_GID_lookup));
  if (!lu){
    return NULL;
  }

  lu->ht = (struct Hash_Node **)ZOLTAN_CALLOC(sizeof(struct Hash_Node*) , tsize);
  hn = lu->htTop = (struct Hash_Node *)ZOLTAN_MALLOC(sizeof(struct Hash_Node) * ngids);

  if (tsize && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu);
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    return NULL;
  }

  lu->lenGID = lenGID;
  lu->table_size = tsize;
  lu->numGIDs = 0;

  for (i=0; i<ngids; i++, nextGID += lenGID){

    found = phg_lookup_GID(lu, nextGID);

    if (found >= 0) continue;

    hn->gid = nextUniqueGID;
    hn->gno = lu->numGIDs;

    if (nextUniqueGID < nextGID){
      for (k=0; k<lenGID; k++){
        nextUniqueGID[k] = nextGID[k]; 
      }
    }

    j = Zoltan_Hash(nextGID, lenGID, tsize);

    hn->next = lu->ht[j];
    lu->ht[j] = hn;
   
    hn++;
    nextUniqueGID += lenGID;
    lu->numGIDs++;
  }

  return lu;
}
/****************************************************************************/
int phg_lookup_GID(phg_GID_lookup *lu, ZOLTAN_ID_PTR gid)
{
  struct Hash_Node *hn;
  int i, k, match;

  if (lu->table_size < 1) return -1;
  if (lu->numGIDs < 1) return -1;

  i = Zoltan_Hash(gid, lu->lenGID, (unsigned int) lu->table_size);
  
  for (hn=lu->ht[i]; hn != NULL; hn = hn->next){
    match = 1;
    for (k=0; k<lu->lenGID; k++){
      if (hn->gid[k] != gid[k]){
        match = 0;
        break;
      }
    }
    if (match){
      return (hn->gno);
    }
  }
  return -1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
