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

#include "hg.h"
#include "phg.h"
#include "parmetis_jostle.h"
#include "zz_util_const.h"

/*****************************************************************************/
/* Function prototypes */

static int hash_lookup (ZZ*, ZHG*, ZOLTAN_ID_PTR, int, struct Hash_Node**);
static int Zoltan_HG_Fill_Hypergraph (ZZ*, ZHG*, HGPartParams*);

/*****************************************************************************/

int Zoltan_HG_Build_Hypergraph(
  ZZ *zz,                            /* Zoltan data structure */
  ZHG **zoltan_hg,                   /* Hypergraph to be allocated and built.*/
  HGPartParams *hgp                  /* Parameters for HG partitioning.*/
)
{
/* allocates and builds hypergraph data structure using callback routines */
ZHG *zhg;                     /* Temporary pointer to Zoltan_HGraph. */
HGraph *hgraph;               /* Temporary pointer to HG field */
int ierr = ZOLTAN_OK;
char *yo = "Zoltan_HG_Build_Hypergraph";
int get_geom_data=0; /* Current hg methods don't use geometry. */

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Allocate a Zoltan hypergraph.  */
  zhg = *zoltan_hg = (ZHG*) ZOLTAN_MALLOC (sizeof(ZHG));
  if (zhg == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* Initialize the Zoltan hypergraph data fields. */
  zhg->nObj = 0;
  zhg->GIDs = NULL;
  zhg->LIDs = NULL;
  zhg->Input_Parts = NULL;
  zhg->Output_Parts = NULL;
  zhg->nRemove = 0;
  zhg->Remove_EGIDs = NULL;
  zhg->Remove_ELIDs = NULL;
  zhg->Remove_Esize = NULL;
  zhg->Remove_Ewgt = NULL;
  zhg->nRecv_GNOs = 0;
  zhg->Recv_GNOs = NULL;
  zhg->VtxPlan = NULL;

  hgraph = &(zhg->HG);
  Zoltan_HG_HGraph_Init(hgraph);
  hgraph->VtxWeightDim = zz->Obj_Weight_Dim;
  hgraph->EdgeWeightDim = zz->Edge_Weight_Dim;

  /* Use callback functions to build the hypergraph. */

  ierr = Zoltan_Get_Obj_List(zz, &(zhg->nObj), &(zhg->GIDs),
    &(zhg->LIDs), zz->Obj_Weight_Dim, &(hgraph->vwgt),
    &(zhg->Input_Parts));
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
    goto End;
  }
  hgraph->nVtx = zhg->nObj;

  ierr = Zoltan_HG_Fill_Hypergraph(zz, zhg, hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph");
    goto End;
  }

  /* Work-around for Trilinos/EpetraExt, which registers dummy query 
     functions even when they don't exist. */
  if (get_geom_data && (zz->Get_Num_Geom != NULL) && 
      (zz->Get_Geom != NULL || zz->Get_Geom_Multi != NULL)) {
     /* Geometric callbacks are registered;       */
     /* get coordinates for hypergraph objects.   */
     ZOLTAN_TRACE_DETAIL(zz, yo, "Getting Coordinates.");
     ierr = Zoltan_Get_Coordinates(zz, hgraph->nVtx, zhg->GIDs,
      zhg->LIDs, &(hgraph->nDim), &(hgraph->coor));
  }

  if (hgp->check_graph) {
    ierr = Zoltan_HG_Check(zz, hgraph);
    if (ierr == ZOLTAN_WARN) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning returned from Zoltan_HG_Check");
    }
    else if (ierr != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_HG_Check");
      goto End;     
    }
  }

  if (hgp->output_level >= HG_DEBUG_PRINT)
    Zoltan_HG_HGraph_Print(zz, zhg, &(zhg->HG), zhg->Input_Parts, stdout);

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    /* Return NULL zhg */
    Zoltan_HG_HGraph_Free(&(zhg->HG));
    Zoltan_Multifree(__FILE__, __LINE__, 4, &(zhg->GIDs),
     &(zhg->LIDs), &(zhg->Input_Parts), zoltan_hg);
  }
    
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/


static int Zoltan_HG_Fill_Hypergraph(
  ZZ *zz,
  ZHG *zhg,
  HGPartParams *hgp                  /* Parameters for HG partitioning.*/
)
{
/* Routine to call HG query function and build HG data structure.  */
/* Also builds Zoltan_HGraph vtxdist array.                        */

char *yo = "Zoltan_HG_Fill_Hypergraph";
ZOLTAN_ID_PTR edge_verts = NULL;  /* Object (vtx) GIDs belonging to hyperedges*/
int *edge_sizes = NULL;           /* # of GIDs in each hyperedge            */
int *edge_procs = NULL;           /* Processor owning each GID of hyperedge */
ZOLTAN_ID_PTR edge_gids = NULL, edge_lids = NULL;  /* Edge GID and LID */

struct Hash_Node *hash_nodes = NULL;  /* Hash table variables for mapping   */
struct Hash_Node **hash_tab = NULL;   /* GIDs to global numbering system.   */

int i, j;
int cnt;
int ierr = ZOLTAN_OK;

ZOLTAN_ID_PTR global_ids = zhg->GIDs;  
HGraph *hg = &(zhg->HG);
int nVtx = hg->nVtx;                     
int num_gid_entries = zz->Num_GID;
static PHGComm scomm;
static int first_time = 1;

  if (zz->Get_Num_HG_Edges && zz->Get_HG_Edge_List && zz->Get_HG_Edge_Info)
    ierr = Zoltan_HG_Hypergraph_Callbacks(zz, zhg, nVtx, hgp->EdgeSizeThreshold,
                                          1, &(hg->nEdge), 
                                          &edge_gids, &edge_lids, &edge_sizes,
                                          &(hg->ewgt), &(hg->nPins), 
                                          &edge_verts, &edge_procs);
  else if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
           (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) 
    ierr = Zoltan_HG_Graph_Callbacks(zz, zhg, nVtx, hgp->EdgeSizeThreshold,
                                     1, &(hg->nEdge), 
                                     &edge_gids, &edge_lids, &edge_sizes,
                                     &(hg->ewgt), &(hg->nPins), 
                                     &edge_verts, &edge_procs);
  else {
    hg->nEdge = 0;
    hg->nPins = 0;
  }

  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Callbacks");
    goto End;
  }


  /* Build hg->hindex */
  /* KDD -- should we remove HEdges with size 1 from edge lists here? */
  hg->hindex = (int *) ZOLTAN_MALLOC((hg->nEdge + 1) * sizeof(int));
  if (!hg->hindex) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  cnt = 0;
  for (i = 0; i < hg->nEdge; i++) {
    hg->hindex[i] = cnt;
    cnt += edge_sizes[i];
  }
  hg->hindex[hg->nEdge] = cnt;

  /* Sanity check */
  if (cnt != hg->nPins) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Input error:  Number of pins != sum of edge sizes");
    goto End;
  }

  /*  Set up place-holder arrays and communicators that, while not used 
   *  in serial, are needed by some of the utility functions shared by 
   *  HG and PHG.
   */
  hg->dist_x = (int *)ZOLTAN_MALLOC(2 * sizeof(int));
  hg->dist_y = (int *)ZOLTAN_MALLOC(2 * sizeof(int));
  if (!hg->dist_x || !hg->dist_y){
    /* Not enough memory */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    goto End;
  }

  hg->dist_x[0] = 0;
  hg->dist_x[1] = hg->nVtx;
  hg->dist_y[0] = 0;
  hg->dist_y[1] = hg->nEdge;

  if (first_time) {
    scomm.nProc_x = scomm.nProc_y = 1;
    scomm.myProc_x = scomm.myProc_y = 0;
    scomm.Communicator = MPI_COMM_SELF;
    scomm.row_comm = MPI_COMM_SELF;
    scomm.col_comm = MPI_COMM_SELF;
    scomm.myProc = 0;
    scomm.nProc = 1;
    scomm.RNGState = Zoltan_Rand(NULL);
    scomm.RNGState_row = Zoltan_Rand(NULL);
    scomm.RNGState_col = Zoltan_Rand(NULL);
    first_time = 0;
  }
  scomm.zz = zz;

  hg->comm = &scomm;
  
  if (hg->nPins > 0) {
    /* 
     * Correlate GIDs in edge_verts with local indexing in zhg to build the
     * input HG.
     * Use hash table to map global IDs to local position in zhg->GIDs.
     * Based on hashing code in Zoltan_Build_Graph.
     * KDD -- This approach is serial for now; look more closely at 
     * KDD -- Zoltan_Build_Graph when move to parallel.
     */

    /* Construct local hash table */
    hash_nodes = (struct Hash_Node *)ZOLTAN_MALLOC(nVtx * 
                                                   sizeof(struct Hash_Node));
    hash_tab = (struct Hash_Node **) ZOLTAN_MALLOC(nVtx *
                                                   sizeof(struct Hash_Node *));
    if (nVtx && ((!hash_nodes) || (!hash_tab))){
      /* Not enough memory */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      goto End;
    }

    /* Assign consecutive numbers based on the order of the ids */
    for (i=0; i< nVtx; i++) {
      hash_tab[i] = NULL;
      hash_nodes[i].gid = &(global_ids[i*num_gid_entries]);
      hash_nodes[i].gno = hg->dist_x[zz->Proc]+i;
    }

    for (i=0; i< nVtx; i++){
      /* insert hashed elements into hash table */
      j = Zoltan_Hash(&(global_ids[i*num_gid_entries]), num_gid_entries,
                      (unsigned int) nVtx);
      hash_nodes[i].next = hash_tab[j];
      hash_tab[j] = &hash_nodes[i];
    }

    hg->hvertex = (int *) ZOLTAN_MALLOC(hg->nPins * sizeof(int));
    if (!hg->hvertex) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < hg->nPins; i++) {
      hg->hvertex[i] = hash_lookup(zz, zhg, &(edge_verts[i*num_gid_entries]),
                                   nVtx, hash_tab);
      if (hg->hvertex[i] == -1) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Hyperedge GID not found.")
        ierr = ZOLTAN_FATAL;
        goto End;
      }
    }

    Zoltan_Multifree(__FILE__, __LINE__, 2, &hash_nodes,
                                            &hash_tab);
  }
  Zoltan_Multifree(__FILE__, __LINE__, 5, &edge_gids, 
                                          &edge_lids, 
                                          &edge_verts, 
                                          &edge_sizes, 
                                          &edge_procs);

  ierr = Zoltan_HG_Create_Mirror(zz, hg);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error from Zoltan_HG_Create_Mirror");
    goto End;
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    Zoltan_HG_HGraph_Free(hg);
  
    Zoltan_Multifree(__FILE__, __LINE__, 7, &edge_gids,
                                            &edge_lids,
                                            &edge_verts, 
                                            &edge_sizes, 
                                            &edge_procs, 
                                            &hash_nodes,
                                            &hash_tab);
  }
  return ierr;
}

/*****************************************************************************/

static int hash_lookup(
  ZZ *zz,
  ZHG *zhg,
  ZOLTAN_ID_PTR key,
  int nVtx,
  struct Hash_Node **hash_tab
)
{
/* Looks up a key GID in the hash table; returns its gno. */
/* Based on hash_lookup in build_graph.c. */

  int i;
  struct Hash_Node *ptr;

  i = Zoltan_Hash(key, zz->Num_GID, (unsigned int) nVtx);
  for (ptr=hash_tab[i]; ptr != NULL; ptr = ptr->next){
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->gno);
  }
  /* Key not in hash table */
  return -1;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
