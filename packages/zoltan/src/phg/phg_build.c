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

#include "phg.h"
#include "phypergraph.h"
#include "parmetis_jostle.h"
#include "zz_util_const.h"



/*****************************************************************************/
/* Function prototypes */

static int hash_lookup (ZZ*, ZHG*, ZOLTAN_ID_PTR, int, struct Hash_Node**);
static int Zoltan_PHG_Fill_Hypergraph (ZZ*, ZHG*);
/*****************************************************************************/



/* allocates and builds hypergraph data structure using callback routines */
int Zoltan_PHG_Build_Hypergraph(
  ZZ *zz,                            /* Zoltan data structure */
  ZHG **zoltan_hg,                   /* Hypergraph to be allocated and built */
  PHGPartParams *hgp                 /* Parameters for HG partitioning */
)
{
  ZHG *zhg;                         /* Temporary pointer to Zoltan_HGraph */
  PHGraph *hgraph;                  /* Temporary pointer to PHG field */
  int err = ZOLTAN_OK;
  char *yo = "Zoltan_PHG_Build_Hypergraph";

  ZOLTAN_TRACE_ENTER (zz, yo);

  /* Allocate a Zoltan hypergraph.  */
  zhg = *zoltan_hg = (ZHG*) ZOLTAN_MALLOC (sizeof(ZHG));
  if (zhg == NULL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
    err = ZOLTAN_MEMERR;
    goto End;
  }

  /* Initialize the Zoltan hypergraph data fields. */
  zhg->Global_IDs = NULL;
  zhg->Local_IDs  = NULL;
  zhg->Parts      = NULL;

  hgraph = &zhg->HG;
  Zoltan_PHG_HGraph_Init(hgraph);

  /* Use callback functions to build the hypergraph. */
  if (zz->Get_Num_HG_Edges && zz->Get_HG_Edge_List && zz->Get_Num_HG_Pins) {
    /* Hypergraph callback functions exist; call them and build the PHG */
    ZOLTAN_TRACE_DETAIL(zz, yo, "Using Parallel Hypergraph Callbacks.");

    err = Zoltan_Get_Obj_List(zz, &hgraph->nVtx, &zhg->Global_IDs,
     &zhg->Local_IDs, zz->Obj_Weight_Dim, &hgraph->vwgt, &zhg->Parts);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
      goto End;
    }

    err = Zoltan_PHG_Fill_Hypergraph(zz, zhg);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph");
      goto End;
    }
  }

  else if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL)
        && (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
        
    /* Hypergraph callback functions don't exist, but graph functions do;     */
    /* call the graph callback, build a graph, and convert it to a hypergraph */
    PGraph graph;             /* Temporary graph. */

    ZOLTAN_TRACE_DETAIL(zz, yo, "Using Graph Callbacks.");
    Zoltan_PHG_Graph_Init (&graph);
    err = Zoltan_Get_Obj_List(zz, &graph.nVtx, &zhg->Global_IDs,
     &zhg->Local_IDs, zz->Obj_Weight_Dim, &graph.vwgt, &zhg->Parts);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
      Zoltan_PHG_Graph_Free (&graph);
      goto End;
    }

    err = Zoltan_Build_Graph (zz, 1, hgp->check_graph, graph.nVtx,
     zhg->Global_IDs, zhg->Local_IDs, zz->Obj_Weight_Dim, zz->Edge_Weight_Dim,
     &graph.vtxdist, &graph.nindex, &graph.neigh, &graph.ewgt);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building graph");
      Zoltan_PHG_Graph_Free (&graph);
      goto End;
    }

    graph.nEdge = graph.nindex[graph.nVtx];
    err = Zoltan_PHG_Graph_to_HGraph (zz, &graph, hgraph);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error converting graph to hypergraph");
      Zoltan_PHG_Graph_Free(&graph);
      goto End;
    }
    Zoltan_PHG_Graph_Free(&graph);
  }

  if (zz->Get_Num_Geom != NULL && 
      (zz->Get_Geom != NULL || zz->Get_Geom_Multi != NULL)) {
      
     /* Geometric callbacks are registered;       */
     /* get coordinates for hypergraph objects.   */
     ZOLTAN_TRACE_DETAIL(zz, yo, "Getting Coordinates.");
     err = Zoltan_Get_Coordinates (zz, hgraph->nVtx, zhg->Global_IDs,
      zhg->Local_IDs, &hgraph->nDim, &hgraph->coor);
  }

  if (hgp->check_graph) {
    err = Zoltan_PHG_Check(zz, hgraph);
    if (err == ZOLTAN_WARN) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning returned from Zoltan_PHG_Check");
    }
    else if (err != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_PHG_Check");
      goto End;     
    }
  }

  if (hgp->output_level >= PHG_DEBUG_PRINT)
    Zoltan_PHG_HGraph_Print(zz, zhg, &zhg->HG, stdout);

End:
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    /* Return NULL zhg */
    Zoltan_PHG_HGraph_Free (&zhg->HG);
    Zoltan_Multifree(__FILE__, __LINE__, 4, &zhg->Global_IDs, &zhg->Local_IDs,
     &zhg->Parts, zoltan_hg);
  }
    
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}



/*****************************************************************************/
/* Routine to call HG query function and build HG data structure.  */
/* Also builds Zoltan_HGraph vtxdist array.                        */

static int Zoltan_PHG_Fill_Hypergraph(
  ZZ *zz,
  ZHG *zhg
)
{
  ZOLTAN_ID_PTR edge_verts = NULL;  /* Object GIDs belonging to hyperedges    */
  int *edge_sizes  = NULL;          /* # of GIDs in each hyperedge            */
  int *edge_procs  = NULL;          /* Processor owning each GID of hyperedge */
  float *edge_wgts = NULL;          /* Hyperedge weights                      */

  struct Hash_Node *hash_nodes = NULL;  /* Hash table variables for mapping   */
  struct Hash_Node **hash_tab = NULL;   /* GIDs to global numbering system.   */

  int i, j;
  int npins, cnt;
  int numwgts = 0;
  int err = ZOLTAN_OK;

  ZOLTAN_ID_PTR global_ids = zhg->Global_IDs;  
  PHGraph *hg = &zhg->HG;
  int nVtx = hg->nVtx;                     
  int num_gid_entries = zz->Num_GID;
  char *yo = "Zoltan_PHG_Fill_Hypergraph";
  
  /* Build vtxdist as in Zoltan_Build_Graph. */
  /* KDD -- I am guessing we will need this array in parallel; we may not. */
  hg->vtxdist = (int*) ZOLTAN_MALLOC ((zz->Num_Proc+1) * sizeof(int));
  if (!(hg->vtxdist)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");   /* Not enough memory */
    goto End;
  }

  /* Construct vtxdist[i] = the number of vertices on all procs < i. */
  /* Scan to compute partial sums of the number of objs */
  MPI_Scan (&nVtx, hg->vtxdist, 1, MPI_INT, MPI_SUM, zz->Communicator);

  /* Gather data from all procs */
  MPI_Allgather (&hg->vtxdist[0], 1, MPI_INT, &hg->vtxdist[1], 1, MPI_INT,
   zz->Communicator);
  hg->vtxdist[0] = 0;

  /* Get hyperedge information from application through query functions. */
  hg->nEdge = zz->Get_Num_HG_Edges (zz->Get_Num_HG_Edges_Data, &err);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_HG_Edges");
    goto End;
  }
  
  /* KDD:  question:  How do we compute size to malloc array for HG Edges? 
   * KDD:  We can't have a size function unless we assume the application
   * KDD:  can "name" the hyperedges.
   * KDD:  For now, assume application can return number of pins.
   */

  hg->nInput = npins = zz->Get_Num_HG_Pins(zz->Get_Num_HG_Pins_Data, &err);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Max_PHG_Edge_Size");
    goto End;
  }

  if (hg->nEdge > 0) {
    numwgts = hg->nEdge * zz->Edge_Weight_Dim;
    edge_verts = ZOLTAN_MALLOC_GID_ARRAY(zz, npins);
    edge_sizes = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int));
    edge_procs = (int*) ZOLTAN_MALLOC (npins     * sizeof(int));
    if (numwgts) edge_wgts = (float*) ZOLTAN_MALLOC (numwgts * sizeof(float));
    if (!edge_verts || !edge_sizes || !edge_procs || (numwgts && !edge_wgts)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      err = ZOLTAN_MEMERR;
      goto End;
    }
    
    err = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data, num_gid_entries,
     zz->Edge_Weight_Dim, hg->nEdge, npins, edge_sizes, edge_verts, edge_procs,
     edge_wgts);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_PHG_Edge_List");
      goto End;
    }
  }
  
  /* Build hg->hindex */
  hg->hindex = (int *) ZOLTAN_MALLOC((hg->nEdge + 1) * sizeof(int));
  if (!hg->hindex) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
    err = ZOLTAN_MEMERR;
    goto End;
  }
  cnt = 0;
  for (i = 0; i < hg->nEdge; i++) {
    hg->hindex[i] = cnt;
    cnt += edge_sizes[i];
  }
  hg->hindex[hg->nEdge] = cnt;

  /* Sanity check */
  if (cnt != npins) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
     "Input error: Number of pins != sum of edge sizes");
    goto End;
  }

  /* Correlate GIDs in edge_verts with local indexing in zhg to build the input */
  /* HG.  Use hash table to map global IDs to local position in zhg->Global_IDs */
  /* Based on hashing code in Zoltan_Build_Graph                                */
  /* KDD -- This approach is serial for now; look more closely at               */
  /* KDD -- Zoltan_Build_Graph when move to parallel.                           */

  if (npins > 0) {
    /* Construct local hash table */
    hash_nodes=(struct Hash_Node*) ZOLTAN_MALLOC(nVtx*sizeof(struct Hash_Node));
    hash_tab = (struct Hash_Node**)ZOLTAN_MALLOC(nVtx*sizeof(struct Hash_Node*));
    if (nVtx && ((!hash_nodes) || (!hash_tab))) {  
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");  /* Not enough memory */
      goto End;
    }

    /* Assign consecutive numbers based on the order of the ids */
    for (i=0; i< nVtx; i++) {
      hash_tab[i] = NULL;
      hash_nodes[i].gid = &global_ids[i*num_gid_entries];
      hash_nodes[i].gno = hg->vtxdist[zz->Proc]+i;
    }

    for (i=0; i< nVtx; i++){
      /* insert hashed elements into hash table */
      j = Zoltan_Hash(&global_ids[i*num_gid_entries], num_gid_entries,
       (unsigned int) nVtx);
      hash_nodes[i].next = hash_tab[j];
      hash_tab[j] = &hash_nodes[i];
    }

    hg->hvertex = (int*) ZOLTAN_MALLOC (npins * sizeof(int));
    if (!hg->hvertex) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      err = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < npins; i++) {
      hg->hvertex[i] = hash_lookup(zz, zhg, &edge_verts[i], nVtx, hash_tab);
      if (hg->hvertex[i] == -1) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Hyperedge GID not found.")
        err = ZOLTAN_FATAL;
        goto End;
      }
    }

    if (zz->Edge_Weight_Dim) {
      hg->EdgeWeightDim = zz->Edge_Weight_Dim;
      hg->ewgt = (float*) ZOLTAN_MALLOC (numwgts * sizeof(float));
      memcpy(hg->ewgt, edge_wgts, numwgts * sizeof(float));
    }

    Zoltan_Multifree(__FILE__, __LINE__, 2, &hash_nodes, &hash_tab);
  }
  Zoltan_Multifree(__FILE__, __LINE__, 4, &edge_verts, &edge_sizes, &edge_procs, 
   &edge_wgts);

  err = Zoltan_PHG_Create_Mirror (zz, hg);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Error from Zoltan_PHG_Create_Mirror");
    goto End;
  }

End:
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    Zoltan_PHG_HGraph_Free(hg);
    Zoltan_Multifree(__FILE__, __LINE__, 6, &edge_verts, &edge_sizes, 
     &edge_procs, &edge_wgts, &hash_nodes, &hash_tab);
  }
  return err;
}



/*****************************************************************************/
/* Looks up a key GID in the hash table; returns its gno. */
/* Based on hash_lookup in build_graph.c. */

static int hash_lookup(
  ZZ *zz,
  ZHG *zhg,
  ZOLTAN_ID_PTR key,
  int nVtx,
  struct Hash_Node **hash_tab
)
{
  int i;
  struct Hash_Node *ptr;

  i = Zoltan_Hash (key, zz->Num_GID, (unsigned int) nVtx);
  for (ptr = hash_tab[i]; ptr != NULL; ptr = ptr->next) {
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->gno);
  }
  return -1;   /* Key not in hash table */
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
