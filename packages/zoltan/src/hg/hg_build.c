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
#include "hypergraph.h"
#include "parmetis_jostle.h"

/*****************************************************************************/

int Zoltan_HG_Build_Hypergraph(
  ZZ *zz,
  struct Zoltan_HGraph **zoltan_hg,   /* Hypergraph to be allocated and built.*/
  int check_graph                     /* Parameter for hypergraph checking.   */
)
{
/* Input Zoltan Hypergraph from application */
char *yo = "Zoltan_HG_Build_Hypergraph";
struct Zoltan_HGraph *zhg;           /* Temporary pointer to Zoltan_HGraph. */
HGraph *hgraph;                      /* Temporary pointer to HG field */
int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Allocate a Zoltan hypergraph.  */
  zhg = *zoltan_hg = (struct Zoltan_HGraph *)
                      ZOLTAN_MALLOC(sizeof(struct Zoltan_HGraph));
  if (zhg == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* Initialize the Zoltan hypergraph data fields. */
  zhg->Global_IDs = NULL;
  zhg->Local_IDs = NULL;
  zhg->Parts = NULL;

  hgraph = &(zhg->HG);
  Zoltan_HG_HGraph_Init(hgraph);

  /* Use callback functions to build the hypergraph. */

  if (zz->Get_Num_HG_Edges != NULL && zz->Get_HG_Edge_List != NULL) {
    /* 
     *  Hypergraph callback functions exist; call them and build the HG
     *  directly.
     */

    ierr = Zoltan_Get_Obj_List(zz, &(hgraph->nVtx), 
                               &(zhg->Global_IDs), 
                               &(zhg->Local_IDs),
                               zz->Obj_Weight_Dim, &(hgraph->vwgt), 
                               &(zhg->Parts));
    Zoltan_HG_Fill_Hypergraph(zz, zhg);
  }

  else if (zz->Get_Num_Edges != NULL && zz->Get_Edge_List != NULL) {
    /* 
     *  Hypergraph callback functions don't exist, but graph functions do; 
     *  call the graph callback, build a graph, and convert it to a hypergraph.
     */
    Graph graph;             /* Temporary graph. */

    Zoltan_HG_Graph_Init(&graph);

    ierr = Zoltan_Get_Obj_List(zz, &(graph.nVtx), 
                               &(zhg->Global_IDs), 
                               &(zhg->Local_IDs),
                               zz->Obj_Weight_Dim, &(graph.vwgt),
                               &(zhg->Parts));
    Zoltan_HG_Fill_Hypergraph(zz, zhg);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      Zoltan_HG_Graph_Free(&graph);
      goto End;
    }

    ierr = Zoltan_Build_Graph(zz, 1, check_graph, graph.nVtx,
                       zhg->Global_IDs, zhg->Local_IDs, 
                       zz->Obj_Weight_Dim, zz->Edge_Weight_Dim,
                       &(graph.vtxdist), &(graph.nindex), &(graph.neigh), 
                       &(graph.ewgt));
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      Zoltan_HG_Graph_Free(&graph);
      goto End;
    }

    graph.nEdge = graph.nindex[graph.nVtx];
    ierr = Zoltan_HG_Graph_to_HGraph(zz, &graph, hgraph);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      Zoltan_HG_Graph_Free(&graph);
      goto End;
    }

    Zoltan_HG_Graph_Free(&graph);
  }

  if (zz->Get_Num_Geom != NULL && zz->Get_Geom != NULL) {
    /* 
     *  Geometric callbacks are registered; 
     *  use them to get coordinates for the hypergraph objects.
     */
    ierr = Zoltan_Get_Coordinates(zz, hgraph->nVtx,
                                  zhg->Global_IDs, 
                                  zhg->Local_IDs, &(hgraph->coor));
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error detected.");
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/
int Zoltan_HG_Fill_Hypergraph(
  ZZ *zz, 
  struct Zoltan_HGraph *zhg
)
{
/* KDD  -- Placeholder. ??? */


  return ZOLTAN_OK;
}

/*****************************************************************************/
int Zoltan_Get_Coordinates(
  ZZ *zz, 
  int num_obj,
  ZOLTAN_ID_PTR gids,
  ZOLTAN_ID_PTR lids,
  double *coor
)
{
/* KDD  -- Placeholder; should be in zz directory. ??? */


  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
