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

#include "phg.h"
#include "phg_verbose.h"
#include "phg_lookup.h"

/****************************************************************************/

void print_zoltan_pins(zoltan_pins *z, int me, int ewgt_dim)
{
  int i;

  printf("%d) %d hyperedges\n\n",me, z->nHedges);

  if (z->nHedges == 0) return;

  for (i=0; i<z->nHedges; i++){
    if (z->edgeHash){
      printf("  GID " ZOLTAN_ID_SPEC ", hashed to %d, num pins %d\n", z->edgeGID[i], z->edgeHash[i], z->esizes[i]);
    }
    else{
      printf("  GID " ZOLTAN_ID_SPEC ", num pins locally %d\n", z->edgeGID[i], z->esizes[i]);
    }
  }
  printf("\n");
}

/****************************************************************************/
void print_hypergraph(ZZ *zz, ZHG *zhg, int sumWeight)
{
  int i, j, npins;
  int ewdim = zhg->edgeWeightDim;
  int vwdim = zhg->objWeightDim;
  float sum;
  float *wgt, *vwgt;
  int *owner;
  ZOLTAN_GNO_TYPE *pin;
  int *lno;
  HGraph *hg = &zhg->HG;
  int p = zz->Proc;

  /* The ZHG structure contains the hypergraph returned by the query functions,
   * including modifications based on ADD_OBJ_WEIGHT and PHG_EDGE_WEIGHT_OPERATION.
   * If the PHG hypergraph build has completed, the edge list only contains the removed
   * edges.  If LB_Eval build the ZHG structure, it contains all edges.
   *
   * the HGraph structure contains that hypergraph with modifications made
   * for the PHG algorithm.  This may include addition of repartition
   * vertices and edges, and removal of dense edges.
   */

  wgt = zhg->objWeight;

  printf("(%d) %d INPUT VERTICES (out of " ZOLTAN_GNO_SPEC ") : gno (gid/lid) (weights) nhedges fixed inpart outpart objSize)\n",p, zhg->nObj, zhg->globalObj);

  for (i=0; i<zhg->nObj; i++){

    printf("  " ZOLTAN_GNO_SPEC " (",zhg->objGNO[i]);

    if (zhg->objGID)
      printf(ZOLTAN_ID_SPEC "/",zhg->objGID[i]);
    else
      printf("-/");

    if (zhg->objLID)
      printf(ZOLTAN_ID_SPEC ") (",zhg->objLID[i]);
    else
      printf("/-) (");

    for (j=0; j < vwdim; j++){
      printf("%f",*wgt++);
      if (j < vwdim-1) printf(", ");
    }

    if (zhg->numHEdges)
      printf(") %d ",zhg->numHEdges[i]);
    else
      printf(") - ");

    if (zhg->fixed)
      printf(" %d ",zhg->fixed[i]);
    else
      printf(" - ");

    if (zhg->Input_Parts)
      printf(" %d ",zhg->Input_Parts[i]);
    else
      printf(" - ");

    if (zhg->Output_Parts)
      printf(" %d ",zhg->Output_Parts[i]);
    else
      printf(" - ");

    if (zhg->AppObjSizes)
      printf(" %d ",zhg->AppObjSizes[i]);
    else
      printf(" - ");

    printf("\n");
  }
  printf("\n");
 
  wgt = zhg->Ewgt;
  pin = zhg->pinGNO;
  owner = zhg->Pin_Procs;
   
  printf("(%d) %d INPUT or REMOVED EDGES (out of " ZOLTAN_GNO_SPEC "), %d pins: gno size (weights) (pinGNO/pinProc)\n",
                  p, zhg->nHedges, zhg->globalHedges, zhg->nPins);

  for (i=0; i < zhg->nHedges; i++){

    printf("  " ZOLTAN_GNO_SPEC " %d (", zhg->edgeGNO[i], zhg->Esize[i]);

    if (wgt){
      for (j=0; j < ewdim; j++){
        printf("%f",*wgt++);
        if (j < ewdim - 1) printf(", ");
      }
    }
    printf(") (");

    for (j=0; j < zhg->Esize[i]; j++){
      printf("" ZOLTAN_GNO_SPEC "/%d", *pin++, *owner++);
      if (j < zhg->Esize[i] - 1) printf(" ");
    }

    printf(")\n");
  }
  printf("\n");

  printf("(%d) %d PHG EDGES (%d weights), %d total PHG PINS:\n",
          p, hg->nEdge, ewdim, hg->nPins);

  wgt = hg->ewgt;
  lno = hg->hvertex;
  vwgt = hg->vwgt;

  for (i=0; i<hg->nEdge; i++){
    npins = hg->hindex[i+1] - hg->hindex[i];

    printf(" edge " ZOLTAN_GNO_SPEC ": ",EDGE_LNO_TO_GNO(hg, i));
    for (j=0; j<ewdim; j++){
      printf(" %f",*wgt++);
    }
    printf("\n %d pins: ", npins);
    for (j=0; j<npins; j++){
      printf("%d ", *lno++);
    }
    printf("\n");
  }
  printf("\n");

  printf("(%d) %d PHG PIN global numbers and %d weights:\n", p, hg->nVtx, vwdim);

  sum = 0;

  for (i=0; i<hg->nVtx; i++){
    printf("  %d  " ZOLTAN_GNO_SPEC ": ", i, VTX_LNO_TO_GNO(hg, i));
    for (j=0; j<vwdim; j++){
      if (j==sumWeight) sum += *vwgt;
      printf("%f ", *vwgt++);
    }
    printf("\n");
  }
  printf("\n");
  if (sum > 0.0) printf("(%d) Weight %d sums to %f\n\n",p, sumWeight+1,sum);
}
/****************************************************************************/
void show_edges(char *s, ZZ *zz, int num_lists, int num_pins,
                ZOLTAN_ID_TYPE *edg_GID, int *row_ptr, ZOLTAN_ID_TYPE *vtx_GID)
{
int i, j, size, sumsize=0;
ZOLTAN_ID_TYPE *v = vtx_GID;

  /* helpful in debugging */
  printf("%s> Process %d, %d edges, %d pins\n",s, zz->Proc, num_lists, num_pins);
  for (i=0; i<num_lists; i++){
    size = (i < num_lists-1 ? row_ptr[i+1] : num_pins) - row_ptr[i];
    sumsize += size;
    printf("Edge " ZOLTAN_ID_SPEC ", size %d\n  ", edg_GID[i], size);
    for (j=0; j<size; j++){
      printf(ZOLTAN_ID_SPEC " ",   *v++);
    }
    printf("\n");
  }
  printf("Sum of edge sizes: %d\n",sumsize);
}
/****************************************************************************/
void debug_graph_to_hg(
  int nedges, ZOLTAN_ID_PTR egids, ZOLTAN_ID_PTR elids,
  int *esizes, float *ewgts, int npins,
  ZOLTAN_ID_PTR pins, int *pin_procs, int ewgtdim, int lenGID, int lenLID)
{
  int i,j,k;
  ZOLTAN_ID_PTR nextpin;
  int *nextproc;

  nextpin = pins;
  nextproc = pin_procs;

  printf("%d hyperedges, %d pins\n",nedges,npins);
  for (i=0; i<nedges; i++){
    printf("GID ");
    for (j=0; j<lenGID; j++) printf(ZOLTAN_ID_SPEC " ", egids[i*lenGID+ j]);
    printf(" LID ");
    for (j=0; j<lenLID; j++) printf(ZOLTAN_ID_SPEC " ", elids[i*lenLID+ j]);
    printf(" weights ");
    for (j=0; j<ewgtdim; j++) printf("%f ", ewgts[i*ewgtdim+ j]);
    printf(" size %d\n",esizes[i]);

    for (j=0; j < esizes[i]; j++){
      printf("  ");
      for (k=0; k<lenGID; k++) printf(ZOLTAN_ID_SPEC " ", *nextpin++);
      printf(" (%d), ",*nextproc++);
      if (j && (j%10==0)) printf("\n");
    }
    printf("\n");
  }
}
/****************************************************************************/




#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
