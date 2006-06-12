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
#include "zz_const.h"
#include "zz_util_const.h"

    
#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}
#define FATAL_ERROR(s) { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, s); \
  ierr = ZOLTAN_FATAL; \
  goto End; \
}



/*****************************************************************************/
int Zoltan_PHG_Build_Finegrain(
  ZZ *zz,                      /* In: Zoltan data structure */
  HGraph *hg,                  /* In: Hypergraph */
  HGraph **fg_hg,              /* Out: fine-grain hypergraph */
  PHGPartParams *hgp           /* Parameters for PHG partitioning.*/
)
{
/* allocates and builds fine-grain hypergraph from hg. */
int ierr = ZOLTAN_OK;
int i, j, k;
HGraph *fg;
char *yo = "Zoltan_PHG_Build_Finegrain";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Allocate new hgraph. Initialize relevant fields. */
  fg = *fg_hg = (HGraph *) ZOLTAN_MALLOC(sizeof(HGraph));
  Zoltan_HG_HGraph_Init(fg);
  fg->comm = hg->comm;
  fg->nVtx = hg->nPins;
  fg->nEdge = hg->nVtx + hg->nEdge;
  fg->nPins = 2*hg->nPins;
  fg->VtxWeightDim  = 0;
  fg->EdgeWeightDim = 0;

  /* Allocate space for vertex-based lookup in new fine grain hgraph.  */
  fg->vindex = (int *) ZOLTAN_MALLOC((fg->nVtx+1)*sizeof(int));
  fg->vedge = (int *) ZOLTAN_MALLOC((fg->nPins)*sizeof(int));

  /* Each proc first builds part of the fine-grain hgraph from local data. */

  /* Loop over the pins in the old hgraph by rows (hyperedge) */
  /* Each hyperedge (row) in the old hg becomes an hyperedge */
  /* Each vertex (column) in the old hg becomes an hyperedge */
  k = 0;
  for (i=0; i<hg->nEdge; i++){
    for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++){
      fg->vindex[k] = 2*k;
      fg->vedge[2*k] = i;
      fg->vedge[2*k+1] = hg->nEdge+hg->hvertex[j]; 
      k++;
    }
  }
  fg->vindex[k] = 2*k; /* k = fg->nVtx */

#ifdef DEBUG
  printf("Debug: fine hgraph (vertex-based). nvtx=%d, nedge=%d, nPins=%d\n", fg->nVtx, fg->nEdge, fg->nPins);
  printf("Debug: vindex: ");
  for (i=0; i<=fg->nVtx; i++)
    printf(" %d ", fg->vindex[i]);
  printf("\n");
  printf("Debug: vedge: ");
  for (i=0; i<fg->nPins; i++)
    printf(" %d ", fg->vedge[i]);
  printf("\n");
#endif

  /* Mirror hypergraph. This sets up hindex. */
  ierr = Zoltan_HG_Create_Mirror(zz, fg); 

  /* TODO: Combine local partial hgraphs to a global hgraph. */

End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

#ifdef __cplusplus
} 
#endif

