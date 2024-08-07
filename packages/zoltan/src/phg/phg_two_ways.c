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

#include "zz_const.h"
#include "zz_util_const.h"
#include "phg_tree.h"
#include "phg.h"
#include <limits.h>

/* #define CEDRIC_PRINT */

#define SET_MIN_NODE(ptr, offset, val) (ptr)[2*(offset)]=-(val)
#define SET_MAX_NODE(ptr, offset, val) (ptr)[2*(offset)+1]=(val)


int
Zoltan_PHG_2ways_hyperedge_partition (
  ZZ *zz,                            /* Input : Zoltan data structure */
  HGraph *hg,
  Partition parts,
  Zoltan_PHG_Tree *tree,
  struct Zoltan_DD_Struct * gnoToGID,
  struct Zoltan_DD_Struct **dd,
  int *numParts,
  int **sizeParts
)
{
  int ierr = ZOLTAN_OK;
  char *yo = "Zoltan_PHG_2ways_hyperedge_partition";
  int nEdge, hEdge;
  int *interval;
  int *partnumber = NULL;
  int tree_size;
  ZOLTAN_ID_TYPE *rowpart =NULL;  /* ZOLTAN_ID_TYPE because it's used in Zoltan_DD_* */
  ZOLTAN_GNO_TYPE *rowGNO = NULL;
  ZOLTAN_ID_PTR rowGID=NULL;
  int index;
  int offset;

  ZOLTAN_TRACE_ENTER(zz, yo);

  nEdge = hg->nEdge;
  fprintf (stderr, "HG (%d %d x %d) : %d %d\n", hg->comm->myProc, hg->comm->myProc_x, hg->comm->myProc_y,  hg->nVtx, nEdge);

  interval = (int*)ZOLTAN_MALLOC(nEdge*2*sizeof(int));
  if ((nEdge > 0 ) && (interval == NULL)) MEMORY_ERROR;

  tree_size = get_tree_size(tree) + 1;
  for (index = 0 ; index < nEdge ; ++index){
    SET_MIN_NODE(interval, index, tree_size);
    SET_MAX_NODE(interval, index, -1);
  }

  /* Update interval with the local knowledge */
  /* XXX: I loop on the hyperedges, as I think it's more cache friendly
   * and it allows me to discoupled the computation if a non blocking MPI_Reduce is
   * available
   */
  for (hEdge = 0 ; hEdge < nEdge ; ++hEdge){
    int part;
    int max = -1;                     /* Trick : we use the initialized values */
    int min = tree_size + 1;

    for (index = hg->hindex[hEdge] ; index < hg->hindex[hEdge+1] ; ++ index) {
      part = parts[hg->hvertex[index]];

      max = MAX(max, part);
      min = MIN(min, part);
    }
    SET_MIN_NODE(interval, hEdge, min);
    SET_MAX_NODE(interval, hEdge, max);
  }

  /* Update results to view the complete hyperedges */
  Zoltan_AllReduceInPlace (interval, 2*nEdge, MPI_INT, MPI_MAX, hg->comm->row_comm);

  /* Now I have to compute the partition of hyperedges according to the "interval"
   * and the tree */

  /* First, compute the partition number corresponding to the nodes in the tree */
  partnumber = compute_part_number(tree);
  if (partnumber == NULL) {
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  (*numParts) = get_tree_size(tree);

  rowpart = (ZOLTAN_ID_TYPE*) ZOLTAN_MALLOC(nEdge*sizeof(ZOLTAN_ID_TYPE));
  if ((nEdge > 0) && (rowpart == NULL)) MEMORY_ERROR;

  rowGNO = (ZOLTAN_GNO_TYPE*) ZOLTAN_MALLOC(nEdge*sizeof(ZOLTAN_GNO_TYPE));
  if ((nEdge > 0) && (rowGNO == NULL)) MEMORY_ERROR;

  (*sizeParts) = (int*)ZOLTAN_CALLOC((*numParts), sizeof(int));
  if (*numParts && (*sizeParts) == NULL) MEMORY_ERROR;

  offset = hg->dist_y[hg->comm->myProc_y];
  (void) offset; /* forestall warning for "variable set but unused" */
  /* Then we search we is the hyperedge in the tree */
  for (hEdge = 0 ; hEdge < nEdge ; ++hEdge) {
    int node;
    node = find_interval_in_tree(tree, interval+2*hEdge);
    rowpart[hEdge] = partnumber[node];
    (*sizeParts)[rowpart[hEdge]] ++;
    rowGNO[hEdge] = EDGE_LNO_TO_GNO(hg, hEdge);
#if 0
    fprintf (stderr, "%zd : " ZOLTAN_ID_SPEC " (%d : %d - %d)\n", rowGNO[hEdge], rowpart[hEdge], node, -interval[2*hEdge], interval[2*hEdge+1]);
#endif
  }

  partnumber += 1;
  ZOLTAN_FREE(&partnumber);
  ZOLTAN_FREE(&interval);

  /* Compute number of elements per parts */
  /* TODO: support processor which are not part of the distribution */

  /* Update results to view the complete hyperedges */
  Zoltan_AllReduceInPlace ((*sizeParts), (*numParts), MPI_INT, MPI_SUM, hg->comm->col_comm);


  /* Export results to data directory */
  /* First, get the GIDs of our edges */
  rowGID = ZOLTAN_MALLOC_GID_ARRAY(zz, nEdge);
  if (nEdge && rowGID == NULL) MEMORY_ERROR;
  ierr = Zoltan_DD_Find (gnoToGID , (ZOLTAN_ID_PTR)rowGNO, rowGID, NULL, NULL,
                         nEdge, NULL);
  ZOLTAN_FREE(&rowGNO);

  ierr = Zoltan_DD_Create (dd, zz->Communicator, zz->Num_GID, 1, 0, nEdge, 0);
  CHECK_IERR;

  /* Make our new numbering public */
  Zoltan_DD_Update (*dd, (ZOLTAN_ID_PTR)rowGID, rowpart, NULL,  NULL, nEdge);

#ifdef CEDRIC_PRINT
  for (hEdge = 0 ; hEdge < nEdge ; ++hEdge) {
    fprintf (stderr, "%d : %d\n", rowGID[hEdge], rowpart[hEdge]);
  }
#endif


 End:
  ZOLTAN_FREE(&rowGID);
  ZOLTAN_FREE(&rowGNO);
  ZOLTAN_FREE(&rowpart);

  if (partnumber != NULL)
    partnumber += 1;
  ZOLTAN_FREE(&partnumber);
  ZOLTAN_FREE(&interval);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}


#ifdef __cplusplus
}
#endif
