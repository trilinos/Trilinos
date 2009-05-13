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

#include "zz_const.h"
#include "zz_util_const.h"
#include "phg_tree.h"
#include "phg.h"
#include <limits.h>


#define SET_MIN_NODE(ptr, offset, val) (ptr)[2*(offset)]=-(val)
#define SET_MAX_NODE(ptr, offset, val) (ptr)[2*(offset)+1]=(val)


int *
Zoltan_PHG_2ways_hyperedge_partition (
  ZZ *zz,                            /* Input : Zoltan data structure */
  HGraph *hg,
  Partition parts,
  Zoltan_PHG_Tree *tree
  )
{
  int nEdge, hEdge;
  int *interval;
  int *partnumber;
  int tree_size;
  int *rowpart;
  int index;
  int offset;

  nEdge = hg->nEdge;
  fprintf (stderr, "HG (%d %d x %d) : %d %d\n", hg->comm->myProc, hg->comm->myProc_x, hg->comm->myProc_y,  hg->nVtx, nEdge);

  interval = (int*)ZOLTAN_MALLOC(nEdge*2*sizeof(int));
  if ((nEdge > 0 ) && (interval == NULL))
    return NULL;

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
    ZOLTAN_FREE(&interval);
    return NULL;
  }

  rowpart = (int*) ZOLTAN_MALLOC(nEdge*sizeof(int));
  if ((nEdge > 0) && (rowpart == NULL)) {
    partnumber += 1;
    ZOLTAN_FREE(&partnumber);
    ZOLTAN_FREE(&interval);
    return (NULL);
  }

  offset = hg->dist_y[hg->comm->myProc_y];
  /* Then we search we is the hyperedge in the tree */
  for (hEdge = 0 ; hEdge < nEdge ; ++hEdge) {
    int node;
    node = find_interval_in_tree(tree, interval+2*hEdge);
    rowpart[hEdge] = partnumber[node];
    fprintf (stderr, "%d : %d (%d : %d - %d)\n", EDGE_LNO_TO_GNO(hg, hEdge), rowpart[hEdge], node, -interval[2*hEdge], interval[2*hEdge+1]);
  }

  partnumber += 1;
  ZOLTAN_FREE(&partnumber);
  ZOLTAN_FREE(&interval);

  return (rowpart);
}


#ifdef __cplusplus
}
#endif
