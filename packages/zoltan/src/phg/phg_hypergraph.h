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

#ifndef __PHG_HYPERGRAPH_H
#define __PHG_HYPERGRAPH_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


  
typedef struct {
  int info;    /* primarily for debugging recursive algorithms;initially 0 */
  int RowBlocks;  /* number of processor blocks containing row information */
  int ColBlocks;  /* number of processor blocks containing column information */
                  /* RowBlocks * ColBlocks should equal number of processors! */
  int myRowBlock; /* my processor's row block number in [0,RowBlock-1] */
  int myColBlock; /* my processor's column block number in [0,ColBlock-1] */
  int *vtxdist;   /* distributions of vertices to processor blocks. Vertices
                 vtxdist[n] to vtxdist[n+1]-1 are stored on processor block n */
  int *hedgedist; /* distribution of hyperedges to processor blocks.   */                  
   
  int global_nVtx;  /* global number of vertices, |V| */
  int global_nEdge; /* global number of hyperedges, |E| */
  int nVtx;         /* number of vertices on this processor */
  int nEdge;        /* number of hyperedges on this processor */
  int nInput;       /* number of inputs, |I|  on this processor */
  
  int nDim;        /* Number of dimensions of a vertex's coordinate */
  int VertexWeightDim;  /* number of weight dimensions for a vertex */
  int EdgeWeightDim;    /* number of weight dimensions for an edge */
  int redl;             /* Working Reduction limit. */

  /* physical coordinates of each vertex, optional */
  double *coor; /*  |V| long by CoordinateDim */

  /* arrays with vertex and edge weights */
  float *vwgt;  /* weights of vertices, |V| long by VtxWeightDim */
  float *ewgt;  /* weights of hypergraph edges, |E| long by EdgeWeightDim */

  /* arrays to look up vertices given a hyperedge */
  int *hindex;  /* length |E|+1 index into hvertex, last is |P| */
  int *hvertex; /* length |P| array containing associated vertices */

  /* arrays to look up hyperedges given a vertex */
  int *vindex;  /* length |V|+1 index into vedge, last is |P| */
  int *vedge;   /* length |P| array containing associated hyperedges */

 
                     
  int *vmap;     /* used when recursively dividing for p > 2 */
  double ratio;  /* split when recursively dividing for p > 2 */
  } PHGraph;


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HYPERGRAPH_H */
