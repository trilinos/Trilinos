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

#ifndef __HG_HYPERGRAPH_H
#define __HG_HYPERGRAPH_H

#include "phg_comm.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

    
typedef struct {
   int info;    /* primarily for debugging recursive algorithms;initially 0 */
   int nVtx;    /* number of vertices, |V| */
   int nEdge;   /* Size of neigh array; 2|E| */
   int nDim;    /* Number of dimensions for a vertex's coordinate */
   int VtxWeightDim;  /* number of weight dimensions for a vertex */
   int EdgeWeightDim;    /* number of weight dimensions for an edge */
   int redl;             /* Working Reduction limit. */

   int *vtxdist;  /* distributions of vertices to processors, as in ParMETIS.
                     Vertices vtxdist[n] to vtxdist[n+1]-1 are stored on
                     processor n.   KDD:  temporary; may change later. */

   /* physical coordinates of each vertex, optional */
   double *coor; /*  |V| long by CoordinateDim */

   /* arrays with vertex and edge weights */
   float *vwgt;  /* weights of vertices, |V| long by VtxWeightDim */
   float *ewgt;  /* weights of hypergraph edges, 2|E| long by EdgeWeightDim */

   /* arrays to look up the neighbors of a vertex */
   int *nindex;  /* length |V|+1 index to neigh, last is 2|E| */
   int *neigh;   /* length 2|E|, list of neighbors for each vertex */
   } Graph;


typedef struct {
  int info;             /* primarily for debugging recursive algorithms;initially 0 */
  int nVtx;             /* number of vertices on this processor */
  int nEdge;            /* number of hyperedges on this processor */
  int nPins;            /* number of pins (nonzeros) on this processor */
  
  int VtxWeightDim;     /* number of weight dimensions for a vertex */
  int EdgeWeightDim;    /* number of weight dimensions for a hyperedge */

    /* arrays with vertex and edge weights */
  float *vwgt;    /* weights of vertices, nVtx long by VtxWeightDim */
  float *ewgt;    /* weights of hypergraph edges, nEdge long by EdgeWeightDim */
    
  /* physical coordinates of each vertex, optional */
  int nDim;         /* number of coordinate dimensions for a vertex */
  double *coor;     /* |V| long by CoordinateDim */

  /* arrays to look up vertices given a hyperedge */
  int *hindex;      /* length nEdge+1 index into hvertex, last is nPins */
  int *hvertex;     /* length nPins array containing associated vertices */

  /* arrays to look up hyperedges given a vertex */
  int *vindex;      /* length nVtx+1 index into vedge, last is nPins */
  int *vedge;       /* length nPins array containing associated hyperedges */

    /* UVCUVC: todo vmap, ratio and redl should be removed from HGraph
       and some could go to VCycle struct in RB code */
  int *vmap;        /* used when recursively dividing for p > 2 */
  double ratio;     /* split when recursively dividing for p > 2 */
  int redl;         /* working reduction limit */

    
  PHGComm *comm;  /* this is a pointer to storage PHGPartParamsStruct: (set in phg_build)
                     UVCUVC: I've included here because nProc_x, nProc_y was here
                     for convenience.
                   */
  int *dist_x;    /* distributions of vertices to processor columns. Vertices
                   * dist_x[n] to dist_x[n+1]-1 are stored in col block n */
  int *dist_y;    /* distribution of hyperedges to processor rows as above */                  
    
} HGraph;


/* Matching, Packing, and Grouping arrays.  Technically, they are all the same;
 * the different typedefs are not needed.  In the description below, Matching is
 * used; the same description applies to Packing and Grouping.  If a vertex i is
 * not being contracted with other vertices,  Matching[i] == i.  If vertices i,
 * j, and k are being contracted together to form one new vertex,
 * Matching[i] == j; Matching[j] == k; and Matching[k] == i;
 * The cycle describes the contraction. */
typedef int *Matching;  /* length |V|, matching information of vertices */
typedef int *Packing;   /* length |V|, packing information of vertices */
typedef int *Grouping;  /* length |V|, grouping information of vertices */

typedef int *LevelMap;  /* length |V|, mapping of fine vertices onto coarse vertices */
typedef int *Partition; /* length |V|, partition ID for each vertex */


extern void Zoltan_HG_Graph_Init  (Graph*);
extern int Zoltan_HG_Graph_Free   (Graph*);

/* Hypergraph utilities */
extern void Zoltan_HG_HGraph_Init (HGraph*);
extern int Zoltan_HG_HGraph_Free  (HGraph*);
extern int Zoltan_HG_Create_Mirror(ZZ*, HGraph*);
extern void Zoltan_HG_Mirror(int, int*, int*, int, int*, int*);

extern int Zoltan_HG_Info         (ZZ*, HGraph*);
extern int Zoltan_HG_Check        (ZZ*, HGraph*);
extern int Zoltan_HG_HGraph_to_Graph(ZZ*, HGraph*, Graph*);
extern int Zoltan_HG_Graph_to_HGraph(ZZ*, Graph*,  HGraph*);
extern void Zoltan_HG_Print(ZZ*, HGraph*, Partition, FILE*, char*);


    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HYPERGRAPH_H */
