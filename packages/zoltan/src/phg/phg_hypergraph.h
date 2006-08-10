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

#include "phg_comm.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/********************************************************************
 * A seldom-used graph data structure.
 ********************************************************************/
    
typedef struct {
   int info;    /* depth of V-cycle for this hypergraph; initially 0 */
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


/********************************************************************
 * Hypergraph data structure.  Supports both serial hypergraph and 
 * 2D data distribution (where each block is, in fact, a serial 
 * hypergraph).
 ********************************************************************/
typedef struct {
  int info;             /* depth of recursion in V-cycle; initially 0 */
  int nVtx;             /* number of vertices on this processor */
  int nEdge;            /* number of hyperedges on this processor */
  int nPins;            /* number of pins (nonzeros) on this processor */
  int nRepartVtx;       /* number of repartition vertices added in PHG_REPART */
  int nRepartEdge;      /* number of repartition edges added in PHG_REPART. */
  int nRepartPin;       /* number of repartition pins added in PHG_REPART. */
  
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

  int *fixed;       /* an array of size |V| containing part assignments of fixed
                       vertices. If it is NULL no vertex is fixed in any part; 
                       otherwise 
                         < 0 (negative) : vertex is free,
                         [0,..., p-1]: part that vertex is fixed
                         >= p is invalid. */
  int bisec_split;  /* For fixed vertex partitioning via Recursive Bisection
                       treat vertices fixed in parts < bisec_split as
                       they were in the part 0 and the others in part 1.
                    if bisec_split < 0 
                       if fixed!=NULL  it is k-way partitioning use fixed
                                      parts as they apear
                       if fixed==NULL not fixed vertex partitioning. */ 
    
    
  PHGComm *comm;  /* this is a pointer to storage PHGPartParamsStruct:
                     (set in phg_build)
                     UVCUVC: I included here because nProc_x, nProc_y was here
                     for convenience.  */
  int *dist_x;    /* distributions of vertices to processor columns. Vertices
                   * dist_x[n] to dist_x[n+1]-1 are stored in col block n */
  int *dist_y;    /* distribution of hyperedges to processor rows as above */    
} HGraph;



/************************************************/
/* Mappings supporting the 2D data distribution */
/************************************************/

/* Mapping of global number to local number           */
/* Code should call VTX_GNO_TO_LNO or EDGE_GNO_TO_LNO */

#define GNO_TO_LNO(gno, dist, myblock) \
    ((gno) - (dist)[(myblock)])
#define VTX_GNO_TO_LNO(phg, gno) \
    GNO_TO_LNO(gno, (phg)->dist_x, (phg)->comm->myProc_x)
#define EDGE_GNO_TO_LNO(phg, gno) \
    GNO_TO_LNO(gno, (phg)->dist_y, (phg)->comm->myProc_y)


/* Mapping of local number to global number           */
/* Code should call VTX_LNO_TO_GNO or EDGE_LNO_TO_GNO */

#define LNO_TO_GNO(lno, dist, myblock) \
    ((lno) + (dist)[(myblock)])
#define VTX_LNO_TO_GNO(phg, lno) \
    LNO_TO_GNO(lno, (phg)->dist_x, (phg)->comm->myProc_x)
#define EDGE_LNO_TO_GNO(phg, lno) \
    LNO_TO_GNO(lno, (phg)->dist_y, (phg)->comm->myProc_y)


/* Mapping of global number to processor block.     */
/* Code should call EDGE_TO_PROC_Y or VTX_TO_PROC_X */

#define EDGE_TO_PROC_Y(phg, gno) \
    Zoltan_PHG_Gno_To_Proc_Block((gno), (phg)->dist_y, (phg)->comm->nProc_y)

#define VTX_TO_PROC_X(phg, gno) \
    Zoltan_PHG_Gno_To_Proc_Block(gno, (phg)->dist_x, (phg)->comm->nProc_x)
    

/********************************************************************
 * Data structure for Zoltan's base hypergraph. Includes Zoltan IDs 
 * corresponding to local objects (vertices) and a HGraph as used 
 * by the partitioning algorithms. 
 ********************************************************************/

struct Zoltan_HGraph {
  int nObj;                 /* Number of on-processor objects. */
  int GnObj;                /* Global number of objects over all procs. */
  ZOLTAN_ID_PTR GIDs;       /* Global IDs for on-processor objects.  */
  ZOLTAN_ID_PTR LIDs;       /* Local IDs for on-processor objects.   */
  int GnInputParts;         /* Global number of input parts == 
                               max(Input_Parts[i]) + 1 for i=0,...,nProcs-1 */
  int *Input_Parts;         /* Initial partition #s for on-processor objects */
  int *Output_Parts;        /* Final partition #s for on-processor objects */
  int nRemove;              /* # of input hyperedges removed */
  ZOLTAN_ID_PTR Remove_EGIDs;/* GIDs of removed hyperedges */
  ZOLTAN_ID_PTR Remove_ELIDs;/* LIDs of removed hyperedges */
  int *Remove_Esize;    /* local size on this proc of each removed hyperedge */
  int *Remove_GEsize;   /* global number of vtx in each removed hyperedge */
  float *Remove_Ewgt;       /* Edge weights for each removed hyperedge */
  ZOLTAN_ID_PTR Remove_Pin_GIDs; /* GIDs of vertices */
  int *Remove_Pin_Procs;         /* Process owning each pin (Pin callbacks only) */
  int nRecv_GNOs;           /* Number of GNOs in Recv_GNOs. */
  int *Recv_GNOs;           /* Vertex GNOs of vtxs in 2D decomposition
                               received from other processors in row.
                               Used to fill buffer for Comm_Do_Reverse
                               with VtxPlan in building return lists. */
  ZOLTAN_COMM_OBJ *VtxPlan; /* Communication plan mapping GIDs to GNOs 
                               within row communicators. */
  HGraph HG;                /* Hypergraph for initial objects.       */
};
typedef struct Zoltan_HGraph ZHG;

/******************************************************************************/
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

typedef int *LevelMap;  /* length |V|, mapping of fine vtxs onto coarse vtxs */
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
extern void Zoltan_HG_HGraph_Print(ZZ *zz, ZHG *, HGraph *, Partition, FILE *fp);
    
struct PHGPartParamsStruct;  /* Forward declaration */

extern int Zoltan_HG_Graph_Callbacks(ZZ *, ZHG *, struct PHGPartParamsStruct *,
  int, float, int, int *,
  ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, float **, int *, ZOLTAN_ID_PTR *,
  int **);
extern int Zoltan_Call_Hypergraph_Pin_Query(ZZ *zz, int *num_lists,
   int *num_pins, ZOLTAN_ID_PTR *edg_GID, int **row_ptr, 
   ZOLTAN_ID_PTR *vtx_GID);
int Zoltan_HG_ignore_some_edges(ZZ *, ZHG *, int, float, int, int *,
  ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *, float *, ZOLTAN_ID_PTR, int *);

    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __PHG_HYPERGRAPH_H */
