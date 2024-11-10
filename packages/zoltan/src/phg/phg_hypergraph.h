// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PHG_HYPERGRAPH_H
#define __PHG_HYPERGRAPH_H

#include "phg_comm.h"

#ifdef CEDRIC_2D_PARTITIONS
#include "zoltan_dd.h"
#endif

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

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
  int nRepartVtx;       /* number of repartition vertices added 
                           when LB_APPROACH=repartition */
  int nRepartEdge;      /* number of repartition edges added 
                           when LB_APPROACH=repartition. */
  int nRepartPin;       /* number of repartition pins added 
                           when LB_APPROACH=repartition. */
  
  int VtxWeightDim;     /* number of weight dimensions for a vertex */
  int EdgeWeightDim;    /* number of weight dimensions for a hyperedge */

  /* arrays with vertex and edge weights */
  float *vwgt;    /* weights of vertices, nVtx long by VtxWeightDim */
  float *ewgt;    /* weights of hypergraph edges, nEdge long by EdgeWeightDim */
    
  /* physical coordinates of each vertex, optional */
  int nDim;         /* number of coordinate dimensions for a vertex */
  double *coor;     /* |V| long by CoordinateDim */
  ZOLTAN_GNO_TYPE *esize;      /* global edge size; nEdge long */

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
    
  int *fixed_part;       /* an array of size |V| containing part assignments of fixed
                       vertices. If it is NULL no vertex is fixed in any part; 
                       otherwise 
                         < 0 (negative) : vertex is free,
                         [0,..., p-1]: part that vertex is fixed
                         >= p is invalid. */
  int *pref_part;   /* parallel to fixed array with size of |V| containing the
                       preferred parts for each vertex; again -1 is free.
                       If both a vertex has both prefered and fixed part;
                       fixed part takes precedence. */ 
  int bisec_split;  /* For fixed vertex partitioning via Recursive Bisection
                       treat vertices fixed in parts < bisec_split as
                       they were in the part 0 and the others in part 1.
                    if bisec_split < 0 
                       it is k-way partitioning use fixed parts as they apear
                    */ 
    
    
  PHGComm *comm;  /* this is a pointer to storage PHGPartParamsStruct:
                     (set in phg_build)
                     UVCUVC: I included here because nProc_x, nProc_y was here
                     for convenience.  */
  ZOLTAN_GNO_TYPE *dist_x;    /* distributions of vertices to processor columns. Vertices
                               * dist_x[n] to dist_x[n+1]-1 are stored in col block n */
  ZOLTAN_GNO_TYPE *dist_y;    /* distribution of hyperedges to processor rows as above */    
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
 * Data structure for hypergraph returned by query functions. 
 * Includes Zoltan IDs corresponding to local objects (vertices).
 * Includes HGraph built by Zoltan for use by the partitioning algorithms. 
 ********************************************************************/

struct Zoltan_HGraph {
  int nObj;                 /* Number of on-processor objects. */
  ZOLTAN_GNO_TYPE globalObj;       /* Global number of objects over all procs. */
  int objWeightDim;         /* Number of weights per object (incl ADD_OBJ_WEIGHT) */
  float *objWeight;         /* weight */
  ZOLTAN_GNO_TYPE *objGNO;   /* object global number */
  ZOLTAN_ID_PTR objGID;     /* user's object global ID */
  ZOLTAN_ID_PTR objLID;     /* user's object local ID */
  int *numHEdges;           /* number of hyperedges containing object */

  int *fixed;               /* size nObj, part assignments for fixed vertices */

  int GnRepartVtx;          /* Global number of repartition vtxs added for
                               LB_APPROACH=repartition. */
  ZOLTAN_GNO_TYPE GnRepartEdge;  /* Global number of repartition edges added for
                               LB_APPROACH=repartition. */

  int *Input_Parts;         /* Initial partition #s for on-processor objects */
  int *Output_Parts;        /* Final partition #s for on-processor objects */

  int *AppObjSizes;         /* Object sizes for on-processor objects */
  int showMoveVol;          /* compute and show move (migration) volume */

  double *coor;             /* Array of gathered coordinates returned from */
                            /* Zoltan_Get_Coordinates */
  
  /* This hyperedge list includes all hyperedges when LB_Eval uses ZHG, and
     it includes only removed edges when Zoltan_PHG uses ZHG.
   */
  int nHedges;              /* # of hyperedges */
  ZOLTAN_GNO_TYPE globalHedges;    /* global number of hyperedges listed here */
  ZOLTAN_GNO_TYPE *edgeGNO;        /* edge global number */
  int *Esize;               /* number of vertices in hyperedge */
  int edgeWeightDim;        /* Number of weights stored for each hyperedge */
  float *Ewgt;              /* Edge weights for each hyperedge */
  ZOLTAN_GNO_TYPE *pinGNO;        /* pin global number NEW */
  int *Pin_Procs;           /* Process owning each pin  */
  int nPins;                /* total number of pins in listed edges */
  ZOLTAN_GNO_TYPE globalPins;  /* global number of pins */

  int nRecv_GNOs;              /* Number of GNOs in Recv_GNOs. */
  ZOLTAN_GNO_TYPE *Recv_GNOs;  /* Vertex GNOs of vtxs in 2D decomposition
                               received from other processors in row.
                               Used to fill buffer for Comm_Do_Reverse
                               with VtxPlan in building return lists. */
  ZOLTAN_COMM_OBJ *VtxPlan; /* Communication plan mapping GIDs to GNOs 
                               within row communicators. */
#ifdef CEDRIC_2D_PARTITIONS
  struct Zoltan_DD_Struct *ddHedge;
#endif /* CEDRIC_2D_PARTITIONS */

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
typedef int *Partition; /* length |V|, partition ID for each vertex */


/* Hypergraph utilities */
extern void Zoltan_HG_HGraph_Init (HGraph*);
extern int Zoltan_HG_HGraph_Free  (HGraph*);
extern int Zoltan_HG_Create_Mirror(ZZ*, HGraph*);
extern void Zoltan_HG_Mirror(int, int*, int*, int, int*, int*);

extern int Zoltan_HG_Info         (ZZ*, HGraph*);
extern int Zoltan_HG_Check        (ZZ*, HGraph*);
extern void Zoltan_HG_Print(ZZ*, HGraph*, Partition, FILE*, char*);
extern void Zoltan_HG_HGraph_Print(ZZ *zz, ZHG *, HGraph *, Partition, FILE *fp);
    
extern void Zoltan_Input_HG_Init (ZHG *);
extern int Zoltan_Input_HG_Free  (ZHG *);
    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __PHG_HYPERGRAPH_H */
