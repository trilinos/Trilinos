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

#ifndef ZOLTAN_HYPERGRAPH_H
#define ZOLTAN_HYPERGRAPH_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <string.h>

#ifdef HTEST
typedef struct
  {
  int Num_GID;
  int Num_LID;
  int Debug_Level;
  int Obj_Weight_Dim;
  int Edge_Weight_Dim;
  int Proc;
  } ZZ ;
#define ZOLTAN_DEBUG_LIST            1
#define ZOLTAN_DEBUG_ALL             1
#define ZOLTAN_TRACE_ENTER(a,b)      {}
#define ZOLTAN_TRACE_EXIT(a,b)       {}
#define ZOLTAN_TRACE_DETAIL(a,b,c)   {}
#define ZOLTAN_MALLOC                malloc
#define ZOLTAN_FREE(a)               {free(*(a)); *a = NULL;}
#define ZOLTAN_OK                    0
#define ZOLTAN_WARN                  1
#define ZOLTAN_FATAL                 -1
#define ZOLTAN_MEMERR                -2
#define ZOLTAN_PRINT_ERROR(a,b,c)    printf("Zoltan Error %s in %s\n",c,b)
#define ZOLTAN_PRINT_WARN(a,b,c)     printf("Zoltan Warning %s\n",c)
#define MAX_PARAM_STRING_LEN 50
#else
#include "zz_const.h"
#include "params_const.h"
#endif

#define MIN(A,B)                (((A)<(B))?(A):(B))
#define MAX(A,B)                (((A)>(B))?(A):(B))

#define RANDOM_SEED             123456789   /* (time ((time_t*)NULL)) */
#define EPS                     1e-6


typedef int* Matching ;  /* length |V|, matching information of vertices */
typedef int* Packing ;   /* length |V|, packing information of vertices */
typedef int* Grouping ;  /* length |V|, grouping information of vertices */
typedef int* LevelMap ;  /* length |V|, mapping of fine vertices onto coarse vertices */
typedef int* Partition ; /* length |V|, partition ID for each vertex */


typedef struct
{
   int info ;    /* primarily for debugging recursive algorithms;initially 0 */
   int nVtx ;    /* number of vertices, |V| */
   int nEdge ;   /* Size of neigh array; 2|E| */
   int nDim ;    /* Number of dimensions for a vertex's coordinate */
   int VertexWeightDim ;  /* number of weight dimensions for a vertex */
   int EdgeWeightDim ;    /* number of weight dimensions for an edge */

   int *vtxdist;  /* distributions of vertices to processors, as in ParMETIS.
                     Vertices vtxdist[n] to vtxdist[n+1]-1 are stored on
                     processor n.   KDD:  temporary; may change later. */

   /* physical coordinates of each vertex, optional */
   double *coor ; /*  |V| long by CoordinateDim */
  
   /* arrays with vertex and edge weights */
   float *vwgt ;  /* weights of vertices, |V| long by VtxWeightDim */
   float *ewgt ;  /* weights of hypergraph edges, 2|E| long by EdgeWeightDim */

   /* arrays to look up the neighbors of a vertex */
   int *nindex ;  /* length |V|+1 index to neigh, last is 2|E| */
   int *neigh ;   /* length 2|E|, list of neighbors for each vertex */
} Graph ;


typedef struct
{
   int info;     /* primarily for debugging recursive algorithms;initially 0 */
   int nVtx ;    /* number of vertices, |V| */
   int nEdge ;   /* number of hyperedges, |E| */
   int nPin ;    /* number of pins, |P| */
   int nDim ;    /* Number of dimensions of a vertex's coordinate */
   int VertexWeightDim ;  /* number of weight dimensions for a vertex */
   int EdgeWeightDim ;    /* number of weight dimensions for an edge */

   /* physical coordinates of each vertex, optional */
   double *coor ; /*  |V| long by CoordinateDim */

   /* arrays with vertex and edge weights */
   float *vwgt ;  /* weights of vertices, |V| long by VtxWeightDim */
   float *ewgt ;  /* weights of hypergraph edges, |E| long by EdgeWeightDim */

   /* arrays to look up vertices given a hyperedge */
   int *hindex ;  /* length |E|+1 index into hvertex, last is |P| */
   int *hvertex ; /* length |P| array containing associated vertices */

   /* arrays to look up hyperedges given a vertex */
   int *vindex ;  /* length |V|+1 index into vedge, last is |P| */
   int *vedge ;   /* length |P| array containing associated hyperedges */
} HGraph ;


typedef struct
{
   HGraph *hg ;
   Partition part ;
   Packing pack ;
   LevelMap lmap ;
} HGraphLevel ;


/* Hypergraph utilities */
void Zoltan_HG_HGraph_Init (HGraph *);
void Zoltan_HG_Graph_Init  (Graph *);
int Zoltan_HG_HGraph_Free  (HGraph *);
int Zoltan_HG_Graph_Free   (Graph *);
int Zoltan_HG_Info         (ZZ *, HGraph *);
int Zoltan_HG_Create_Mirror(ZZ *, HGraph *);
int Zoltan_HG_Check        (ZZ *, HGraph *);
int Zoltan_HG_HGraph_to_Graph(ZZ *, HGraph *, Graph *);
int Zoltan_HG_Graph_to_HGraph(ZZ *, Graph *, HGraph *);

/* Hypergraph read from file */
int Zoltan_HG_Readfile     (ZZ *, HGraph *, char *hgraphfile);

/* Hypergraph Partitioning */
/* Function types for options to hypergraph partitioning */
typedef int ZOLTAN_HG_MATCHING_FN (ZZ *, Graph *,  Matching, int);
typedef int ZOLTAN_HG_PACKING_FN  (ZZ *, HGraph *, Packing,  int);
typedef int ZOLTAN_HG_GROUPING_FN (ZZ *, HGraph *, Grouping, int);
typedef int ZOLTAN_HG_GLOBAL_PART_FN(ZZ *, HGraph *, int, Partition);
typedef int ZOLTAN_HG_LOCAL_REF_FN(ZZ *, HGraph *);
/* Parameters to the hypergraph functions */
typedef struct {
  int redl;                              /* Reduction limit. */
  char redm_str[MAX_PARAM_STRING_LEN];   /* Reduction method string. */
  ZOLTAN_HG_MATCHING_FN *matching;       /* Pointers to Matching, Packing and */
  ZOLTAN_HG_PACKING_FN  *packing;        /*  Grouping fn specified by */
  ZOLTAN_HG_GROUPING_FN *grouping;       /*  redm_str; NULL if not used */
  char global_str[MAX_PARAM_STRING_LEN]; /* Global partitioning string and */
  ZOLTAN_HG_GLOBAL_PART_FN *global_part; /*  pointer to Global partitioning fn */
  char local_str[MAX_PARAM_STRING_LEN];  /* Local refinement string and */
  ZOLTAN_HG_LOCAL_REF_FN *local_ref;     /*  pointer to Local refinement fn */
  int check_graph;                       /* Flag to indicate input hgraph
                                            checking level*/
} HGPartParams;
int Zoltan_HG_Set_Options  (ZZ *, HGPartParams *);
int Zoltan_HG_HPart_Lib    (ZZ *, HGraph *, int, Partition, HGPartParams *);
int Zoltan_HG_HPart_Info   (ZZ *, HGraph *, int, Partition);
int Zoltan_HG_Scale_Graph_Weight (ZZ *, Graph *);
int Zoltan_HG_Scale_HGraph_Weight (ZZ *, HGraph *, float *);

/* Matching functions */
int Zoltan_HG_Matching (ZZ *, Graph *, Matching, HGPartParams *, int);
ZOLTAN_HG_MATCHING_FN *Zoltan_HG_Set_Matching_Fn(char *);

/* Packing */
int Zoltan_HG_Packing (ZZ *, HGraph *, Packing, HGPartParams *);
ZOLTAN_HG_PACKING_FN *Zoltan_HG_Set_Packing_Fn(char *);

/* Grouping */
int Zoltan_HG_Grouping (ZZ *, HGraph *, Packing, HGPartParams *);
ZOLTAN_HG_GROUPING_FN *Zoltan_HG_Set_Grouping_Fn(char *);

/* Coarsening */
int Zoltan_HG_Coarsening   (ZZ *, HGraph *, Packing, HGraph *, int *);

/* Global Partitioning functions */
int Zoltan_HG_Global (ZZ *, HGraph *, int, Partition, HGPartParams *);
ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char *);

/* Local refinement functions */ /* KDD Placeholder for later. */
ZOLTAN_HG_LOCAL_REF_FN *Zoltan_HG_Set_Local_Ref_Fn(char *);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HYPERGRAPH_H */
