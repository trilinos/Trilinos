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
#include "zoltan_mem.h"
typedef struct
  {
  int Debug_Level;
  } ZZ ;
#define ZOLTAN_DEBUG_LIST            1
#define ZOLTAN_DEBUG_ALL             2
#define ZOLTAN_TRACE_ENTER(a,b)      {}
#define ZOLTAN_TRACE_EXIT(a,b)       {}
#define ZOLTAN_TRACE_DETAIL(a,b,c)   {}
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

unsigned long Zoltan_HG_Rand (void) ;
void          Zoltan_HG_Srand (unsigned long) ;

/* Hypergraph read from file */
int HG_Readfile (ZZ *, HGraph *, char *hgraphfile);

/* Hypergraph Partitioning */
/* Function types for options to hypergraph partitioning */
typedef int ZOLTAN_HG_MATCHING_FN   (ZZ *, HGraph *, Graph *,  Matching, int*);
typedef int ZOLTAN_HG_PACKING_FN    (ZZ *, HGraph *, Packing,  int*);
typedef int ZOLTAN_HG_GROUPING_FN   (ZZ *, HGraph *, Grouping, int*);
typedef int ZOLTAN_HG_GLOBAL_PART_FN(ZZ *, HGraph *, int, Partition);
typedef int ZOLTAN_HG_LOCAL_REF_FN  (ZZ *, HGraph *, int, Partition, float);

/* Function types for edge-weight scaling functions. */
/* Placeholder for now; if these functions end up having the same argument */
/* list for each type, do not need separate types here or separate pointers */
/* in HGPartParams.  KDD */
typedef int ZOLTAN_HG_MATCHING_EWS_FN (ZZ *, Graph *);
typedef int ZOLTAN_HG_PACKING_EWS_FN  (ZZ *, HGraph *);
typedef int ZOLTAN_HG_GROUPING_EWS_FN (ZZ *, HGraph *);


/* Parameters to the hypergraph functions */
typedef struct {
  float bal_tol;                        /*Balance tolerance in % of average */
  int redl;                             /*Reduction limit. */

  char redm_str[MAX_PARAM_STRING_LEN];  /*Reduction method string. */
  ZOLTAN_HG_MATCHING_FN *matching;      /*Pointers to Matching, Packing and */
  ZOLTAN_HG_PACKING_FN  *packing;       /* Grouping fn specified by */
  ZOLTAN_HG_GROUPING_FN *grouping;      /* redm_str; NULL if not used */

  char redmo_str[MAX_PARAM_STRING_LEN]; /* Matching optimization string*/
  ZOLTAN_HG_MATCHING_FN *matching_opt;  /* Pointers to Matching, Packing and */
  ZOLTAN_HG_PACKING_FN  *packing_opt;   /* Grouping optimization fn specified*/
  ZOLTAN_HG_GROUPING_FN *grouping_opt;  /* by redmo_str; NULL if not used */

  int ews;                              /* Flag for Edge weight scaling */

  char global_str[MAX_PARAM_STRING_LEN];/*Global partitioning string and */
  ZOLTAN_HG_GLOBAL_PART_FN *global_part;/* pointer to Global partitioning fn */

  char local_str[MAX_PARAM_STRING_LEN]; /*Local refinement string and */
  ZOLTAN_HG_LOCAL_REF_FN *local_ref;    /* pointer to Local refinement fn */

  int check_graph;                      /* Flag indicating whether the input
                                           hypergraph should be checked for 
                                           errors. */
} HGPartParams;

int Zoltan_HG_Set_Part_Options  (ZZ *, HGPartParams *);
int Zoltan_HG_HPart_Lib    (ZZ *, HGraph *, int, Partition, HGPartParams *);
int Zoltan_HG_HPart_Info   (ZZ *, HGraph *, int, Partition);
float hcut_size_total (HGraph *, Partition);
float hcut_size_links (ZZ *, HGraph *, int, Partition);

/* Scale Edge Weight */
int Zoltan_HG_Scale_Graph_Weight  (ZZ *, Graph *, float *, int);
int Zoltan_HG_Scale_HGraph_Weight (ZZ *, HGraph *, float *);

/* Matching functions */
int Zoltan_HG_Matching (ZZ *, HGraph *, Graph *, Matching, HGPartParams *, int*);
int Zoltan_HG_Set_Matching_Fn(HGPartParams *);

/* Packing */
int Zoltan_HG_Packing (ZZ *, HGraph *, Packing, HGPartParams *, int*);
int Zoltan_HG_Set_Packing_Fn(HGPartParams *);

/* Grouping */
int Zoltan_HG_Grouping (ZZ *, HGraph *, Packing, HGPartParams *, int*);
int Zoltan_HG_Set_Grouping_Fn(HGPartParams *);

/* Coarsening */
int Zoltan_HG_Coarsening   (ZZ *, HGraph *, Packing, HGraph *, int *);

/* Global Partitioning functions */
int Zoltan_HG_Global (ZZ *, HGraph *, int, Partition, HGPartParams *);
ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char *);

/* Local refinement functions */ /* KDD Placeholder for later. */
int Zoltan_HG_Local (ZZ *, HGraph *, int, Partition, HGPartParams *);
ZOLTAN_HG_LOCAL_REF_FN *Zoltan_HG_Set_Local_Ref_Fn(char *);

/* Sorting */
void quicksort_pointer_dec_float_int (int*, float*, int*, int, int);
void quicksort_pointer_dec_float     (int*, float*, int, int);
void quicksort_pointer_inc_int_int   (int*, int*, int*, int, int);
void quicksort_list_inc_int          (int*, int, int);
void quicksort_pointer_inc_int_mult  (int *, int, int, int*, int*);

extern int Zoltan_HG_Readfile ( int, FILE *, int *, int *, int *,
 int **, int **, int *, float **, int *, float **, int *);
extern void Zoltan_HG_Print(ZZ *, HGraph *);

/* Heap datastructure */
typedef struct
{ int space;
  int n;
  int *ele;
  int *pos;
  float *d;
} HEAP;
#define heap_empty(H)         (((H)->n)==0)
#define heap_not_empty(H)     (((H)->n)!=0)
extern int  heap_init         (HEAP*, int);
extern void heap_free         (HEAP*);
extern int  heap_check        (HEAP*);
extern void heap_input        (HEAP*, int, float);
extern void heap_make         (HEAP*);
extern void heapify           (HEAP*, int);
extern void heap_change_key   (HEAP*, int, float);
extern int  heap_extract_max  (HEAP*);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HYPERGRAPH_H */
