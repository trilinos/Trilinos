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

#ifndef __ZOLTAN_PHG_H
#define __ZOLTAN_PHG_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "phg_const.h"
#include "phg_hypergraph.h"
#include "phg_util.h"
#include "params_const.h"
#include "zoltan_comm.h"

/********************************************************************
 * Data structure for Zoltan's base hypergraph. Includes Zoltan IDs 
 * corresponding to local objects (vertices) and a PHGraph as used 
 * by the partitioning algorithms. 
 ********************************************************************/

struct Zoltan_PHGraph {
  int nObj;                 /* Number of on-processor objects. */
  ZOLTAN_ID_PTR Global_IDs; /* Global IDs for on-processor objects.  */
  ZOLTAN_ID_PTR Local_IDs;  /* Local IDs for on-processor objects.   */
  int *Parts;               /* Initial partition #s for on-processor objects */
                            /* KDD In parallel version Part may be in HG */
  PHGraph PHG;              /* Hypergraph for initial objects.       */
};
typedef struct Zoltan_PHGraph ZPHG;

/************************************************/
/* Mappings supporting the 2D data distribution */
/************************************************/

/* Mapping of global number to local number           */
/* Code should call VTX_GNO_TO_LNO or EDGE_GNO_TO_LNO */

#define GNO_TO_LNO(gno, dist, myblock) \
    (gno) - (dist)[(myblock)]
#define VTX_GNO_TO_LNO(phg, gno) \
    GNO_TO_LNO(gno, (phg)->dist_x, (phg)->myProc_x)
#define EDGE_GNO_TO_LNO(phg, gno) \
    GNO_TO_LNO(gno, (phg)->dist_y, (phg)->myProc_y)


/* Mapping of local number to global number           */
/* Code should call VTX_LNO_TO_GNO or EDGE_LNO_TO_GNO */

#define LNO_TO_GNO(lno, dist, myblock) \
    (lno) + (dist)[(myblock)]
#define VTX_LNO_TO_GNO(phg, lno) \
    LNO_TO_GNO(lno, (phg)->dist_x, (phg)->myProc_x)
#define EDGE_LNO_TO_GNO(phg, lno) \
    LNO_TO_GNO(lno, (phg)->dist_y, (phg)->myProc_y)


/* Mapping of global number to processor block.     */
/* Code should call EDGE_TO_PROC_Y or VTX_TO_PROC_X */

#define EDGE_TO_PROC_Y(phg, gno) \
    Zoltan_PHG_Gno_To_Proc_Block((gno), (phg)->dist_y, (phg)->nProc_y)

#define VTX_TO_PROC_X(phg, gno) \
    Zoltan_PHG_Gno_To_Proc_Block(gno, (phg)->dist_x, (phg)->nProc_x)


/*****************/
/* Useful macros */
/*****************/

#define MIN(A,B)  (((A) < (B)) ? (A) : (B))
#define MAX(A,B)  (((A) > (B)) ? (A) : (B))

#define RANDOM_SEED   123456789   /* (time ((time_t*)NULL)) */
#define EPS           1e-6        /* small value, like a machine epsilon */

/********************/
/* Type definitions */
/********************/

typedef int *Matching;  /* length |V|, matching information of vertices */
/* If a vertex i is not being contracted (matched) with other vertices,  
 * Matching[i] == i.  If vertices i, j, and k are being contracted together to
 * form one new vertex, Matching[i] == j; Matching[j] == k; and 
 * Matching[k] == i. The cycle (of any length) describes the contraction. 
 */
 
typedef int *LevelMap;  /* length |V|, mapping of fine vertices onto coarse 
                           vertices */
typedef int *Partition; /* length |V|, partition ID for each vertex */

typedef struct {
   PHGraph   *hg;
   Partition part;
   Matching  match;
   LevelMap  lmap;
} PHGraphLevel;

/* Function types for options to hypergraph partitioning */

struct PHGPartParamsStruct;  /* Forward declaration */

typedef int ZOLTAN_PHG_MATCHING_FN(ZZ*, PHGraph*, Matching);
typedef int ZOLTAN_PHG_COARSEPARTITION_FN(ZZ*, PHGraph*, int, Partition,
                                          struct PHGPartParamsStruct*);
typedef int ZOLTAN_PHG_REFINEMENT_FN(ZZ*, PHGraph*, int, Partition,
                                     struct PHGPartParamsStruct*, float);

/* Function types for edge-weight scaling functions. Placeholder for now; */
/* if these functions end up having the same argument list for each type, */
/* do not need separate types here or separate pointers in HGPartParams.  KDD */
typedef int ZOLTAN_PHG_MATCHING_EWS_FN(ZZ*, PGraph*);

/******************************************/
/* Parameters to the hypergraph functions */
/******************************************/
struct PHGPartParamsStruct {
  float bal_tol;                 /* Balance tolerance in % of average */
  int kway;                      /* 1 -> direct kway, 0->recursive bisection */
  int redl;                      /* Reduction limit (constant). */
  char redm_str[MAX_PARAM_STRING_LEN];   /* Reduction method string. */
  ZOLTAN_PHG_MATCHING_FN *matching;      /* Pointers to Matching function */
  ZOLTAN_PHG_MATCHING_FN *matching_opt;  /* Pointers to Matching opt function */
  int ews;                               /* type of hyperedge weight scaling */
  char coarsepartition_str[MAX_PARAM_STRING_LEN]; 
                                         /* Coarse partitioning string */
  ZOLTAN_PHG_COARSEPARTITION_FN *CoarsePartition;
                                         /* pointer to coarse partitioning fn */
  char refinement_str[MAX_PARAM_STRING_LEN]; /* Refinement string and */
  ZOLTAN_PHG_REFINEMENT_FN *Refinement;      /* pointer to refinement fn */
  int check_graph;      /* Flag indicating whether the input hypergraph should 
                         * be checked for errors. */
  int output_level;     /* Flag indicating amount of output from HG algorithms.
                         * See levels PHG_DEBUG_* below.  */
  int nProc_x;    /* number of processors in x-direction of 2D data distrib.  */
  int nProc_y;    /* number of processors in y-direction of 2D data distrib.  */
                  /* nProc_x * nProc_y should equal number of processors!     */
  int myProc_x;   /* my processor's row block number in [0,nProc_x-1] */
  int myProc_y;   /* my processor's column block number in [0,nProc_y-1] */
};

typedef struct PHGPartParamsStruct PHGPartParams;

/*****************************/
/* Hypergraph output levels: */
/*****************************/
#define PHG_DEBUG_NONE 0
#define PHG_DEBUG_LIST 1
#define PHG_DEBUG_ALL  2
#define PHG_DEBUG_PRINT 3
#define PHG_DEBUG_PLOT 4


/*********************/
/* Scale Edge Weight */
/*********************/
int Zoltan_PHG_Scale_HGraph_Weight (ZZ*, PHGraph*, float*, int);
int Zoltan_PHG_Scale_Graph_Weight  (ZZ*, PGraph*,  float*, int);

/**********************/
/* Matching functions */
/**********************/
int Zoltan_PHG_Matching (ZZ*, PHGraph*, Matching, PHGPartParams*);
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams*);

/**************/
/* Coarsening */
/**************/
int Zoltan_PHG_Coarsening   (ZZ*, PHGraph*, Matching, PHGraph*, int*);

/*********************************/
/* Coarse Partitioning functions */
/*********************************/
int Zoltan_PHG_CoarsePartition (ZZ*, PHGraph*, int, Partition, PHGPartParams*);
ZOLTAN_PHG_COARSEPARTITION_FN *Zoltan_PHG_Set_CoarsePartition_Fn(char*);

/************************/
/* Refinement functions */ 
/************************/
int Zoltan_PHG_Refinement (ZZ*, PHGraph*, int, Partition, PHGPartParams*);
ZOLTAN_PHG_REFINEMENT_FN *Zoltan_PHG_Set_Refinement_Fn(char*);

/***********/
/* Sorting */
/***********/
void Zoltan_quicksort_pointer_dec_float_int (int*, float*, int*, int, int);
void Zoltan_quicksort_pointer_dec_float     (int*, float*, int,  int);
void Zoltan_quicksort_pointer_inc_int_int   (int*, int*,   int*, int, int);
void Zoltan_quicksort_list_inc_int          (int*, int,    int);
void Zoltan_quicksort_pointer_inc_int_mult  (int*, int,    int,  int*, int*);

/**************************************/
/* Heap data structure and operations */
/**************************************/
typedef struct {
   int    space;
   int    n;
   int   *ele;
   int   *pos;
   float *value;
} HEAP;

#define Zoltan_PHG_heap_empty(H)         (((H)->n)==0)
#define Zoltan_PHG_heap_not_empty(H)     (((H)->n)!=0)
#define Zoltan_PHG_heap_max_value(H)     ((H)->value[(H)->ele[0]])
#define Zoltan_PHG_heap_peek_max(H)      ((H)->ele[0])
#define Zoltan_PHG_heap_count(H)         ((H)->n)

int  Zoltan_PHG_heap_init         (ZZ*, HEAP*, int);
void Zoltan_PHG_heap_clear        (HEAP*);
void Zoltan_PHG_heap_free         (HEAP*);
int  Zoltan_PHG_heap_check        (HEAP*);
int  Zoltan_PHG_heap_input        (HEAP*, int, float);
int  Zoltan_PHG_heap_make         (HEAP*);
int  Zoltan_PHG_heap_change_value (HEAP*, int, float);
int  Zoltan_PHG_heap_extract_max  (HEAP*);
int  Zoltan_PHG_heap_extract      (HEAP*, int);

/*****************/
/* Communication */
/*****************/
extern int Zoltan_PHG_gather_slice_root(
        int, int, int, int, int, char*, int*, char**, MPI_Comm*, int);
extern int Zoltan_PHG_gather_slice(
        int, int, int, int, int, char*, int*, char**, MPI_Comm*, int);

/*****************************/
/* Other Function Prototypes */
/*****************************/

extern int Zoltan_PHG_Set_Part_Options(ZZ*, PHGPartParams*);
extern int Zoltan_PHG_HPart_Lib(ZZ*, PHGraph*, int, Partition, PHGPartParams*, 
                                int);
extern int Zoltan_PHG_HPart_Info(ZZ*, PHGraph*, int, Partition, PHGPartParams*);
extern double Zoltan_PHG_hcut_size_total(PHGraph*, Partition);
extern double Zoltan_PHG_hcut_size_links(ZZ*, PHGraph*, Partition);
extern double Zoltan_PHG_HPart_balance(ZZ*, PHGraph*, int, Partition);

extern int Zoltan_PHG_move_vertex(PHGraph*, int, int, int, int*, int**, double*,
                                  HEAP*);

extern int Zoltan_PHG_Build_Hypergraph(ZZ*, ZPHG**, PHGPartParams*);
extern void Zoltan_PHG_HGraph_Print(ZZ*, ZPHG*,  PHGraph*, FILE*);
extern void Zoltan_PHG_Plot(int, int, int, int*, int*, int*, char*);
extern void Zoltan_PHG_Plot_2D_Distrib(ZZ*, PHGraph*);

extern int Zoltan_PHG_Gno_To_Proc_Block(int gno, int*, int);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_H */
