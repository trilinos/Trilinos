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
#include <math.h>

#ifdef HGEXEC
#include "zoltan_mem.h"

typedef struct
  {
  int Proc;
  int Num_Proc;
  int Debug_Level;
  } ZZ;
  
#define ZOLTAN_DEBUG_LIST            8
#define ZOLTAN_DEBUG_ALL             10
#define ZOLTAN_TRACE_ENTER(a,b)      {strlen(b);}   /* NOP using variable to avoid compiler warning */
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
#include "hg_hypergraph.h"
#include "hg_util.h"

#define MIN(A,B)  (((A) < (B)) ? (A) : (B))
#define MAX(A,B)  (((A) > (B)) ? (A) : (B))

#define RANDOM_SEED   123456789   /* (time ((time_t*)NULL)) */
#define EPS           1e-6


/* Matching, Packing, and Grouping arrays.  Technically, they are all the same;
 * the different typedefs are not needed.  In the description below, Matching is
 * used; the same description applies to Packing and Grouping.  If a vertex i is
 * not being contracted with other vertices,  Matching[i] == i.  If vertices i,
 * j, and k are being contracted together to form one new vertex,
 * Matching[i] == j; Matching[j] == k; and Matching[k] == i;  The cycle describes
 * the contraction. */
typedef int *Matching;  /* length |V|, matching information of vertices */
typedef int *Packing;   /* length |V|, packing information of vertices */
typedef int *Grouping;  /* length |V|, grouping information of vertices */

typedef int *LevelMap;  /* length |V|, mapping of fine vertices onto coarse vertices */
typedef int *Partition; /* length |V|, partition ID for each vertex */

typedef struct {
   HGraph *hg;
   Partition part;
   Packing pack;
   LevelMap lmap;
} HGraphLevel;

/* Hypergraph Partitioning */
/* Function types for options to hypergraph partitioning */
struct HGPartParamsStruct;  /* Forward declaration */
typedef int ZOLTAN_HG_MATCHING_FN   (ZZ*, HGraph*, Matching, int*);
typedef int ZOLTAN_HG_PACKING_FN    (ZZ*, HGraph*, Packing,  int*);
typedef int ZOLTAN_HG_GROUPING_FN   (ZZ*, HGraph*, Grouping, int*);
typedef int ZOLTAN_HG_GLOBAL_PART_FN(ZZ*, HGraph*, int, Partition,
                                     struct HGPartParamsStruct *);
typedef int ZOLTAN_HG_LOCAL_REF_FN  (ZZ*, HGraph*, int, Partition,
                                     struct HGPartParamsStruct *, float);

/* Function types for edge-weight scaling functions. */
/* Placeholder for now; if these functions end up having the same argument */
/* list for each type, do not need separate types here or separate pointers */
/* in HGPartParams.  KDD */
typedef int ZOLTAN_HG_MATCHING_EWS_FN (ZZ*, Graph*);
typedef int ZOLTAN_HG_PACKING_EWS_FN  (ZZ*, HGraph*);
typedef int ZOLTAN_HG_GROUPING_EWS_FN (ZZ*, HGraph*);


/* Parameters to the hypergraph functions */
struct HGPartParamsStruct {
  float bal_tol;                        /* Balance tolerance in % of average */
  int kway;                             /* 1 -> kway, 0->recursive bisection */
  int redl;                             /* Reduction limit (constant). */

  char redm_str[MAX_PARAM_STRING_LEN];  /* Reduction method string. */
  ZOLTAN_HG_MATCHING_FN *matching;      /* Pointers to Matching, Packing and */
  ZOLTAN_HG_PACKING_FN  *packing;       /* Grouping fn specified by */
  ZOLTAN_HG_GROUPING_FN *grouping;      /* redm_str; NULL if not used */

  char redmo_str[MAX_PARAM_STRING_LEN]; /* Matching optimization string*/
  ZOLTAN_HG_MATCHING_FN *matching_opt;  /* Pointers to Matching, Packing and */
  ZOLTAN_HG_PACKING_FN  *packing_opt;   /* Grouping optimization fn specified*/
  ZOLTAN_HG_GROUPING_FN *grouping_opt;  /* by redmo_str; NULL if not used */

  int ews;                              /* Flag for Edge weight scaling */

  char global_str[MAX_PARAM_STRING_LEN];/* Global partitioning string and */
  ZOLTAN_HG_GLOBAL_PART_FN *global_part;/* pointer to Global partitioning fn */

  char local_str[MAX_PARAM_STRING_LEN]; /* Local refinement string and */
  ZOLTAN_HG_LOCAL_REF_FN *local_ref;    /* pointer to Local refinement fn */

  int check_graph;                      /* Flag indicating whether the input
                                           hypergraph should be checked for
                                           errors. */
  int output_level;                     /* Flag indicating amount of output
                                           from HG algorithms.  See levels
                                           HG_DEBUG_* below.  */
};
typedef struct HGPartParamsStruct HGPartParams;

/*
 * Hypergraph output levels.
 */
#define HG_DEBUG_NONE 0
#define HG_DEBUG_LIST 1
#define HG_DEBUG_ALL  2
#define HG_DEBUG_PRINT 3
#define HG_DEBUG_PLOT 4

int Zoltan_HG_rdivide (int, int, Partition, ZZ*, HGraph*, HGPartParams*, int);

int Zoltan_HG_Set_Part_Options  (ZZ*, HGPartParams*);
int Zoltan_HG_HPart_Lib    (ZZ*, HGraph*, int, Partition, HGPartParams*, int);
int Zoltan_HG_HPart_Info   (ZZ*, HGraph*, int, Partition, HGPartParams*);
double Zoltan_HG_hcut_size_total (HGraph*, Partition);
double Zoltan_HG_hcut_size_links (ZZ*, HGraph*, Partition);
double Zoltan_HG_HPart_balance   (ZZ*, HGraph*, int, Partition);

/* Scale Edge Weight */
int Zoltan_HG_Scale_HGraph_Weight (ZZ*, HGraph*, float*, int);
int Zoltan_HG_Scale_Graph_Weight  (ZZ*, Graph*,  float*, int);

/* Matching functions */
int Zoltan_HG_Matching (ZZ*, HGraph*, Matching, HGPartParams*, int*);
int Zoltan_HG_Set_Matching_Fn(HGPartParams*);

/* Packing */
int Zoltan_HG_Packing (ZZ*, HGraph*, Packing, HGPartParams*, int*);
int Zoltan_HG_Set_Packing_Fn(HGPartParams*);

/* Grouping */
int Zoltan_HG_Grouping (ZZ*, HGraph*, Packing, HGPartParams*, int*);
int Zoltan_HG_Set_Grouping_Fn(HGPartParams*);

/* Coarsening */
int Zoltan_HG_Coarsening   (ZZ*, HGraph*, Packing, HGraph*, int*);

/* Global Partitioning functions */
int Zoltan_HG_Global (ZZ*, HGraph*, int, Partition, HGPartParams*);
ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char*);

/* Local refinement functions */ /* KDD Placeholder for later. */
int Zoltan_HG_Local (ZZ*, HGraph*, int, Partition, HGPartParams*);
ZOLTAN_HG_LOCAL_REF_FN *Zoltan_HG_Set_Local_Ref_Fn(char*);

/* Sorting */
void Zoltan_quicksort_pointer_dec_float_int (int*, float*, int*, int, int);
void Zoltan_quicksort_pointer_dec_float     (int*, float*, int,  int);
void Zoltan_quicksort_pointer_inc_int_int   (int*, int*,   int*, int, int);
void Zoltan_quicksort_list_inc_int          (int*, int,    int);
void Zoltan_quicksort_pointer_inc_int_mult  (int*, int,    int,  int*, int*);

/* Heap datastructure */
typedef struct {
   int    space;
   int    n;
   int   *ele;
   int   *pos;
   float *value;
} HEAP;

#define Zoltan_HG_heap_empty(H)         (((H)->n)==0)
#define Zoltan_HG_heap_not_empty(H)     (((H)->n)!=0)
#define Zoltan_HG_heap_max_value(H)     ((H)->value[(H)->ele[0]])
#define Zoltan_HG_heap_peek_max(H)      ((H)->ele[0])
#define Zoltan_HG_heap_count(H)         ((H)->n)

int  Zoltan_HG_heap_init         (ZZ*, HEAP*, int);
void Zoltan_HG_heap_clear        (HEAP*);
void Zoltan_HG_heap_free         (HEAP*);
int  Zoltan_HG_heap_check        (HEAP*);
int  Zoltan_HG_heap_input        (HEAP*, int, float);
int  Zoltan_HG_heap_make         (HEAP*);
int  Zoltan_HG_heap_change_value (HEAP*, int, float);
int  Zoltan_HG_heap_extract_max  (HEAP*);
int  Zoltan_HG_heap_extract      (HEAP*, int);

int  Zoltan_HG_move_vertex (HGraph*, int, int, int, int*, int**, double*, HEAP*);
void Zoltan_HG_Plot(int, int, int, int*, int*, int*, char*);



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HYPERGRAPH_H */
