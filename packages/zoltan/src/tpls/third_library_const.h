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


#ifndef __THIRD_LIBRARY_CONST_H
#define __THIRD_LIBRARY_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdint.h>
#include <float.h>
#include "zoltan_util.h"

/*****************************************************************************/
/* Include appropriate files for TPLs */
#ifdef ZOLTAN_METIS
  #include "metis.h"
  #define __metis__ 1
#else
  #define __metis__ 0
#endif

#ifdef ZOLTAN_PARMETIS
  #include "parmetis.h"
  #define __parmetis__ 1
#else
  #define __parmetis__ 0
#endif

#ifdef ZOLTAN_PTSCOTCH
  #include "ptscotch.h"
  #define __ptscotch__ 1
#else
  #define __ptscotch__ 0
#endif

#ifdef ZOLTAN_SCOTCH
  #ifndef ZOLTAN_PTSCOTCH
    #include "scotch.h"
  #endif
  #define __scotch__ 1
#else
  #define __scotch__ 0
#endif

/****************************************************************************/
/*  TPL-Specific settings for data types                                    */
/*
 * "indextype" is the type used for global numbers and indices in the graph
 *  data structure.
 * "weighttype" is the type used for weights.
 * 
 * If there are no third party graph/ordering libraries, let indextype be
 * ZOLTAN_GNO_TYPE and let "weighttype" be float.
 *
 * If there is only one third party library used for graph algorithms,
 * define indextype and weighttype to match the types used by that library.
 *
 * If more than one library is linked in, arbitrarily choose one.
 * Check for compatibility between the libraries here; all should use the same
 * size integer for indices.
 *
 * At runtime, if
 * either the indextype or weighttype is not compatible with the graph library
 * API, return an error.
 */

#define TPL_SCOTCH_DATATYPES   1
#define TPL_METIS_DATATYPES    2
#define TPL_ZOLTAN_DATATYPES   3

#undef TPL_USE_DATATYPE
#undef indextype
#undef weighttype

/* Select the data types to use */
#if __parmetis__ + __metis__ + __ptscotch__ + __scotch__ == 0
  /* No graph TPLs used; use Zoltan values */
  #define TPL_USE_DATATYPE TPL_ZOLTAN_DATATYPES
  #define indextype ZOLTAN_GNO_TYPE
  #define weighttype float
  #define realtype float
  #define TPL_FLOAT_WEIGHT
  #define MAX_WGT_SUM (FLT_MAX/8)
  #define TPL_IDX_SPEC ZOLTAN_GNO_SPEC
  #define TPL_WGT_SPEC "%f"
#elif (__ptscotch__ + __scotch__ > 0) && (__parmetis__ + __metis__  == 0)
  /* Using only Scotch/PTScotch */
  #define TPL_USE_DATATYPE TPL_SCOTCH_DATATYPES
#elif (__parmetis__ + __metis__ > 0) && (__ptscotch__ + __scotch__ == 0)
  /* Using only METIS/ParMETIS */
  #define TPL_USE_DATATYPE TPL_METIS_DATATYPES

#else
  /* Using both METIS/ParMETIS and Scotch/PTScotch; let METIS datatypes rule */
  #define TPL_USE_DATATYPE TPL_METIS_DATATYPES
#endif


#if TPL_USE_DATATYPE == TPL_METIS_DATATYPES

  #if PARMETIS_MAJOR_VERSION == 3
    /* Assume IDXTYPE_INT in ParMETIS v3.x */
    #ifndef IDXTYPE_INT
      /* typedef short idxtype; IDXTYPE_INT is not defined in parmetis.h */
      #error "ParMETIS short idxtype is not supported in Zoltan; define IDXTYPE_INT in parmetis.h."
    #endif
    #define indextype idxtype
    #define weighttype idxtype
    #define realtype float
    #define TPL_INTEGRAL_WEIGHT
    #define MAX_WGT_SUM (INT_MAX/8)
    #define TPL_IDX_SPEC "%d"
    #define TPL_WGT_SPEC "%d"
    #define IDXTYPEWIDTH 32

  #elif PARMETIS_MAJOR_VERSION == 4
    #define indextype idx_t
    #define weighttype idx_t
    #define realtype real_t
    #define TPL_INTEGRAL_WEIGHT
    #if IDXTYPEWIDTH == 32  /* defined in parmetis.h */
      #define MAX_WGT_SUM (INT32_MAX/8)
      #define TPL_IDX_SPEC "%d"
      #define TPL_WGT_SPEC "%d"
    #elif IDXTYPEWIDTH == 64 /* defined in parmetis.h */
      #define MAX_WGT_SUM (INT64_MAX/8)
      #define TPL_IDX_SPEC "%lld"
      #define TPL_WGT_SPEC "%lld"
    #endif
  #else
    #error "Unsupported version of ParMETIS; use ParMETIS 3.1 or 4."
  #endif

#elif TPL_USE_DATATYPE == TPL_SCOTCH_DATATYPES

  #define indextype SCOTCH_Num
  #define weighttype SCOTCH_Num
  #define realtype float
  #define MAX_WGT_SUM (SCOTCH_NUMMAX/8)
  #define TPL_IDX_SPEC SCOTCH_NUMSTRING
  #define TPL_WGT_SPEC SCOTCH_NUMSTRING
  #define TPL_INTEGRAL_WEIGHT

#endif

/**************************************************************************/


/* Graph types, used as mask to set bit in graph_type */
#define NO_GRAPH     0
#define LOCAL_GRAPH  1
#define TRY_FAST     2
#define FORCE_FAST   3
#define UNSYMMETRIC  4
  /* At this time, means A+At */
#define SYMMETRIZE   5

#define SET_NO_GRAPH(gtype) do { (*(gtype)) &= ~(1<<NO_GRAPH); (*(gtype)) &= ~(1<<LOCAL_GRAPH); } while (0)
#define SET_GLOBAL_GRAPH(gtype) do { (*(gtype)) &= ~(1<<LOCAL_GRAPH); (*(gtype)) &= ~(1<<NO_GRAPH); } while (0)
#define SET_LOCAL_GRAPH(gtype) do { (*(gtype)) |= (1<<LOCAL_GRAPH); (*(gtype)) &= ~(1<<NO_GRAPH); } while (0)
#define IS_NO_GRAPH(gtype) ((!((gtype)&(1<<LOCAL_GRAPH))) && (((gtype)&(1<<NO_GRAPH))))
#define IS_GLOBAL_GRAPH(gtype) ((!((gtype)&(1<<NO_GRAPH))) && (!((gtype)&(1<<LOCAL_GRAPH))))
#define IS_LOCAL_GRAPH(gtype) ((!((gtype)&(1<<NO_GRAPH))) && (((gtype)&(1<<LOCAL_GRAPH))))


/* Misc. defs to be used with MPI */
#define TAG1  32001
#define TAG2  32002
#define TAG3  32003
#define TAG4  32004
#define TAG5  32005
#define TAG6  32006
#define TAG7  32007


/* Zoltan function prototypes */
extern int Zoltan_Graph_Package_Set_Param(char *, char *);
#ifdef ZOLTAN_PARMETIS
extern int Zoltan_ParMetis_Set_Param(char *, char *);
#endif /* ZOLTAN_PARMETIS */
#ifdef ZOLTAN_SCOTCH
extern int Zoltan_Scotch_Set_Param(char *, char *);
#endif /* ZOLTAN_SCOTCH */
extern int Zoltan_Third_Set_Param(char *, char *);

extern int Zoltan_Build_Graph(struct Zoltan_Struct *zz, int *graph_type, int check_graph,
       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
       int obj_wgt_dim, int * edge_wgt_dim,
       ZOLTAN_GNO_TYPE **vtxdist, int **xadj, ZOLTAN_GNO_TYPE **adjncy, float **ewgts,
       int **adjproc);

extern int Zoltan_Get_Num_Edges_Per_Obj(struct Zoltan_Struct *, int, ZOLTAN_ID_PTR,
       ZOLTAN_ID_PTR, int **, int *, int *);


/*==========================================================================
 * The ZOS structure copied from order/order_const.h, but using TPL datatypes.
 */

struct Zoltan_TPL_Order_Struct {
  indextype needfree;
  indextype nbr_objects;              /* # of objects (local) */
  ZOLTAN_ID_PTR gids;           /* ptr to list of global ids */
  ZOLTAN_ID_PTR lids;           /* ptr to list of local ids */
  indextype *rank;            /* rank[i] is the rank of gids[i] */
  ZOLTAN_ID_PTR gidrank;
  indextype *iperm;
  indextype  start_index;
  char method[MAX_PARAM_STRING_LEN+1]; /* Ordering method used */
  char order_type[MAX_PARAM_STRING_LEN+1]; /* Ordering method used */

  /* Elimination Tree */
  indextype nbr_blocks;               /* Out: number of ordering blocks */
  indextype *start;                   /* Out: start[i] is the first vertex of block i */
  indextype *ancestor;                /* Out: father of block i */
  indextype *leaves;                  /* Out: list of all leaves */
  indextype nbr_leaves;               /* Number of leaves */

  indextype *vtxdist;                 /* How vertices are distributed accross processors */

  /* Deprecated */
  indextype  num_separators;          /* Optional: # of separators. */
  indextype *sep_sizes;               /* Optional: Separator sizes. */
};

typedef struct Zoltan_TPL_Order_Struct ZTPL_OS;

int  Zoltan_TPL_Order_Init_Tree (struct Zoltan_TPL_Order_Struct *order, indextype blocknbr, indextype leavesnbr);
void Zoltan_TPL_Order_Free_Struct(struct Zoltan_TPL_Order_Struct *order);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
