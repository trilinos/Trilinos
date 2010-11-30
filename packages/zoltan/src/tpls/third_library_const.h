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

#include "zoltan_util.h"

/*
 * "indextype" is the type used for global IDs and indices in the graph data structure.
 * "weighttype" is the type used for weights.  If there is only one third party library
 * linked in for graph algorithms, define indextype and weighttype to match the types
 * used by that library.
 *
 * If more than one such library is linked in, arbitrarily choose one.  At runtime if
 * either the indextype or weighttype is not compatible with the graph library API, then
 * print an error message.
 *
 * If there are no third party graph/ordering libraries then let indextype be ZOLTAN_ID_TYPE 
 * and let "weighttype" be float.
 *
 */

#define TPL_SCOTCH_DATATYPES   1
#define TPL_METIS_DATATYPES    2
#define TPL_ZOLTAN_DATATYPES   3

#undef TPL_USE_DATATYPE

#undef indextype
#undef weighttype

#ifdef ZOLTAN_PARMETIS
#define __parmetis__ 1
#else
#define __parmetis__ 0
#endif

#ifdef ZOLTAN_METIS
#define __metis__ 1
#else
#define __metis__ 0
#endif

#ifdef ZOLTAN_SCOTCH
#define __scotch__ 1
#else
#define __scotch__ 0
#endif

#ifdef ZOLTAN_PTSCOTCH
#define __ptscotch__ 1
#else
#define __ptscotch__ 0
#endif

#if __parmetis__ + __metis__ + __ptscotch__ + __scotch__ == 1

  #ifdef ZOLTAN_PARMETIS
    #include<parmetis.h>
    #define TPL_USE_DATATYPE TPL_METIS_DATATYPES
  #endif
  #ifdef ZOLTAN_METIS
    #include<metis.h>
    #define TPL_USE_DATATYPE TPL_METIS_DATATYPES
  #endif
  #ifdef ZOLTAN_SCOTCH
    #include<scotch.h>
    #define TPL_USE_DATATYPE TPL_SCOTCH_DATATYPES
  #endif
  #ifdef ZOLTAN_PTSCOTCH
    #include<ptscotch.h>
    #define TPL_USE_DATATYPE TPL_SCOTCH_DATATYPES
  #endif

#else

  #if __parmetis__ + __metis__ > 0
    #if __parmetis__ == 1
       #include<parmetis.h>
    #else
       #include<metis.h>
    #endif
    #if __parmetis__ + __metis__ ==  __parmetis__ + __metis__ + __ptscotch__ + __scotch__
      #define TPL_USE_DATATYPE TPL_METIS_DATATYPES
    #endif
  #endif

  #if __ptscotch__ + __scotch__ > 0
    #if __ptscotch__ == 1
       #include<ptscotch.h>
    #else
       #include<scotch.h>
    #endif
    #if __ptscotch__ + __scotch__ ==  __parmetis__ + __metis__ + __ptscotch__ + __scotch__
      #define TPL_USE_DATATYPE TPL_SCOTCH_DATATYPES
    #endif
  #endif

#endif

#ifndef TPL_USE_DATATYPE
  #if __parmetis__ + __metis__ + __ptscotch__ + __scotch__ == 0
    #define TPL_USE_DATATYPE TPL_ZOLTAN_DATATYPES
  #else
    #define TPL_USE_DATATYPE TPL_METIS_DATATYPES
  #endif
#endif

#if TPL_USE_DATATYPE == TPL_METIS_DATATYPES

  #define indextype idxtype
  #define weighttype idxtype

  #ifdef IDXTYPE_INT
    #define TPL_IDX_SPEC "%d"
    #define TPL_WGT_SPEC "%d"
  #else
    #define TPL_IDX_SPEC "%hd"
    #define TPL_WGT_SPEC "%hd"
  #endif

  #define TPL_INTEGRAL_WEIGHT

#else
  #if TPL_USE_DATATYPE == TPL_SCOTCH_DATATYPES

    #define indextype SCOTCH_Num
    #define weighttype SCOTCH_Num
    #define TPL_IDX_SPEC SCOTCH_NUMSTRING
    #define TPL_WGT_SPEC SCOTCH_NUMSTRING

    #define TPL_INTEGRAL_WEIGHT

  #else

    #define indextype ZOLTAN_ID_TYPE
    #define weighttype float
    #define TPL_IDX_SPEC ZOLTAN_ID_SPEC
    #define TPL_WGT_SPEC "%f"

    #define TPL_FLOAT_WEIGHT

  #endif

#endif

/* ParMETIS data types and definitions. */
/* #define IDXTYPE_IS_SHORT in order to use short as the idxtype.
 * Make sure these defs are consistent with those in your
 * ParMetis installation ! It is strongly recommended to use
 * integers, not shorts, if you load balance with weights.
*/

#ifdef IDXTYPE_IS_SHORT
/* typedef short idxtype; This should have been done in parmetis.h */
#define IDX_DATATYPE    MPI_SHORT
#define MAX_WGT_SUM (SHRT_MAX/8)
#else /* the default for idxtype is int; this is recommended */
/* typedef int idxtype; This should have been done in parmetis.h */
#define IDX_DATATYPE    MPI_INT
#define MAX_WGT_SUM (INT_MAX/8)
#endif


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
