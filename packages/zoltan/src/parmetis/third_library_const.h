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

#ifdef indextype
#undef indextype
#endif

/* Include ParMetis and/or Scotch header files if necessary.
 * These include files must be available in the include path set in the
 * Zoltan configuration file.
 */
#ifdef ZOLTAN_PARMETIS
#include <parmetis.h>
#define indextype idxtype
#define weighttype idxtype
#endif /* ZOLTAN_PARMETIS */

#ifdef ZOLTAN_METIS
#include <metis.h>
#ifndef indextype
#define indextype idxtype
#define weighttype idxtype
#endif
#endif /* ZOLTAN_METIS */

#ifdef ZOLTAN_SCOTCH
#include <ptscotch.h>
#ifndef indextype
#define indextype SCOTCH_Num
#define weighttype SCOTCH_Num
#endif /*def indextype */
#endif /* ZOLTAN_SCOTCH */

#ifndef indextype
#define indextype int
#define weighttype int
#endif /* indextype */

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


/* Graph types */
#define NO_GRAPH     0
#define GLOBAL_GRAPH 1
#define LOCAL_GRAPH  2


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

extern int Zoltan_Build_Graph(ZZ *zz, int graph_type, int check_graph,
       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
       int obj_wgt_dim, int edge_wgt_dim,
       indextype **vtxdist, indextype **xadj, indextype **adjncy, float **ewgts,
       int **adjproc);

/* extern int Zoltan_Symmetrize_Graph( */
/*        const ZZ *zz, int graph_type, int check_graph, int num_obj, */
/*        ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, */
/*        int obj_wgt_dim, int edge_wgt_dim, */
/*        const indextype * const * vtxdist, indextype **xadj, indextype **adjncy, */
/*        float **ewgts, const int * const *adjproc); */

extern int Zoltan_Get_Num_Edges_Per_Obj(ZZ *, int, ZOLTAN_ID_PTR,
       ZOLTAN_ID_PTR, int **, int *, int *);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
