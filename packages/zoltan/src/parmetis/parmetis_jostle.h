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

#ifndef __PARMETIS_JOSTLE_H
#define __PARMETIS_JOSTLE_H

#include <limits.h>
#include "zoltan_comm.h"
#include "parmetis_jostle_const.h"

/* Include ParMetis and/or Jostle header files if necessary. 
 * These include files must be available in the include path set in the 
 * Zoltan configuration file.
 */
#ifdef ZOLTAN_PARMETIS
#include "parmetis.h"
#else 
typedef int idxtype; 
#endif

#ifdef ZOLTAN_JOSTLE
#include "jostle.h"
#endif

/* Guess the version number of ParMetis if not defined */
/* PARMETIS_MAJOR_VERSION is only defined in version 3.0 and higher */
#if (!defined(PARMETIS_MAJOR_VERSION))
#define PARMETIS_MAJOR_VERSION 2
#define PARMETIS_MINOR_VERSION 0
#endif 

/* ParMetis option defs. These must be identical to the defs
 * in defs.h in the version of ParMetis you are using!
 * Both ParMetis 2.0 and 3.0 defs are included below.
 */
#define OPTION_IPART            1
#define OPTION_FOLDF            2
#define OPTION_DBGLVL           3
#define PMV3_OPTION_DBGLVL      1
#define PMV3_OPTION_SEED        2
#define PMV3_OPTION_IPART       3
#define PMV3_OPTION_PSR         3
#define PMV3_OPT_USE_OBJ_SIZE   9  /* Added by EB, not in ParMetis */
#define MAX_OPTIONS             10 /* Max number of options +1 */
/* Other ParMetis constants we may need */
#define GLOBAL_DBGLVL		0  /* Default debug level */
#define GLOBAL_SEED		15 /* Default random seed */
#define COUPLED                 1  /* Processors coupled to partitions? */
#define DISCOUPLED              2  /* Processors coupled to partitions? */


/* Misc. defs to be used with MPI */
#define TAG1  32001
#define TAG2  32002
#define TAG3  32003
#define TAG4  32004
#define TAG5  32005
#define TAG6  32006
#define TAG7  32007

/* Misc. local constants */
#define CHUNKSIZE 20  /* Number of nodes to allocate in initial chunk. */
#define REALLOC_FACTOR 1.5  /* Increase size by this factor if too small. */

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


/* Data structures used in ParMetis interface routines */
/* An array of this data structure works with a parallel array of
 * ZOLTAN_ID_PTR called proc_list_nbor containing the global IDs of the
 * neighboring object.
 * This separate array is needed to prevent individual mallocs of
 * neighboring global IDs.
 */
struct Edge_Info {
  ZOLTAN_ID_PTR my_gid;  /* Pointer to the Global id of local vtx */
  int my_gno;        /* Global number of local vtx */
  int nbor_proc;     /* Proc id for the neighboring proc */
  int *adj;          /* Pointer to adjcny array */
};

struct Hash_Node {
  ZOLTAN_ID_PTR gid;     /* Pointer to a Global id */
  int gno;           /* Global number */
  struct Hash_Node * next;
};


/* Extern function prototypes. Should be in a separate header file? */
extern int Zoltan_Verify_Graph(MPI_Comm comm, idxtype *vtxdist, idxtype *xadj,
       idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
       int vwgt_dim, int ewgt_dim, int check_graph, int debug_level);
extern int Zoltan_Scatter_Graph(idxtype **vtxdist, idxtype **xadj, 
       idxtype **adjncy, idxtype **vwgt, idxtype **vsize, idxtype **adjwgt, 
       float **xyz, int ndims, ZZ *zz, ZOLTAN_COMM_OBJ **plan);
extern int Zoltan_Build_Graph( ZZ *zz, int get_graph, int check_graph,
       ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
       int obj_wgt_dim, int edge_wgt_dim,
       idxtype **vtxdist, idxtype **xadj, idxtype **adjncy, 
       float **ewgts);



#endif
