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

#ifndef __PARMETIS_JOSTLE_CONST_H
#define __PARMETIS_JOSTLE_CONST_H

#include <limits.h>
#include "comm_const.h"

/* Do we have the multiconstraint beta version of ParMetis? */
/* #define BETA_PARMETIS */

/* Include ParMetis and/or Jostle header files if necessary. 
 * These include files must be available in the include path set in the 
 * Zoltan configuration file.
 */
#ifdef LB_PARMETIS
#include "parmetis.h"
#else 
typedef int idxtype; 
#endif

#ifdef LB_JOSTLE
#include "jostle.h"
#endif

/* ParMetis option defs. These must be identical to the defs
 * in defs.h in the version of ParMetis you are using!
 */
#define OPTION_IPART            1
#define OPTION_FOLDF            2
#define OPTION_DBGLVL           3
#define MAX_OPTIONS             4  /* Total number of options +1 */

/* Misc. defs to be used with MPI */
#define TAG1  32001
#define TAG2  32002
#define TAG3  32003
#define TAG4  32004
#define TAG5  32005

/* Misc. local constants */
#define CHUNKSIZE 20  /* Number of nodes to allocate in one chunk. */

/* Data structures used in ParMetis interface routines */

/* An array of this data structure works with a parallel array of 
 * LB_ID_PTR called proc_list_nbor containing the global IDs of the
 * neighboring object. 
 * This separate array is needed to prevent individual mallocs of
 * neighboring global IDs.
 */
struct LB_edge_info {
  LB_ID_PTR my_gid;  /* Pointer to the Global id of local vtx */
  int my_gno;        /* Global number of local vtx */
  int nbor_proc;     /* Proc id for the neighboring proc */
  int *adj;          /* Pointer to adjcny array */
};

struct LB_hash_node {
  LB_ID_PTR gid;     /* Pointer to a Global id */
  int gno;           /* Global number */
  struct LB_hash_node * next;
};


/* ParMETIS data types and definitions. */

/* Undefine the following #define in order to use short as the idxtype.
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

/* Zoltan function prototypes */
extern int LB_Set_ParMetis_Param(char *, char *);
extern int LB_Set_Jostle_Param(char *, char *);
extern int LB_Verify_Graph(MPI_Comm comm, idxtype *vtxdist, idxtype *xadj, 
              idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
              int vwgt_dim, int ewgt_dim, int check_graph, int debug_level);
extern int LB_Scatter_Graph(int have_graph, idxtype **vtxdist, idxtype **xadj, idxtype **adjncy,
              idxtype **vwgt, idxtype **adjwgt, float   **xyz, int     ndims,
              LB      *lb, struct Comm_Obj **plan);


#endif
