/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#ifndef lint
static char *cvs_ch_driver_h = "$Id$";
#endif

/* Prototypes */

LB_NUM_OBJ_FN  num_obj_fn;
LB_OBJ_LIST_FN  list_obj_fn;
LB_OBJ_WEIGHT_FN  obj_weight_fn;
LB_NUM_EDGES_FN  num_edges_fn;
LB_EDGE_LIST_FN  edge_list_fn;


/* Data structures */

struct GraphStruct {
  int lnvtxs;     /* # of vertices on local proc */
  int *vtxdist; 
  int *start;
  int *adjncy;
  int *vwgts;
  float *ewgts;
  int nprocs;
  int myproc;
};

typedef struct GraphStruct Graph;
