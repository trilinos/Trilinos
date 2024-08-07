// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stddef.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_err_const.h"
#include "ch_init_dist_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*
 *  Routines to perform initial data distributions for Chaco files.
 *  To add a new initial distribution, add a new INITIAL_* value to
 *  ch_init_dist_const.h and add code handling the new case to each
 *  of the switch statements in this file.
 *
 *  Methods currently implemented:
 *    INITIAL_FILE:    takes assignments from user input file 
 *                     (graphname.assign).
 *    INITIAL_LINEAR:  assigns the first n/p vertices to proc 0, the
 *                     next n/p vertices to proc 1, etc.
 *    INITIAL_OWNER:   same as INITIAL_LINEAR; takes on different meaning 
 *                     only for hyperedges.
 *    INITIAL_CYCLIC:  deals out the vertices like cards; i.e., proc 0
 *                     gets the first vertex, proc 1 gets the second vertex,
 *                     ..., proc p-1 gets the pth vertex, proc 0 gets the
 *                     (p+1)th vertex, proc 1 gets the (p+2)th vertex, etc.
 *
 *  It is possible to distribute objects to a subset of the
 *  the procs. The variable Num_Proc_Dist indicates that only
 *  procs 0 through (Num_Proc_Dist-1) should be assigned any
 *  vertices, the rest will have no objects.
 */

/*****************************************************************************/
/* Global variables for all routines in this file.  These variables 
 * are set in ch_dist_init and are used in all subsequent queries of
 * the chaco distribution.
 */
static ZOLTAN_GNO_TYPE Gnvtxs;     /* Global number of vertices (across all procs)  */
static int Num_Proc;        /* Global number of processors                   */
static int Num_Proc_Dist;   /* Number of processors to distribute data; 
                               Num_Proc_Dist may be < Num_Proc!              */
static int Initial_Method;  /* Flag indicating which initial decomposition
                               method to use.                                */

/*****************************************************************************/
/* Data structures specifically for INITIAL_LINEAR */
static ZOLTAN_GNO_TYPE *Vtxdist;   /* Array of size Num_Proc+1 indicating range of
                               vertices assigned to a processor.  It is 
                               assumed that vertices are numbered 
                               consecutively.  The values stored in Vtxdist
                               assume vertex IDs are zero-based (e.g., the
                               lowest-numbered vertex has ID 0).
                               Processor p is assigned vertices Vtxdist[p]
                               through (Vtxdist[p+1]-1), inclusive.          */

static void create_Vtxdist();

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ch_dist_init(int num_proc, int gnvtxs, PARIO_INFO_PTR pio_info,
                  short **assignments, int host_proc, MPI_Comm comm)
{
/* Routine to initialize the decomposition.   This routine sets the
 * global variables to be used in the query routines included in this
 * file.  It also performs any initializations that a particular method
 * may need.
 * All processors call this routine prior to the actual distribution of 
 * graph data.
 */
int proc, i;
int max_assignment, have_assignments;

  Num_Proc = num_proc; 
  Gnvtxs = gnvtxs;
  Num_Proc_Dist  = pio_info->init_dist_procs;
  /* If Num_Proc_Dist has an invalid value, set it to Num_Proc by default */ 
  if ((Num_Proc_Dist <= 0) || (Num_Proc_Dist > Num_Proc)) 
     Num_Proc_Dist = Num_Proc;

  /* Broadcast initial assignments if they exist.  
   * Assignments can be used for partitions and/or processors.
   */
  MPI_Comm_rank(zoltan_get_global_comm(), &proc);

  /* First, tell other processors whether the assignments array is NULL. */
  if (proc == host_proc) 
    have_assignments = (*assignments != NULL);
  MPI_Bcast(&have_assignments, 1, MPI_INT, host_proc, comm);

  /* Now send the assignments if they exists. */
  if (have_assignments) {
    /* Allocate and broadcast the assignments array. */
    if (proc != host_proc) {
      *assignments = (short *) malloc(gnvtxs * sizeof(short));
      if (*assignments == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return;
      }
    }
    MPI_Bcast( *assignments, gnvtxs, MPI_SHORT, host_proc, comm);

    for (max_assignment = 0, i = 0; i < gnvtxs; i++)
      if ((*assignments)[i] > max_assignment)
        max_assignment = (*assignments)[i];
    if (max_assignment >= Num_Proc) {
      /* Highest partition number in assignment array is greater than
       * the number of processors.  Assume the assignments are to be
       * used ONLY as partition assignments; use INITIAL_LINEAR
       * assignment to put vertices on processors.
       */
      pio_info->init_dist_type = INITIAL_LINEAR;
    }
  }

  Initial_Method = pio_info->init_dist_type;
  switch(Initial_Method) {
  case INITIAL_LINEAR:
  case INITIAL_OWNER:
    /* create the vtxdist array. */
    create_Vtxdist();
    break;
  case INITIAL_FILE:
  case INITIAL_CYCLIC:
    /* no initialization needed for this method. */
    break;
  default: 
    Gen_Error(0, "Invalid Initial Distribution Type in ch_dist_init");
    error_report(proc);
    return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ch_dist_num_vtx(
  int target_proc,
  short *assignments
)
{
/* Function that returns the number of vertices assigned to processor
 * target_proc.
 */
int num;

  switch(Initial_Method) {
  case INITIAL_FILE: {
    int i;
    num = 0;
    for (i = 0; i < Gnvtxs; i++)
      if (assignments[i] == target_proc) num++;
    break;
  }
  case INITIAL_LINEAR:
  case INITIAL_OWNER:
    num = Vtxdist[target_proc+1] - Vtxdist[target_proc];
    break;
  case INITIAL_CYCLIC:
    if (target_proc < Num_Proc_Dist){
      num = Gnvtxs / Num_Proc_Dist;
      if ((Gnvtxs % Num_Proc_Dist) > target_proc)
        num++;
    }
    else
      num = 0;
    break;
  default:
    Gen_Error(0, "Invalid Initial Distribution Type in ch_dist_num_vtx");
    return(-1);
  }
  return(num);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ch_dist_max_num_vtx(short *assignments)
{
/* Function that returns the maximum number of vertices assigned to any 
 * processor.
 */
int i;
int tmp, max = 0;

  for (i = 0; i < Num_Proc; i++)
    if ((tmp = ch_dist_num_vtx(i,assignments)) > max) max = tmp;

  return max;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ch_dist_vtx_list(
  int *vtx_list,
  int *nvtx,
  int target_proc,
  short *assignments          
)
{
/* Function that returns a list of vertices assigned to proc target_proc.
 * The list is returned in vtx_list.  The function assumes vertex ID
 * numbering is zero-based; e.g., the lowest-numbered vertex has ID 0.
 * This convention allows the list entries to be used as array indices
 * in the arrays containing chaco input file data.
 * The number of entries in the list is returned in nvtx.
 */
int i, proc;

  *nvtx = 0;
  proc = (ZOLTAN_GNO_TYPE)target_proc;

  switch(Initial_Method) {
  case INITIAL_FILE:
    for (i = 0; i < Gnvtxs; i++)
      if (assignments[i] == target_proc)
        vtx_list[(*nvtx)++] = i;
    break;
  case INITIAL_LINEAR:
  case INITIAL_OWNER:
    for (i = Vtxdist[target_proc]; i < Vtxdist[target_proc+1]; i++)
      vtx_list[(*nvtx)++] = i;
    break;
  case INITIAL_CYCLIC:
    if (target_proc < Num_Proc_Dist){
      for (i = proc; i < Gnvtxs; i+=Num_Proc_Dist) 
        vtx_list[(*nvtx)++] = i;
    }
    break;
  default:
    Gen_Error(0, "Invalid Initial Distribution Type in ch_dist_vtx_list");
    return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ch_dist_proc(int v, short *assignments, int base)
{
/* Function that returns the processor to which a vertex v is assigned.
 * The function assumes the vertex numbering is "base"-based (i.e., lowest 
 * numbered vertex is vertex base); this convention is used since the function
 * is called primarily to find adjacent vertices' processor assignments and
 * the read-in adjacency information is "base"-based.
 * E.g., for Chaco input files, base == 1; vertex comparisons are one-based.
 *       for HG input files, base may be 0 or 1.
 */
int p;
ZOLTAN_GNO_TYPE b = (ZOLTAN_GNO_TYPE)base;

  switch(Initial_Method) {
  case INITIAL_FILE:
    /* return the appropriate entry from the assignments array. */
    p = assignments[v-b];
    break;
  case INITIAL_LINEAR:
  case INITIAL_OWNER:
    for (p = 0; p < Num_Proc_Dist; p++)
      /* Since v is "base"-based and Vtxdist is 0-based, add base to 
       * Vtxdist[p+1]. */
      if (v < Vtxdist[p+1]+b) break;
    break;
  case INITIAL_CYCLIC:
    /* test for (v-base) as v is "base"-based and 
     * INITIAL_CYCLIC equations are 0-based */
    p = (v-b) % Num_Proc_Dist;
    break;
  default:
    Gen_Error(0, "Invalid Initial Distribution Type in ch_dist_proc");
    return -1;
  }
  return p;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void create_Vtxdist()
{
/* Function that creates Vtxdist for the INITIAL_LINEAR distribution.
 * Vtxdist is an array of size Num_Proc+1 indicating range of
 * vertices assigned to a processor.  It is assumed that vertices are numbered 
 * consecutively.  The values stored in Vtxdist assume vertex IDs are 0-based 
 * (e.g., the lowest-numbered vertex has ID 0).  Proc p is assigned vertices 
 * Vtxdist[p] through (Vtxdist[p+1]-1), inclusive.         
 */
int rest, i, n;

  /* Set up Vtxdist data */
  if (Vtxdist == NULL){
    Vtxdist = (ZOLTAN_GNO_TYPE *) malloc((Num_Proc+1) * sizeof(ZOLTAN_GNO_TYPE));
    if (Vtxdist == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return;
    }
  }
  /* Calculate uniform vertex distribution */
  Vtxdist[0] = 0;
  rest = Gnvtxs;
  for (i=0; i<Num_Proc_Dist; i++){
    n = rest/(Num_Proc_Dist-i);
    Vtxdist[i+1] = Vtxdist[i] + n;
    rest -= n;
  }
  /* Procs >= Num_Proc_Dist get no vertices */
  for (i=Num_Proc_Dist; i<Num_Proc; i++){
    Vtxdist[i+1] = Vtxdist[i];
  }

}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
