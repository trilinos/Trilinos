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
 * $Name$
 *====================================================================*/

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_err_const.h"
#include "ch_init_dist_const.h"

/*
 *  Routines to perform initial data distributions for Chaco files.
 *  To add a new initial distribution, add a new INITIAL_* value to
 *  ch_init_dist_const.h and add code handling the new case to each
 *  of the switch statements in this file.
 *
 *  Methods currently implemented:
 *    INITIAL_LINEAR:  assigns the first n/p vertices to proc 0, the
 *                     next n/p vertices to proc 1, etc.
 *    INITIAL_CYCLIC:    deals out the vertices like cards; i.e., proc 0
 *                     gets the first vertex, proc 1 gets the second vertex,
 *                     ..., proc p-1 gets the pth vertex, proc 0 gets the
 *                     (p+1)th vertex, proc 1 gets the (p+2)th vertex, etc.
 */

/*****************************************************************************/
/* Global variables for all routines in this file.  These variables 
 * are set in ch_dist_init and are used in all subsequent queries of
 * the chaco distribution.
 */
static int Gnvtxs;          /* Global number of vertices (across all procs)  */
static int Num_Proc;        /* Global number of processors                   */
static int Initial_Method;  /* Flag indicating which initial decomposition
                               method to use.                                */

/*****************************************************************************/
/* Data structures specifically for INITIAL_LINEAR */
static int *Vtxdist;        /* Array of size Num_Proc+1 indicating range of
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

void ch_dist_init(int num_proc, int gnvtxs, PARIO_INFO_PTR pio_info)
{
/* Routine to initialize the decomposition.   This routine sets the
 * global variables to be used in the query routines included in this
 * file.  It also performs any initializations that a particular method
 * may need.
 * All processors call this routine prior to the actual distribution of 
 * graph data.
 */

  Num_Proc = num_proc;
  Gnvtxs = gnvtxs;
  Initial_Method = pio_info->init_dist_type;

  switch(Initial_Method) {
  case INITIAL_LINEAR:
    /* create the vtxdist array. */
    create_Vtxdist();
    break;
  case INITIAL_CYCLIC:
    /* no initialization needed for this method. */
    break;
  default: {
    int proc;
    Gen_Error(0, "Invalid Initial Distribution Type in ch_dist_init");
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    error_report(proc);
    return;
  }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ch_dist_num_vtx(
  int target_proc
)
{
/* Function that returns the number of vertices assigned to processor
 * target_proc.
 */
int num;

  switch(Initial_Method) {
  case INITIAL_LINEAR:
    num = Vtxdist[target_proc+1] - Vtxdist[target_proc];
    break;
  case INITIAL_CYCLIC:
    num = Gnvtxs / Num_Proc;
    if ((Gnvtxs % Num_Proc) > target_proc)
      num++;
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

int ch_dist_max_num_vtx()
{
/* Function that returns the maximum number of vertices assigned to any 
 * processor.
 */
int i;
int tmp, max = 0;

  for (i = 0; i < Num_Proc; i++)
    if ((tmp = ch_dist_num_vtx(i)) > max) max = tmp;

  return max;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ch_dist_vtx_list(
  int *vtx_list,
  int *nvtx,
  int target_proc
)
{
/* Function that returns a list of vertices assigned to proc target_proc.
 * The list is returned in vtx_list.  The function assumes vertex ID
 * numbering is zero-based; e.g., the lowest-numbered vertex has ID 0.
 * This convention allows the list entries to be used as array indices
 * in the arrays containing chaco input file data.
 * The number of entries in the list is returned in nvtx.
 */
int i;

  *nvtx = 0;

  switch(Initial_Method) {
  case INITIAL_LINEAR:
    for (i = Vtxdist[target_proc]; i < Vtxdist[target_proc+1]; i++)
      vtx_list[(*nvtx)++] = i;
    break;
  case INITIAL_CYCLIC:
    for (i = target_proc; i < Gnvtxs; i+=Num_Proc) 
        vtx_list[(*nvtx)++] = i;
    break;
  default:
    Gen_Error(0, "Invalid Initial Distribution Type in ch_dist_vtx_list");
    return;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ch_dist_proc(int v)
{
/* Function that returns the processor to which a vertex v is assigned.
 * The function assumes the vertex numbering is one-based (i.e., lowest 
 * numbered vertex is vertex 1); this convention is used since the function
 * is called primarily to find adjacent vertices' processor assignments and
 * the read-in adjacency information is one-based.
 */
int p;

  switch(Initial_Method) {
  case INITIAL_LINEAR:
    for (p = 0; p < Num_Proc; p++)
      /* Compare with <= since v is 1-based and Vtxdist is 0-based. */
      if (v <= Vtxdist[p+1]) break;
    break;
  case INITIAL_CYCLIC:
    /* test for (v-1) as v is 1-based and INITIAL_CYCLIC equations are 0-based */
    p = (v-1) % Num_Proc;
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
    Vtxdist = (int *) malloc((Num_Proc+1) * sizeof(int));
    if (Vtxdist == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return;
    }
  }
  /* Calculate uniform vertex distribution */
  (Vtxdist)[0] = 0;
  rest = Gnvtxs;
  for (i=0; i<Num_Proc; i++){
    n = rest/(Num_Proc-i);
    (Vtxdist)[i+1] = (Vtxdist)[i] + n;
    rest -= n;
  }
}
