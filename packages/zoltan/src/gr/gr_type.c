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
static char *cvs_gr_type_c = "$Id$";
#endif

#include "all_const.h"
#include "gr_const.h"
#include "gr_type_const.h"
#include "gr_type.h"
#include "gr_tree_const.h"
#include "gr_hash_const.h"

/*****************************************************************************/
/*
 *  Functions that check value of BL_Graph_Type in determining what to do.
 */
/*****************************************************************************/

GRAPH *LB_new_graph()
{
/*
 *  Function to allocate and initialize a new graph data structure.
 *  A pointer to the new graph is returned.
 */
char *yo = "LB_new_graph";
GRAPH *graph;

  /*  
   *  For now, assuming size of GRAPH is same for both graph types. 
   */

  graph = (GRAPH *) LB_Malloc(sizeof(GRAPH));

  switch (BL_Graph_Type) {
  case GRAPH_BINARY_TREE:
    LB_initialize_tree_graph(graph);
    break;
  case GRAPH_HASH_TABLE:
    LB_initialize_hash_graph(graph);
    break;
  default:
    fprintf(stderr, "Error (%s):  Invalid graph type %d\n", yo, BL_Graph_Type);
    exit(-1);
  }

  return(graph);
}
