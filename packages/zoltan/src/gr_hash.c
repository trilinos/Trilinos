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
static char *cvs_gr_hash_c = "$Id$";
#endif

#include "all_const.h"
#include "gr_const.h"
#include "gr_hash_const.h"
#include "id_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Functions implementing the graph as hash tables.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * Prototypes
 */

static GRAPH_ID_FN find_vertex_hash_graph;
static GRAPH_VERTEX_FN add_vertex_hash_graph;
static GRAPH_VERTEX_FN delete_vertex_hash_graph;
static GRAPH_FREE_FN free_hash_graph;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void initialize_hash_graph(GRAPH *graph)
{
/*
 *  Function that initializes function pointers for graphs implemented as 
 *  hash tables.
 */

  graph->Find_Vertex     = find_vertex_hash_graph;
  graph->Add_Vertex      = add_vertex_hash_graph;
  graph->Delete_Vertex   = delete_vertex_hash_graph;
  graph->Free_Graph      = free_hash_graph;
  graph->First_Vertex    = NULL;        /* KDD -- Need routines       */
  graph->Next_Vertex     = NULL;        /* KDD -- Need routines       */
  graph->Graph_Data      = NULL;        /* KDD -- Need data structure */
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *find_vertex_hash_graph(GRAPH *graph, ID *id)
{
/*
 *  Function to find a vertex in a graph.  Returns a pointer to the vertex.
 */

  return NULL;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void add_vertex_hash_graph(GRAPH *graph, VERTEX *vertex)
{
/*
 *  Function to add a vertex to the graph.
 */

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void delete_vertex_hash_graph(GRAPH *graph, VERTEX *vertex)
{
/*
 *  Function to remove a vertex from the graph.
 */

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void free_hash_graph(GRAPH **graph)
{

}
