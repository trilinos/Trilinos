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
static char *cvs_gr_c = "$Id$";
#endif

#include "all_const.h"
#include "gr_const.h"
#include "id_const.h"
#include "id_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

GRAPH *build_graph_ala_chaco(
  int num_vertices,          /* Number of vertices in the processor's graph. */
  int vertex_weights[],      /* Array of vertex weights (ordered in same 
                                order as vertex_list).                       */
  int start_nbor_list[],     /* Array of indices of first nbor in nbor_list.
                                (ordered in same order as vertex_list)       */
  ID nbor_list[],            /* Array of neighboring vertex IDs.  The nbors
                                of vertex i are stored in nbor_list[sum_i] to
                                nbor_list[sum_i+num_nbors[i]-1], where
                                sum_i = SUM(num_nbors[j], j=0 to i-1).       */
  float edge_weights[],      /* Array of edge_weights for each edge of a 
                                vertex.  The weights for edges from vertex i
                                are stored in edge_weights[sum_i] to 
                                edge_weights[sum_i+num_nbors[i]-1], where
                                sum_i = SUM(num_nbors[j], j=0 to i-1).
                                The order of the edge_weights in this 
                                sub_vector corresponds with the order of
                                vertex neighbors in nbor_list.               */
  int num_dim,               /* Number of dimensions in the geometry; 
                                num_dim == 0 if geometry information will 
                                not be used and stored.                      */
  float *x,                  /* x-coordinates of the vertices (if using
                                geometry info; NULL otherwise).              */
  float *y,                  /* y-coordinates of the vertices (if using 
                                geometry info; NULL otherwise).              */
  float *z                   /* z-coordinates of the vertices (if using
                                geometry info; NULL otherwise).              */
)
{
/*
 *  Routine to build the graph of a processor's vertices.
 */

int i; 
int start_i;
double coor[3] = {0., 0., 0.};
VERTEX *vertex;
GRAPH *graph;
ID vertex_ID;

  graph = new_graph();


  for (i = 0; i < num_vertices; i++) {

    /*
     * Chaco numbers vertices from 1 to num_vertices; 
     * increment i for building ID.
     */

    BL_ID_Util.New_ID(&vertex_ID, i+1);

    start_i = start_nbor_list[i];

    if (num_dim > 0) {
      coor[0] = x[i];
      if (num_dim > 1) {
        coor[1] = y[i];
        if (num_dim > 2) {
          coor[2] = z[i];
        }
      }
    }

    vertex = new_vertex(&vertex_ID, 
                 (vertex_weights != NULL ?  vertex_weights[i] : 1),
                 start_nbor_list[i+1]-start_i, &(nbor_list[start_i]), 
                 (edge_weights != NULL ? &(edge_weights[start_i]) : NULL), 
                  num_dim, coor);

    graph->Add_Vertex(graph, vertex);
  }

  return(graph);
}
