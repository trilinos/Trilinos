
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
static char *cvs_all_main = "$Id$";
#endif

#include "all_const.h"
#include "gr_const.h"
#include "id_const.h"
#include "id_util_const.h"
#include "ch_input_const.h"
#include "ch_input.h"
#include "gr_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Functions to run the dynamic load-balancer in stand-alone mode.
 *  The input files read are in Chaco format.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_chaco_graph_stuff(char *, int *, int **, ID **, int **, 
                                   float **,
                                   int *, float **, float **, float **); 

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

main(int *argv, char *argc[])
{
int num_vertices;             /* Number of vertices in problem to be solved  */
ID *edge_list;                /* List of edges for vertices of graph.        */
int *start_edge_list;         /* Array of ptrs to first edge for each vertex */
int *vertex_weights;          /* Array of vertex weights for each vertex.    */
float *edge_weights;          /* Array of edge weights for each edge of 
                                 each vertex.  Ordered same as edge_list.    */
int num_dim;                  /* Dimensionality of the problem.              */
float *x, *y, *z;             /* Geometric coordinates for each vertex.      */
GRAPH *graph;                 /* Data structure for graph to be processed.   */


  read_chaco_graph_stuff(argc[1], &num_vertices, &start_edge_list, &edge_list,
                   &vertex_weights, &edge_weights, &num_dim, &x, &y, &z);

  graph = LB_build_graph_ala_chaco(num_vertices, vertex_weights, 
                                start_edge_list, edge_list, 
                                edge_weights, num_dim, x, y, z);

  LB_print_graph(graph);

  graph->Free_Graph(&graph);

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_chaco_graph_stuff(
  char *name,                  /* Base name of chaco input files             */
  int *num_vertices,           /* Number of vertices read.                   */
  int **start_edge_list,       /* Array of indices into edge list; ordered
                                  by vertex number.                          */
  ID **edge_list,              /* Array of neighbor nodes for all vertices.  */
  int **vertex_weights,        /* Array of vertex weights; ordered by vertex 
                                  number.                                    */
  float **edge_weights,        /* Array of edge weights of each edge; 
                                  ordered in same way as edge_list.          */
  int *num_dim,                /* Number of dimensions in the mesh.          */
  float **x,                   /* x-coordinates for each vertex; ordered by
                                  vertex number.                             */
  float **y,                   /* y-coordinates for each vertex; ordered by
                                  vertex number.                             */
  float **z                    /* z-coordinates for each vertex; ordered by
                                  vertex number.                             */
) 
{
char filename[30];
FILE *fp;
int *edge_list_ints;
int i;
int num_edges;

  /* 
   *  Read the Chaco graph file.
   */

  printf("Reading Chaco graph file...%s\n", name);
  sprintf(filename, "%s.graph", name);
  fp = fopen(filename, "r");
  LB_chaco_input_graph(fp, filename, start_edge_list, &edge_list_ints,
                    num_vertices, vertex_weights, edge_weights);
  fclose(fp);

  /*
   *  Convert from Chaco's integer ids to ID types.
   */

  num_edges = (*start_edge_list)[*num_vertices];
  *edge_list = (ID *) Array_Alloc(1, num_edges, sizeof(ID));
  for (i = 0; i < num_edges; i++)
    BL_ID_Util.New_ID(&((*edge_list)[i]), edge_list_ints[i]);
  LB_FREE(&edge_list_ints);

  /*
   *  Read the Chaco geometry file.
   */

  printf("Reading Chaco geometry file...%s \n", name);
  sprintf(filename, "%s.coords", name);
  fp = fopen(filename, "r");
  LB_chaco_input_geom(fp, filename, *num_vertices, num_dim, x, y, z);

  fclose(fp);
}
