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
#ifndef __GR_CONST_H
#define __GR_CONST_H

#ifndef lint
static char *cvs_gr_const_h = "$Id$";
#endif

typedef void GRAPH_VERTEX_FN(GRAPH*, VERTEX *);
typedef VERTEX *GRAPH_ID_FN(GRAPH*, ID *);
typedef void GRAPH_FREE_FN(GRAPH **);
typedef VERTEX *GRAPH_ITERATOR_FN(GRAPH *, LOOP_CONTROL *);

struct Graph_Struct {

  /* PUBLIC */

  GRAPH_ID_FN *Find_Vertex;            /* Routine to find a vertex with a 
                                          given ID in the graph.             */
  GRAPH_VERTEX_FN *Add_Vertex;         /* Routine to add a given vertex to
                                          the graph.                         */
  GRAPH_VERTEX_FN *Delete_Vertex;      /* Routine to delete a given vertex
                                          from the graph.                    */
  GRAPH_FREE_FN *Free_Graph;           /* Routine to deallocate a graph.     */
  GRAPH_ITERATOR_FN *First_Vertex;     /* Routine that returns the first
                                          vertex in the graph. Used in the
                                          FOR_EACH_VERTEX macro.             */
  GRAPH_ITERATOR_FN *Next_Vertex;      /* Routine that returns the next
                                          vertex in the graph; returns NULL 
                                          if all vertices have been processed.
                                          Used in the FOR_EACH_VERTEX macro. */
  /* PRIVATE */

  GRAPH_DATA *Graph_Data;              /* Data that actually describes the 
                                          graph.                             */
};

#endif
