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

#ifndef __VX_CONST_H
#define __VX_CONST_H


#ifndef lint
static char *cvs_vx_const_h = "$Id$";
#endif

#include "id_const.h"

/*
 * KDD -- for now, include both geometric and graph information in vertex.
 * KDD -- later, may prefer to make included data depend on choice of 
 * KDD -- method.
 */

struct Vertex_Struct {
  /* PUBLIC */
  ID Id;                              /* Unique global identifier for vertex */
  int Weight;                         /* Vertex weight for this vertex       */

  /* PRIVATE */
  int Num_Nbors;                      /* Degree of this vertex               */
  ID *Nbor_List;                      /* Array of neighboring vertices' IDs  */
  float *Edge_Weights;                /* Array of edge weights (optional)    */
  double *Coor;                       /* Geometry information for vertex.    */
};

/*
 *  Structure for public access of vertex information.  The implementation
 *  of the functions will depend on the actual implementation of a vertex.
 *  However, since all vertices are stored the same way, we do not need
 *  to keep a copy of the vertex function pointers for each vertex.
 */

#ifdef KDD
/* NOT READY FOR THIS YET! */
struct Vertex_Utils {
  /* PUBLIC */
  VERTEX_EDGE_FN *Get_Edge_List;      /* Function that returns the edge
                                         list for a vertex.  This function
                                         either accesses the storage
                                         in the vertex or calls a user 
                                         function that provides the list.    */
  VERTEX_GEOM_FN *Get_Geom;           /* Function that returns the geometry
                                         info for a vertex.  This function
                                         either accesses the storage in the
                                         vertex or calls a user function that
                                         provides the geometry.              */
};
#endif

#endif
