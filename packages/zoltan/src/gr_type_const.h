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

#ifndef __GR_TYPE_CONST_H
#define __GR_TYPE_CONST_H

#ifndef lint
static char *cvs_gr_type_const_h = "$Id$";
#endif

enum Graph_Data_Structures {
  GRAPH_BINARY_TREE,
  GRAPH_HASH_TABLE
};

extern enum Graph_Data_Structures BL_Graph_Type;

#endif
