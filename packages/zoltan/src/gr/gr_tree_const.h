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

#ifndef __GR_TREE_CONST_H
#define __GR_TREE_CONST_H

#ifndef lint
static char *cvs_gr_tree_const_h = "$Id$";
#endif

#include "vx_const.h"
/*
 * PROTOTYPES 
 */

extern void initialize_tree_graph(GRAPH *);

/*
 *  Tree data structures
 */

typedef struct Tree_Struct {
  struct Tree_Struct *Left;     /* pointer to left subtree                  */
  struct Tree_Struct *Right;    /* pointer to right subtree                 */
  struct Tree_Struct *Parent;   /* pointer to parent node in tree           */
                                /* This pointer is used for the 
                                   FOR_EACH_VERTEX traversal of the tree.
                                   It should be replaced if a different
                                   traversal (e.g., a doubly-linked list)
                                   is implemented.                          */
  int Bal;                      /* AVL balance factor                       */
  VERTEX *Vertex;               /* Pointer to VERTEX.                       */
} TREE;

/*
 *  Flags controlling the FOR_EACH_VERTEX loops over trees.
 *  This structure could be replaced by something much simpler if
 *  the First_Vertex and Next_Vertex routines are placed by a 
 *  different traversal (e.g., a doubly-linked list).
 */

typedef struct Tree_Loop_Struct { 
  TREE *Current_Node;
  unsigned int Go_Left  :1;
  unsigned int Go_Mid   :1;
  unsigned int Go_Right :1;
  unsigned int Go_Up    :1;
} TREE_LOOP;

#endif
