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

#ifndef __GR_HASH_CONST_H
#define __GR_HASH_CONST_H

#ifndef lint
static char *cvs_gr_hash_const_h = "$Id$";
#endif

/*
 * PROTOTYPES 
 */

extern void initialize_hash_graph(GRAPH *);

/*
 *  Hash table data structures.
 */

typedef struct Hash_Table_Struct {
  int Num_Buckets;
  VERTEX **Buckets;
} HASH_TABLE;

#endif
