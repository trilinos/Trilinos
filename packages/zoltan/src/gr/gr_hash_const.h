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


#define NUM_HASH_BUCKETS 128

/*
 * PROTOTYPES 
 */

extern void LB_initialize_hash_graph(GRAPH *);

/*
 *  Hash table data structures.
 *  This implementation uses a LIST for each bucket.
 */

typedef int HASH_FN(ID *);

typedef struct List_Entry_Struct {
  struct List_Entry_Struct *Next;
  struct List_Entry_Struct *Prev;
  VERTEX *Vertex;
} LIST_ENTRY;

typedef struct Hash_Table_Struct { 
  int Num_Buckets;
  LIST_ENTRY *Bucket[NUM_HASH_BUCKETS];
  HASH_FN *Hash_Fn;
} HASH_TABLE;

typedef struct Hash_Table_Loop_Struct {
  int Bucket;
  LIST_ENTRY *Entry;
} HASH_TABLE_LOOP;

#endif
