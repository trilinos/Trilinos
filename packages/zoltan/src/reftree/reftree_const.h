/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __REFTREE_CONST_H
#define __REFTREE_CONST_H

/* Some constants */

/* Maximum number of vertices per element */
/* Used for dimensioning space for a query function to return vertices */
#define MAXVERT 8

/* Default dimension of the hash table */
#define DEFAULT_HASH_TABLE_SIZE 16384

/* Data structures for refinement tree */

/* The main refinement tree structure */

struct LB_Reftree_Struct {
   LB_ID_PTR global_id;  /* global ID of the corresponding element */
   LB_ID_PTR local_id;   /* local ID of the corresponding element */
   struct LB_Reftree_Struct *children; /* array of the children in the tree */
   int num_child;        /* number of children */
   float *weight;        /* weight of the node; dimension Obj_Weight_Dim */
   float *summed_weight; /* sum of the weights in the subtree rooted at
                            this node */
   float *my_sum_weight; /* sum of weights of nodes assigned to this proc */
   int num_vertex;       /* the number of vertices in the corresponding
                            element */
   int *vertices;        /* the vertices of the corresponding element;
                            local to this processor */
   int in_vertex;        /* starting vertex for determining the path through
                            the children */
   int out_vertex;       /* ending vertex for determining the path through
                            the children */
   int assigned_to_me;   /* for a leaf, 1 if this element is assigned to
                            this processor, 0 if not.  for nonleaves, 1 if
                            the entire subtree is assigned to this proc,
                            0 if none of the subtree, -1 if part */
   int partition;        /* partition to which this node is assigned;
                            meaningful only during the partition algorithm */
};

typedef struct LB_Reftree_Struct LB_REFTREE;

/* Hash table structures */

struct LB_reftree_hash_node {
  LB_ID_PTR gid;            /* Global id */
  LB_REFTREE *reftree_node; /* pointer to a node of the refinement tree */
  struct LB_reftree_hash_node *next;
};

/* data structure pointed to by lb->Data_Structure */

struct LB_reftree_data_struct {
  LB_REFTREE *reftree_root;
  struct LB_reftree_hash_node **hash_table;
  int hash_table_size;
};

/* Prototypes */

extern int LB_Set_Reftree_Param(char *name, char *val);
extern int LB_Reftree_Init(LB *lb);
extern int LB_Reftree_Build(LB *lb);
extern void LB_Reftree_Free_Structure(LB *lb);
extern void LB_Reftree_Print(LB *lb,LB_REFTREE *subroot, int level);

extern LB_REFTREE* LB_Reftree_hash_lookup(LB *lb, 
                                          struct LB_reftree_hash_node **hashtab,
                                          LB_ID_PTR key, int n);
extern void LB_Reftree_Hash_Insert(LB *lb, LB_REFTREE *reftree_node,
                            struct LB_reftree_hash_node **hashtab, int size);
extern void LB_Reftree_Hash_Remove(LB *lb, LB_REFTREE *reftree_node,
                            struct LB_reftree_hash_node **hashtab, int size);
extern void LB_Reftree_Clear_Hash_Table(struct LB_reftree_hash_node **hashtab,
                                 int size);

extern void LB_Get_Child_Order(LB *lb, int *order, int *ierr); /* TEMP child_order */

#endif /* __REFTREE_CONST_H */
