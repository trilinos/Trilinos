// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "zz_util_const.h"
#include "reftree.h"

/* Zoltan_Reftree_hash_lookup uses Zoltan_Hash to lookup a key 
 *
 * Input:
 *   zz, a Zoltan structure
 *   hashtab, pointer to the hash table
 *   key, a key to look up of type GID (any data type)
 *   n,   dimension of the hash table
 *
 * Return value:
 *   a pointer to the refinement tree node with GID key
 *   or NULL if the key is not in the hash table
 *
 * Author: Erik Boman, eboman@cs.sandia.gov (9226) (for parmetis/jostle)
 * Modified for refinement tree nodes by william.mitchell@nist.gov
 */

ZOLTAN_REFTREE* Zoltan_Reftree_hash_lookup (ZZ *zz, 
                                    struct Zoltan_Reftree_hash_node **hashtab,
                                    ZOLTAN_ID_PTR key, int n)
{
  int i;
  struct Zoltan_Reftree_hash_node *ptr;

  i = Zoltan_Hash(key, zz->Num_GID, (unsigned int)n);
  for (ptr=hashtab[i]; ptr != NULL; ptr = ptr->next){
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->reftree_node);
  }
  /* Key not in hash table */
  return (ZOLTAN_REFTREE *)NULL;
}

/* second version for int instead of refinement tree node */

int Zoltan_Reftree_inthash_lookup (ZZ *zz, 
                                   struct Zoltan_Reftree_inthash_node **hashtab,
                                   ZOLTAN_ID_PTR key, int n)
{
  int i;
  struct Zoltan_Reftree_inthash_node *ptr;

  i = Zoltan_Hash(key, zz->Num_GID, (unsigned int)n);
  for (ptr=hashtab[i]; ptr != NULL; ptr = ptr->next){
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->lid);
  }
  /* Key not in hash table */
  return -1;
}

/* Zoltan_Reftree_Hash_Insert adds an entry to the hash table
 *
 * Input:
 *   zz, a Zoltan structure
 *   reftree_node, pointer to a node of the refinement tree
 *   hashtab, pointer to the hash table
 *   size, dimension of the hash table
 *
 * Author: William Mitchell, william.mitchell@nist.gov
 */

void Zoltan_Reftree_Hash_Insert(ZZ *zz, ZOLTAN_REFTREE *reftree_node,
                            struct Zoltan_Reftree_hash_node **hashtab, int size)
{
int i;
struct Zoltan_Reftree_hash_node *new_entry;

  i = Zoltan_Hash(reftree_node->global_id, zz->Num_GID, (unsigned int)size);

  new_entry = (struct Zoltan_Reftree_hash_node *)
              ZOLTAN_MALLOC(sizeof(struct Zoltan_Reftree_hash_node));
  new_entry->gid = ZOLTAN_MALLOC_GID(zz);
  ZOLTAN_SET_GID(zz, new_entry->gid,reftree_node->global_id);
  new_entry->reftree_node = reftree_node;
  new_entry->next = hashtab[i];
  hashtab[i] = new_entry;
}

/* second version for int instead of refinement tree node */

void Zoltan_Reftree_IntHash_Insert(ZZ *zz, ZOLTAN_ID_PTR gid, int lid,
                                   struct Zoltan_Reftree_inthash_node **hashtab,
                                   int size)
{
int i;
struct Zoltan_Reftree_inthash_node *new_entry;

  i = Zoltan_Hash(gid, zz->Num_GID, (unsigned int)size);

  new_entry = (struct Zoltan_Reftree_inthash_node *)
              ZOLTAN_MALLOC(sizeof(struct Zoltan_Reftree_inthash_node));
  new_entry->gid = ZOLTAN_MALLOC_GID(zz);
  ZOLTAN_SET_GID(zz, new_entry->gid,gid);
  new_entry->lid = lid;
  new_entry->next = hashtab[i];
  hashtab[i] = new_entry;
}

/* Zoltan_Reftree_Hash_Remove removes a key from the hash table
 *
 * Input:
 *   zz, a Zoltan structure
 *   hashtab, pointer to the hash table
 *   key, a key to look up of type GID (any data type)
 *   n,   dimension of the hash table
 *
 * Author: William Mitchell, william.mitchell@nist.gov
 */

void Zoltan_Reftree_Hash_Remove (ZZ *zz, ZOLTAN_REFTREE *reftree_node,
                             struct Zoltan_Reftree_hash_node **hashtab, int n)
{
  int i;
  struct Zoltan_Reftree_hash_node *ptr, *prev, *next;

  i = Zoltan_Hash(reftree_node->global_id, zz->Num_GID, (unsigned int)n);
  ptr = hashtab[i];
  prev = NULL;
  while (ptr != NULL) {
    if (ZOLTAN_EQ_GID(zz, ptr->gid, reftree_node->global_id)) {
      next = ptr->next;
      ZOLTAN_FREE(&(ptr->gid));
      ZOLTAN_FREE(&ptr);
      if (prev == NULL) {
        hashtab[i] = next;
      } else {
        prev->next = next;
      }
      ptr = NULL;
    } else {
      prev = ptr;
      ptr = ptr->next;
    }
  }
}

/* Zoltan_Reftree_Clear_Hash_Table empties a hash table and frees the
 * memory, except for the memory of the table itself
 *
 * Input:
 *   hashtab, pointer to the hash table
 *   size, dimension of the hash table
 *
 * Author: William Mitchell, william.mitchell@nist.gov
 */

void Zoltan_Reftree_Clear_Hash_Table(struct Zoltan_Reftree_hash_node **hashtab,
                                 int size)
{
  int i;
  struct Zoltan_Reftree_hash_node *ptr, *next;

  for (i=0; i<size; i++) {
    ptr = hashtab[i];
    while (ptr != NULL) {
      next = ptr->next;
      ZOLTAN_FREE(&(ptr->gid));
      ZOLTAN_FREE(&ptr);
      ptr = next;
    }
    hashtab[i] = (struct Zoltan_Reftree_hash_node *)NULL;
  }

}

/* second version for int instead of refinement tree node */

void Zoltan_Reftree_Clear_IntHash_Table(
                         struct Zoltan_Reftree_inthash_node **hashtab, int size)
{
  int i;
  struct Zoltan_Reftree_inthash_node *ptr, *next;

  for (i=0; i<size; i++) {
    ptr = hashtab[i];
    while (ptr != NULL) {
      next = ptr->next;
      ZOLTAN_FREE(&(ptr->gid));
      ZOLTAN_FREE(&ptr);
      ptr = next;
    }
    hashtab[i] = (struct Zoltan_Reftree_inthash_node *)NULL;
  }

}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
