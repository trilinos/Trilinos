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

#include "lb_const.h"
#include "lb_util_const.h"
#include "reftree_const.h"
#include "all_allo_const.h"

/* LB_Reftree_hashf is a hash function for global ids. 
 *
 * Input:
 *   key: a key to hash of type LB_GID (any data type)
 *   n: the desired range of the hash function is 0..n-1
 *
 * Return value:
 *   the hash value, an unsigned integer between 0 and n-1
 *
 * Algorithm: 
 *   This hash function uses bitwise xor to hash into an unsigned int,
 *   and then finally employs Don Knuth's golden ratio multiplicative
 *   method. This hash function is good for int-sized keys
 *   (integers and pointers) but may be poor for longer keys.
 *   Feel free to replace it with a more sophisticated method.
 *
 * Author: Erik Boman, eboman@cs.sandia.gov (9226) (for parmetis/jostle)
 * Modified for refinement tree nodes by william.mitchell@nist.gov
 */


unsigned int LB_Reftree_hashf(LB_GID key, int n)
{
  unsigned int h, rest, *p;
  char *byteptr;
  int bytes;

  /* First xor the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)&key, bytes=sizeof(LB_GID); bytes >= sizeof(int); 
       bytes-=sizeof(int), p++){
    h ^= *p;
  }
  /* Then take care of the remaining bytes, if any */
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
    rest = (rest<<8) | (*byteptr);
  }

  /* Compute hash value based on Knuth's golden ratio mult. method */
  h = (h ^ rest) * 2654435761U;

  /* Take the hash value mod n */
  h = h%n;

  return h;
}

/* LB_Reftree_hash_lookup uses LB_hashf to lookup a key 
 *
 * Input:
 *   hashtab, pointer to the hash table
 *   key, a key to look up of type LB_GID (any data type)
 *   n,   dimension of the hash table
 *
 * Return value:
 *   a pointer to the refinement tree node with LB_GID key
 *   or NULL if the key is not in the hash table
 *
 * Author: Erik Boman, eboman@cs.sandia.gov (9226) (for parmetis/jostle)
 * Modified for refinement tree nodes by william.mitchell@nist.gov
 */

LB_REFTREE* LB_Reftree_hash_lookup (struct LB_reftree_hash_node **hashtab,
                                    LB_GID key, int n)
{
  int i;
  struct LB_reftree_hash_node *ptr;

  i = LB_hashf(key, n);
  for (ptr=hashtab[i]; ptr != NULL; ptr = ptr->next){
    if (LB_EQ_GID(ptr->gid, key))
      return (ptr->reftree_node);
  }
  /* Key not in hash table */
  return (LB_REFTREE *)NULL;
}

/* LB_Reftree_Hash_Insert adds an entry to the hash table
 *
 * Input:
 *   reftree_node, pointer to a node of the refinement tree
 *   hashtab, pointer to the hash table
 *   size, dimension of the hash table
 *
 * Author: William Mitchell, william.mitchell@nist.gov
 */

void LB_Reftree_Hash_Insert(LB_REFTREE *reftree_node,
                            struct LB_reftree_hash_node **hashtab, int size)
{
int i;
struct LB_reftree_hash_node *new_entry;

  i = LB_hashf(reftree_node->global_id, size);

  new_entry = (struct LB_reftree_hash_node *)
              LB_MALLOC(sizeof(struct LB_reftree_hash_node));
  LB_SET_GID(new_entry->gid,reftree_node->global_id);
  new_entry->reftree_node = reftree_node;
  new_entry->next = hashtab[i];
  hashtab[i] = new_entry;
}

/* LB_Reftree_Clear_Hash_Table empties a hash table and frees the
 * memory, except for the memory of the table itself
 *
 * Input:
 *   hashtab, pointer to the hash table
 *   size, dimension of the hash table
 *
 * Author: William Mitchell, william.mitchell@nist.gov
 */

void LB_Reftree_Clear_Hash_Table(struct LB_reftree_hash_node **hashtab,
                                 int size)
{
  int i;
  struct LB_reftree_hash_node *ptr, *next;

  for (i=0; i<size; i++) {
    ptr = hashtab[i];
    while (ptr != NULL) {
      next = ptr->next;
      LB_FREE(&ptr);
      ptr = next;
    }
    hashtab[i] = (struct LB_reftree_hash_node *)NULL;
  }

}
