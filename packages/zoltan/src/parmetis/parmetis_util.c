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
static char *cvs_parmetis_util = "$Id$";
#endif

#include "lb_const.h"
#include "parmetis_const.h"

/* LB_hashf is a hash function for global ids. 
 *
 * Input:
 *   key, a key to look up of type LB_GID (any data type)
 *   n,   the hash function returns an integer < n
 *
 * Return value:
 *   the hash value, an integer between 0 and n-1
 *
 * Algorithm: This version uses bitwise xor and mod. 
 *            Feel free to replace it with a more sophisticated method.
 *
 * Author: Erik Boman, eboman@cs.sandia.gov (9226)
 */


int LB_hashf(LB_GID key, int n)
{
  unsigned int h, *p;
  int bytes;

  if (sizeof(LB_GID) <= sizeof(int))
    return (((unsigned int) key)%n);

  h = 0;
  for (p = (unsigned int *)&key, bytes=sizeof(LB_GID); bytes >= sizeof(int); 
       bytes-=sizeof(int), p ++){
    h ^= *p;
  }
  return (h%n);
}

/* LB_hash_lookup uses LB_hashf to lookup a key 
 *
 * Input:
 *   hashtab, pointer to the hash table
 *   key, a key to look up of type LB_GID (any data type)
 *   n,   dimension of the hash table
 *
 * Return value:
 *   the global number of the element with key key,
 *   or -1 if the key is not in the hash table
 *
 * Author: Erik Boman, eboman@cs.sandia.gov (9226)
 */

int LB_hash_lookup (struct LB_hash_node **hashtab, LB_GID key, int n)
{
  int i;
  struct LB_hash_node *ptr;

  i = LB_hashf(key, n);
  for (ptr=hashtab[i]; ptr != NULL; ptr = ptr->next){
    if (LB_EQ_GID(ptr->gid, key))
      return (ptr->gno);
  }
  /* Key not in hash table */
  return -1;
}
