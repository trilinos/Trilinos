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
#include "parmetis_jostle_const.h"

#ifndef LB_NO_PARMETIS

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
#endif  /* !LB_NO_PARMETIS */
