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

/* LB_hashf is a hash function for global ids. LB_GID can be any data type. */

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

/* LB_hash_lookup uses LB_hashf to lookup a certain key */

int LB_hash_lookup (struct LB_hash_node **hashtab, int n, LB_GID key)
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
