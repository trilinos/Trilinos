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
#include "parmetis_jostle_const.h"

/* LB_hashf is a hash function for global ids. 
 *
 * Input:
 *   key: a key to hash of type LB_GID (any data type)
 *   n: the range of the hash function is 0..n-1
 *
 * Return value:
 *   the hash value, an unsigned integer between 0 and n-1
 *
 * Algorithm: 
 *   This hash function is based on Don Knuth's golden ratio
 *   multiplicative method. Bitwise xor is used for keys
 *   longer than an int. The method works well for keys
 *   of size one or two ints, which is typically the case.
 *
 *   This hash function should be replaced with a stronger method
 *   if good hashing of a large number of keys is important.
 *
 * Author: 
 *   Erik Boman, eboman@cs.sandia.gov (SNL 9226)
 */


unsigned int LB_hashf(LB_GID key, int n)
{
  unsigned int h, rest, *p;
  char *byteptr;
  int bytes;

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)&key, bytes=sizeof(LB_GID); bytes >= sizeof(int); 
       bytes-=sizeof(int), p++){
    h = (h*2654435761U) ^ (*p);
  }

  /* Then take care of the remaining bytes, if any */
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
    rest = (rest<<8) | (*byteptr);
  }

  /* Merge the two parts */
  h = (h*2654435761U) ^ rest;

  /* Return h mod n */
  h = h%n;
  return h;
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

