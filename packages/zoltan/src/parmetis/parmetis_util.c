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
 * Author: Erik Boman, eboman@cs.sandia.gov (9226)
 */


unsigned int LB_hashf(LB_GID key, int n)
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

