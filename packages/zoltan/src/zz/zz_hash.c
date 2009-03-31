/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <ctype.h>

#include "zz_util_const.h"
#include "zoltan_mem.h"
#include "zz_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef ZZ_NEW_HASH
/* EBEB Slightly improved hash function that fixes some weaknesses
        in the old hash function. An unfortunate side effect is
        that the new hash function changes some answers, therefore
        we have not yet replaced the old version. We should consider
        even stronger hash (like MD5) if we think good hashing is
        very important. 
 */

/* Zoltan_Hash is a hash function for Zoltan ids (local or global). 
 *
 * Input:
 *   key: a key to hash of type ZOLTAN_ID_PTR
 *   num_id_entries: the number of (ZOLTAN_ID_TYPE-sized) entries of the key to use
 *   n: the range of the hash function is 0..n-1
 *      (n=0 returns an unsigned int in 0..INT_MAX; a bit faster)
 *
 * Return value:
 *   the hash value, an unsigned integer between 0 and n-1
 *
 * Algorithm: 
 *   This hash function is a variation of Fibonacci hashing,
 *   a multiplicative method. See e.g. Don Knuth's TAOCP,
 *   volume 3, for background information. Bitwise xor is used for 
 *   keys longer than an int. 
 *
 *   This hash function should be replaced with a stronger method,
 *   like MD5, if high-quality hashing is important.
 *
 * Author: 
 *   Erik Boman, eboman@cs.sandia.gov 
 */

/* Let phi be the golden ratio  = 1.618033887 */
#define MAXINT_DIV_PHI  2654435761U   
   /* =2^32/phi,  best for 32 bit machines */
/* #define MAXINT_DIV_PHI  11400714819323198485U    */
   /* =2^64/phi,  best for 64 bit machines */

unsigned int Zoltan_Hash(ZOLTAN_ID_PTR key, int num_id_entries, unsigned int n)
{
  unsigned int h, rest, *p, bytes, num_bytes;
  char *byteptr;

  num_bytes = (unsigned int) num_id_entries * sizeof(ZOLTAN_ID_TYPE);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)key, bytes=num_bytes;
       bytes >= (unsigned int) sizeof(int); 
       bytes-=sizeof(int), p++){
    h = (h^(*p))*MAXINT_DIV_PHI;
  }

  /* Then take care of the remaining bytes, if any */
  /* This will never be executed when ZOLTAN_ID_PTR points to ints */
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
    rest = (rest<<8) | (*byteptr);
  }

  /* If extra bytes, merge the two parts */
  if (rest)
    h = (h^rest)*MAXINT_DIV_PHI;

  /* At this point h is a good hash value but it's not in the range (0,n-1).
   * We would like to return (h*n/INT_MAX), but this is difficult to compute
   * due to overflow. A simple alternative is to return (h%n) but this
   * is not quite satisfactory since the leading bits are more evenly 
   * distributed than the trailing bits. Therefore we first shift
   * the bits around.
   */

  if (n){
    /* Shift targeted to 32-bit ints but works for other cases too. */
    /* Take h mod n to get a value in [0,n) */
    h = ((h<<15)^(h>>17))%n;
  }

  return h;

}

#else
/* EBEB Old hash function that is compatible with the answer files. */

/* Zoltan_Hash is a hash function for Zoltan ids (local or global). 
 *
 * Input:
 *   key: a key to hash of type ZOLTAN_ID_PTR
 *   num_id_entries: the number of (ZOLTAN_ID_TYPE-sized) entries of the key to use
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


unsigned int Zoltan_Hash(ZOLTAN_ID_PTR key, int num_id_entries, unsigned int n)
{
  unsigned int h, rest, *p, bytes, num_bytes;
  char *byteptr;

  num_bytes = (unsigned int) num_id_entries * sizeof(ZOLTAN_ID_TYPE);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)key, bytes=num_bytes;
       bytes >= (unsigned int) sizeof(int); 
       bytes-=sizeof(int), p++){
    h = (h*2654435761U) ^ (*p);
  }

  /* Then take care of the remaining bytes, if any */
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
    rest = (rest<<8) | (*byteptr);
  }

  /* Merge the two parts */
  if (rest)
    h = (h*2654435761U) ^ rest;

  /* Return h mod n */
  return (h%n);
}
#endif /* ZZ_NEW_HASH */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

