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


#include <stdio.h>
#include <ctype.h>

#include "zz_util_const.h"
#include "zoltan_mem.h"
#include "zz_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

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
 */

#ifndef HAVE_ZOLTAN_KNUTH_HASH
#define ZZ_MURMUR_HASH
#else
#define ZZ_KNUTH_HASH
#endif

#ifdef ZZ_MURMUR_HASH
#include "murmur3.c"

unsigned int Zoltan_Hash(ZOLTAN_ID_PTR key, int num_id_entries, unsigned int n)
{
 /* Using Murmurhash3:
  * MurmurHash3 was written by Austin Appleby, and is placed in the
  * public domain. The author hereby disclaims copyright to this source
  * code.
  */

  uint32_t k;
  MurmurHash3_x86_32((void *)key, sizeof(ZOLTAN_ID_TYPE)*num_id_entries,  
                     1, (void *)&k);
  return(k % n);
}

#endif  /* ZZ_MURMUR_HASH */

#ifdef ZZ_NEW_HASH
/* Let phi be the golden ratio  = 1.618033887 */
#define MAXINT_DIV_PHI  2654435761U   
   /* =2^32/phi,  best for 32 bit machines */
/* #define MAXINT_DIV_PHI  11400714819323198485U    */
   /* =2^64/phi,  best for 64 bit machines */

unsigned int Zoltan_Hash(ZOLTAN_ID_PTR key, int num_id_entries, unsigned int n)
{
/* EBEB Slightly improved hash function that fixes some weaknesses
        in the old hash function. An unfortunate side effect is
        that the new hash function changes some answers, therefore
        we have not yet replaced the old version. We should consider
        even stronger hash (like MD5) if we think good hashing is
        very important. 
 */
/* Algorithm: 
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
  unsigned int h, rest, *p, bytes, num_bytes;
  char *byteptr;
  unsigned int low_order_id[10];

  /* We want the same hash value for 64 bit ZOLTAN_ID_TYPE that we get for 
   * 32 bit ZOLTAN_ID_TYPE so that test answers do not change.  For 64 bit
   * ZOLTAN_ID_TYPE, we still get a good spread of hash values.
   */

  for (h=0; (h < num_id_entries) && (h < 10); h++){
    low_order_id[h] = key[h] & 0xffff;
  }

  num_bytes = (unsigned int) num_id_entries * sizeof(unsigned int);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)low_order_id, bytes=num_bytes;
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

#endif


#ifdef ZZ_KNUTH_HASH

unsigned int Zoltan_Hash(ZOLTAN_ID_PTR key, int num_id_entries, unsigned int n)
{
/* EBEB Old hash function that is compatible with the answer files. */
/* Algorithm: 
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
  unsigned int h, rest, *p, bytes, num_bytes;
  char *byteptr;
  unsigned int low_order_id[10];

  /* We want the same hash value for 64 bit ZOLTAN_ID_TYPE that we get for 
   * 32 bit ZOLTAN_ID_TYPE so that test answers do not change.  For 64 bit
   * ZOLTAN_ID_TYPE, we still get a good spread of hash values.
   */

  for (h=0; (h < num_id_entries) && (h < 10); h++){
    low_order_id[h] = key[h] & 0xffff;
  }

  num_bytes = (unsigned int) num_id_entries * sizeof(unsigned int);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)low_order_id, bytes=num_bytes;
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
#endif /* ZZ_KNUTH_HASH */


/* Zoltan_Recommended_Hash_Size recommends a hash table size that can be used
 * with Zoltan_Hash.
 *
 * Input : #entries that will be stored in the hash table. 
 * 
 * Output : recommended size for the the hash table. This is also the range for 
 * Zoltan_Hash. This is a simple 
 * implementation that tries to bring the collisons to almost zero. The 
 * recommeded hash size is the next largest prime number that is approximately 
 * in the middle of 2^x and 2^(x+1) for some x.
 * */
unsigned int Zoltan_Recommended_Hash_Size (unsigned int n)
{
    /* Prime number approximately in the middle of the range [2^x..2^(x+1)]
     * is in primes[x-1]. Every prime number stored is approximately two times
     * the previous one, so hash table size doubles every time.
     */
    unsigned int primes[] = {
    3, 7, 13, 23, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593,
    49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469,
    12582917, 25165842, 50331653, 100663319, 201326611, 402653189,
    805306457, 1610612741 } ;
    int i, hsize ;

    /* SRSR : err on the side of performance and choose the next largest prime 
     * number. One can also choose primes[i-1] below to cut the memory by half.
     */
    hsize = primes[29] ;
    for (i = 0 ; i < 30 ; i++)
    {
        if (n <= primes[i])
        {
            /*hsize = (i == 0 ? n : primes[i-1]) ;*/
            hsize = primes[i] ;
            break ;
        }
    }

    return hsize ;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

