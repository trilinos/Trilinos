/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#include <stdio.h>
#include <stdlib.h>
#include "DD.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*  NOTE: See file, README, for associated documentation. (RTH) */



/*****************************************************************************/
/* DD_Hash2 is a hash function for Zoltan ids (local or global).
 * It is derived from Zoltan_Hash.
 *
 * Input:
 *   key: a key to hash of type ZOLTAN_ID_PTR
 *   num_id_entries: the number of (ZOLTAN_ID_TYPE-sized) entries of
 *                    the key to use
 *   n: the range of the hash function is 0..n-1
 *
 * Return value:
 *   the hash value, an unsigned integer between 0 and n-1
 */

#define ZZ_MURMUR_HASH
#ifdef ZZ_MURMUR_HASH

#include "murmur3.c"

unsigned int Zoltan_DD_Hash2(ZOLTAN_ID_PTR key, int num_id_entries,
                             unsigned int n,
                             void *hashdata, ZOLTAN_HASH_FN *fn)
{
/*
 * Algorithm:
 *   Murmurhash3 as in Zoltan_Hash, with a different seed value.
 *   MurmurHash3 was written by Austin Appleby, and is placed in the
 *   public domain. The author hereby disclaims copyright to this source
 *   code.
 */

  uint32_t k;
  MurmurHash3_x86_32((void *)key, sizeof(ZOLTAN_ID_TYPE)*num_id_entries,
                     14, (void *)&k);
  return(k % n);
}

#endif  /* ZZ_MURMUR_HASH */

#ifdef ZZ_KNUTH_HASH

#define ZOLTAN_DD_HASH_CONSTANT 2654435761U   /* consider 516595003U */

unsigned int Zoltan_DD_Hash2(ZOLTAN_ID_PTR key, int num_id_entries,
                             unsigned int n,
                             void *hashdata, ZOLTAN_HASH_FN *fn)
{
/*
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
    
 *   Replaced explict constant with #define below.  Consider changing
 *   to new constant. Changed name to DD_Hash2.  RTH
 */
  unsigned int h, rest, *p, bytes, num_bytes;
  char *byteptr;

  num_bytes = (unsigned int) num_id_entries * sizeof(ZOLTAN_ID_TYPE);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)key, bytes=num_bytes;
       bytes >= (unsigned int) sizeof(int);
       bytes-=sizeof(int), p++){
    h = (h*ZOLTAN_DD_HASH_CONSTANT) ^ (*p);
  }

  /* Then take care of the remaining bytes, if any */
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
    rest = (rest<<8) | (*byteptr);
  }

  /* Merge the two parts */
  if (rest)
    h = (h*ZOLTAN_DD_HASH_CONSTANT) ^ rest;

  /* Return h mod n */
  return (h%n);
}

#endif /* ZOLTAN_KNUTH_HASH */


void Zoltan_DD_default_cleanup (void * hashdata)
{
 ZOLTAN_FREE(&hashdata);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
