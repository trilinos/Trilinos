/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "DD.h"



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
     
 *   Replaced explict constant with #define below.  Consider changing
 *   to new constant. Changed name to DD_Hash2.  RTH
 */

#define ZOLTAN_DD_HASH_CONSTANT 2654435761U   /* consider 516595003U */

unsigned int Zoltan_DD_Hash2(ZOLTAN_ID_PTR key, int num_id_entries, unsigned int n)
{
  unsigned int h, rest, *p;
  char *byteptr;
  int bytes;
  int num_bytes = num_id_entries * sizeof(ZOLTAN_ID_TYPE);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)key, bytes=num_bytes;
       bytes >= sizeof(int); 
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
