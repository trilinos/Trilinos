/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* LB_Hash is a hash function for Zoltan ids (local or global). 
 *
 * Input:
 *   key: a key to hash of type LB_ID_PTR
 *   num_id_entries: the number of (LB_ID_TYPE-sized) entries of the key to use
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


int sfc_hash(LB_ID_PTR key, int num_id_entries, unsigned int n)
{
  unsigned int h, rest, *p;
  char *byteptr;
  int bytes;
  int num_bytes = num_id_entries * sizeof(LB_ID_TYPE);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)key, bytes=num_bytes;
       bytes >= sizeof(int); 
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
