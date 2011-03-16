// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_HASH_HPP_
#define _ZOLTAN2_HASH_HPP_

/*! \file Zoltan2_Hash.hpp
*/

#include <climits>
#include <Teuchos_HashUtils.hpp>

namespace Z2
{

/* 
 *  TODO - doxygen'ify this information
 * Input:
 *   key: a key to hash of type std::string, int, double or bool.
 *   n: the hash function returns an int in the range 0..n-1
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

// PHI is golden ratio
#define ZOLTAN2_PHI 1.618033887 
static size_t Zoltan2HashFactor=0;

template <typename T>
unsigned int Zoltan_Hash(T &key, unsigned int n)
{
  unsigned int hashValue;

  if (Zoltan2HashFactor == 0)
    Zoltan2HashFactor = INT_MAX / ZOLTAN2_PHI;

  int ikey = Teuchos::hashCode(key);

  size_t h = ikey * Zoltan2HashFactor;

  if (n){
    hashValue = h%n;
  }
  else{
    hashValue = static_cast<unsigned int>(h);
  }

  return h;
}

}                   // namespace Z2

#endif
