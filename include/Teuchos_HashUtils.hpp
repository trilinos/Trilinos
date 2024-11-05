// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_HASHUTILS_H
#define TEUCHOS_HASHUTILS_H

/*! \file Teuchos_HashUtils.hpp
    \brief Utilities for generating hashcodes
*/

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos
{
  using std::string;

  /**
   * \ingroup Containers
   * \brief Utilities for generating hashcodes.
   * This is not a true hash ! For all ints and types less than ints
   * it returns the i/p typecasted as an int. For every type more than the
   * size of int it is only slightly more smarter where it adds the bits
   * into int size chunks and calls that an hash. Used with a capacity
   * limited array this will lead to one of the simplest hashes.
   * Ideally this should be deprecated and not used at all.
   */

  class TEUCHOSCORE_LIB_DLL_EXPORT HashUtils
    {
    public:
      /* Get the next prime in a sequence of hashtable sizes */
      static int nextPrime(int newCapacity);
      static int getHashCode(const unsigned char *a, size_t len);

    private:

      // sequence of primes generated via mathematica:
      // Table[Prime[Round[1.5^x]], {x, 8, 36}]
      static const int primeCount_;
      static const int primes_[];
      /*={101, 163, 271, 443, 733, 1187, 1907, 3061,
        4919, 7759, 12379, 19543, 30841, 48487, 75989,
        119089, 185971, 290347, 452027, 703657, 1093237,
        1695781, 2627993, 4067599, 6290467, 9718019,
        15000607, 23133937, 35650091};*/
    };


  /** \relates HashUtils
      \brief Standard interface for getting the hash code of an object
  */
  template <class T> int hashCode(const T& x);

  /** \relates HashUtils
      \brief Get the hash code of an int
  */
  template <> inline int hashCode(const int& x)
    {
      return x;
    }

  /** \relates HashUtils
      \brief Get the hash code of an unsigned
  */
  template <> inline int hashCode(const unsigned& x)
    {
      return HashUtils::getHashCode(
        reinterpret_cast<const unsigned char *>(&x), sizeof(unsigned));
    }

  /** \relates HashUtils
      \brief Get the hash code of a double
  */
  template <> inline int hashCode(const double& x)
    {
      return HashUtils::getHashCode(
        reinterpret_cast<const unsigned char *>(&x), sizeof(double));
    }

  /** \relates HashUtils
      \brief Get the hash code of a bool
  */
  template <> inline int hashCode(const bool& x)
    {
      return (int) x;
    }

  /** \relates HashUtils
      \brief Get the hash code of a long long
  */
  template <> inline int hashCode(const long long& x)
    {
      return HashUtils::getHashCode(
        reinterpret_cast<const unsigned char *>(&x), sizeof(long long));
    }

  /** \relates HashUtils
      \brief Get the hash code of a long
  */
  template <> inline int hashCode(const long& x)
    {
      return HashUtils::getHashCode(
        reinterpret_cast<const unsigned char *>(&x), sizeof(long));
    }

  /** \relates HashUtils
      \brief Get the hash code of a std::string
  */
  template <> inline int hashCode(const std::string& x)
    {
      /* This specialization could use the HashUtils::getHashCode as well,
       * but they are both true hashes anyway, so leaving it !
       * */
      const char* str = x.c_str();
      int len = static_cast<int>(x.length());
      int step = len/4 + 1;
      int base = 1;
      int rtn = 0;

      for (int i=0; i<len/2; i+=step)
        {
          rtn += base*(int) str[i];
          base *= 128;
          rtn += base*(int) str[len-i-1];
          base *= 128;
        }

      if (rtn < 0)
      {
          /* Convert the largest -ve int to zero and -1 to
           * std::numeric_limits<int>::max()
           * */
          size_t maxIntBeforeWrap = std::numeric_limits<int>::max();
          maxIntBeforeWrap ++;
          rtn += maxIntBeforeWrap;
      }
      return rtn;
    }

}
#endif // TEUCHOS_HASHUTILS_H
