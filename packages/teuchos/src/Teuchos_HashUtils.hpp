#ifndef HASHUTILS_H
#define HASHUTILS_H

#include "Teuchos_ConfigDefs.hpp"
#include <string>

namespace Teuchos
{
  using std::string;

  /**
   * \ingroup Containers
   * Utilities for generating hashcodes.
   */

  class HashUtils
    {
    public:
      /* get the next prime in a sequence of hashtable sizes */
      static int nextPrime(int newCapacity);

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

  /** \relates HashUtils standard interface for getting the hash code of
   * an object */

  template <class T> int hashCode(const T& x);

  /** \relates HashUtils get the hash code of an int */
  template <> inline int hashCode(const int& x)
    {
      return x;
    }

  /** \relates HashUtils  get the hash code of a double */
  template <> inline int hashCode(const double& x)
    {
      return (int) x;
    }

  /**\relates HashUtils  get the hash code of a bool */
  template <> inline int hashCode(const bool& x)
    {
      return (int) x;
    }


  /** \relates HashUtils get the hash code od a string */
  template <> inline int hashCode(const string& x)
    {
      const char* str = x.c_str();
      int len = x.length();
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

      return rtn;
    }



}
#endif
