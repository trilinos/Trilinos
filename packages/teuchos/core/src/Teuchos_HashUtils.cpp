// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_HashUtils.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CompilerCodeTweakMacros.hpp"

using namespace Teuchos;


const int HashUtils::primeCount_ = 33;
const int HashUtils::primes_[]
= {11, 19, 37, 59, 101, 163, 271, 443, 733, 1187, 1907, 3061,
	 4919, 7759, 12379, 19543, 30841, 48487, 75989,
	 119089, 185971, 290347, 452027, 703657, 1093237,
	 1695781, 2627993, 4067599, 6290467, 9718019,
	 15000607, 23133937, 35650091};


int HashUtils::nextPrime(int newCapacity)
{
	TEUCHOS_TEST_FOR_EXCEPTION(newCapacity > primes_[primeCount_-1],
                     std::logic_error,
                     "HashUtils::nextPrime() overflow");

	for (int i=0; i<primeCount_; i++)
		{
			if (newCapacity <= primes_[i])
				{
					return primes_[i];
				}
		}

  TEUCHOS_TEST_FOR_EXCEPTION(true,
                     std::logic_error,
                     "unexpected case in HashUtils::nextPrime()");
	TEUCHOS_UNREACHABLE_RETURN(0);
}

/** helper to hash values larger than int to an int.
 * This is similar to the version in Zoltan2, not a true hash,
 * but this is an improvement over what was done before which was
 * typecasting everything to ints and returning -ve hash codes which was
 * in turn used to index arrays.
 */
int HashUtils::getHashCode(const unsigned char *a, size_t len)
{
  int total=0;
  unsigned char *to = reinterpret_cast<unsigned char *>(&total);
  int c=0;
  for (size_t i=0; i < len; i++){
    to[c++] += a[i];
    if (c == sizeof(int))
      c = 0;
  }
  if (total < 0)
  {
      /* Convert the largest -ve int to zero and -1 to
       * std::numeric_limits<int>::max()
       * */
      size_t maxIntBeforeWrap = std::numeric_limits<int>::max();
      maxIntBeforeWrap ++;
      total += maxIntBeforeWrap;
  }
  return total;
}
