// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_HashUtils.hpp"
#include "Teuchos_TestForException.hpp"

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
	TEST_FOR_EXCEPTION(newCapacity > primes_[primeCount_-1],
                     length_error,
                     "HashUtils::nextPrime() overflow");

	for (int i=0; i<primeCount_; i++)
		{
			if (newCapacity <= primes_[i])
				{
					return primes_[i];
				}
		}

  TEST_FOR_EXCEPTION(true,
                     logic_error,
                     "unexpected case in HashUtils::nextPrime()");
	return 0;
}
