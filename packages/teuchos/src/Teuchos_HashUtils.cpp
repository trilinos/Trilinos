#include "Teuchos_HashUtils.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;


const int HashUtils::primeCount_ = 29;
const int HashUtils::primes_[] 
= {101, 163, 271, 443, 733, 1187, 1907, 3061, 
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
