#include "Teuchos_DLARAN.hpp"


using namespace Teuchos;

DLARAN::DLARAN()
	: RandomNumberGeneratorBase(), seed48_()
{
	init(microsecondSeed());
}

DLARAN::DLARAN(unsigned int seed)
	: RandomNumberGeneratorBase(), seed48_()
{
	init(seed);
}

void DLARAN::init(unsigned int inputSeed)
{
	unsigned int seed = inputSeed;
	int r = 4096;
  int r2 = r*r;
  seed48_[0] = 0;
  seed48_[1] = seed/r2;
  seed48_[2] = (seed - seed48_[1]*r2)/r;
  seed48_[3] = (seed - seed48_[1]*r2 - seed48_[2]*r);

  if ((seed48_[3] % 2) == 0)
    {
      seed48_[3]++;
    }
}

void DLARAN::generateRandomNumbers(int n, double* x) const
{
	for (int i=0; i<n; i++)
		{
			*x++ = getNext();
		}
}

double DLARAN::getNext() const
{
	static const int m1=494, m2=322, m3=2508, m4=2549;
	static const int ipw2=4096;
	static const double r = 1.0/4096.0;

	int it1, it2, it3, it4;
	
	it4 = seed48_[3]*m4;
	it3 = it4/ipw2;
	it4 = it4 - ipw2*it3;
	it3 = it3 + seed48_[2]*m4 + seed48_[3]*m3;
	it2 = it3/ipw2;
	it3 = it3 - ipw2*it2;
	it2 = it2 + seed48_[1]*m4 + seed48_[2]*m3 + seed48_[3]*m2;
	it1 = it2/ipw2;
	it2 = it2 - ipw2*it1;
	it1 = it1 + seed48_[0]*m4 + seed48_[1]*m3 + seed48_[2]*m2 + seed48_[3]*m1;
	it1 = it1 % ipw2;

	const_cast<int&> (seed48_[0]) = it1;
	const_cast<int&> (seed48_[1]) = it2;
	const_cast<int&> (seed48_[2]) = it3;
	const_cast<int&> (seed48_[3]) = it4;

	return r*(it1 + r*(it2 + r*(it3 + r*it4)));
}

