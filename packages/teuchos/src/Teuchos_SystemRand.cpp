#include "Teuchos_SystemRand.hpp"


using namespace Teuchos;

SystemRand::SystemRand()
{
	srand(RandomNumberGeneratorBase::microsecondSeed());
}

SystemRand::SystemRand(unsigned int seed)
{
	srand(seed);
}

void SystemRand::generateRandomNumbers(int n, double* x) const
{
	for (int i=0; i<n; i++)
		{
			x[i] = ((double) rand()) / ((double) RAND_MAX);
		}
}



