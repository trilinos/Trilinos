#include "Teuchos_RandomNumberGeneratorBase.hpp"

#include <sys/time.h>
#include <unistd.h>

using namespace Teuchos;

unsigned int RandomNumberGeneratorBase::microsecondSeed()
{

	struct timeval tp;
  struct timezone tzp;
  
  gettimeofday(&tp, &tzp);

  unsigned int rtn = tp.tv_usec;
  
  return rtn;
}

