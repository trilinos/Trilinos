#ifndef TEUCHOS_RANDOMNUMBERGENERATORBASE_H
#define TEUCHOS_RANDOMNUMBERGENERATORBASE_H

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos
{
  /** 
   * Base class for random-number generators. 
   */
  class RandomNumberGeneratorBase
    {
    public:
      /** Empty ctor */
      RandomNumberGeneratorBase(){;}

      /** Virtual dtor */
      virtual ~RandomNumberGeneratorBase(){;}

      /** generate a vector of random numbers */
      virtual void generateRandomNumbers(int n, double* x) const = 0 ;

    protected:
      static unsigned int microsecondSeed() ;
    private:
    };

}
#endif
