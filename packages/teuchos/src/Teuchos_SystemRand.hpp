#ifndef TEUCHOS_SYSTEMRAND_H
#define TEUCHOS_SYSTEMRAND_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RandomNumberGeneratorBase.hpp"

namespace Teuchos
{
  /**
   * Use the system's rand() function to fill a random vector.
   *
   * @author Kevin Long
   */
  class SystemRand : public RandomNumberGeneratorBase
    {
    public:
      /** Empty ctor uses system's default seed */
      SystemRand();
      /** Initialize to an integer seed */
      SystemRand(unsigned int seed);

      /** the usual virtual dtor */
      virtual ~SystemRand() {;}

      /** generate a vector of random uniform unit deviates */
      virtual void generateRandomNumbers(int n, double* x) const ;

    private:
    };


}

#endif
