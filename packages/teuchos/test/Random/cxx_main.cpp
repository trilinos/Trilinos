#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DLARAN.hpp"
#include "Teuchos_SystemRand.hpp"
#include "Teuchos_Out.hpp"
#include "Teuchos_Array.hpp"

using namespace Teuchos;
using std::string;

/* Test of Teuchos random number classes */

int main(int argc, char** argv)
{
  try
    {
      int n = 100000;
      Array<double> x(n);

      RefCountPtr<RandomNumberGeneratorBase> dlaran = rcp(new DLARAN());
      RefCountPtr<RandomNumberGeneratorBase> sysRand = rcp(new SystemRand());

      dlaran->generateRandomNumbers(n, &(x[0]));
      double sum = 0.0;
      for (int i=0; i<n; i++)
        {
          sum += x[i];
        }
      sum = sum/((double) n);

      Out::printf("sum = %g\n", sum);
      
      sysRand->generateRandomNumbers(n, &(x[0]));
      sum = 0.0;
      for (int i=0; i<n; i++)
        {
          sum += x[i];
        }
      sum = sum/((double) n);
      Out::printf("sum = %g\n", sum);

      return 0;
    }
  catch(std::exception& e)
    {
      cerr << e.what() << endl;
    }
}
