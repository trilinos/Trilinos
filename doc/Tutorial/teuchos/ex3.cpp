#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "Teuchos_CommandLineProcessor.hpp"

using namespace Teuchos;

int main(int argc, char* argv[])
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

  try {
    
    CommandLineProcessor CLP;
    
    CLP.recogniseAllOptions(false);
    CLP.throwExceptions(false);
      
    int NumIters = 1550;
    CLP.setOption("iterations",&NumIters,"number of iterations");

    double Tolerance = 1e-10;    
    CLP.setOption("tolerance",&Tolerance,"tolerance");
    
    bool Precondition;
    CLP.setOption("precondition","no-precondition",
                  &Precondition,"prec flag");
    
    CLP.parse(argc,argv);
    cout << NumIters << endl;
    cout << Tolerance << endl;
    cout << Precondition << endl;
    
  } catch(std::exception& e) {
    cerr << "caught exception " << e.what() << endl;
  } 
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
