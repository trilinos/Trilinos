#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "Teuchos_XMLObject.hpp"

using namespace Teuchos;

int main(int argc, char* argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  /* create an XML object */
  XMLObject solver("solver");
  XMLObject prec("preconditioner");

  solver.addAttribute("krylov method", "gmres");
  solver.addInt("iterations", 1000);
  solver.addDouble("tolerance", 1.0e-10);

  solver.addChild(prec);

  prec.addAttribute("type", "smoothed aggregation");
  prec.addInt("max levels", 4);

  string str = solver.toString();
  cout << str << endl;
        
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
