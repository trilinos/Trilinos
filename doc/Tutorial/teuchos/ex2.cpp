#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Hashtable.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"
#include "complex.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace Teuchos;
using std::string;

int main(int argc, char* argv[])
{
  try {
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

    Teuchos::SerialDenseMatrix<int,float> Matrix(4,4);
    Teuchos::SerialDenseVector<int,float> Vector(4);

    Matrix.random();
    Vector.random();

    int ipiv[4], info;
    
    Teuchos::LAPACK<int,float> L;

    cout << Matrix;
    
    L.GETRF(4,4,Matrix.values(),4,ipiv,&info);
    
    cout << Matrix;
    
    
  } catch(std::exception& e) {
    cerr << "caught exception " << e.what() << endl;
  } 
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
