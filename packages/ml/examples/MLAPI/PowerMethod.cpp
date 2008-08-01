/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_common.h"
#if defined(HAVE_ML_MLAPI) && defined(HAVE_ML_GALERI)

#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Expressions.h"

using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // initialize the MLAPI workspace
    
    Init();

    // Global dimension of our problem. It must be a square number
    
    int NumGlobalRows = 16;

    // Create the space in which vectors and operators live
    
    Space MySpace(NumGlobalRows);

    // Define vectors and operator, based on `MySpace'. The linear system
    // matrix is created using the MLAPI Gallery, which is only a 
    // convenient wrapper of the Galeri package.
    
    MultiVector x(MySpace), y(MySpace);
    Operator A = Gallery("Tridiag", MySpace);

    // We can now start coding the power method. We want a random vector of
    // unitary 2-norm. First, we set random elements in the vector. Then,
    // we divide each element by the 2-norm of the vector, computed as
    // x * x. For two general vectors, x * y represents the scalar product
    // between x and y (defined on the same space).
    
    x.Random();
    x = x / sqrt(x * x);

    // ==================== //
    // loop of power method //
    // ==================== //
    
    int MaxIters = 10;

    for (int i = 0 ; i < MaxIters ; ++i) {

      y = A * x;  // matrix-vector product

      // note that you need parenthesis in the following expression
      // and that x * y is a global operation (therefore all
      // processes must execute it)
      
      double RQ = (y * x) / (x * x);
      if (GetMyPID() == 0)
        cout << "iter = " << i << ", RQ = " << RQ << endl;

      x = y / sqrt(y * y);

    }

    // finalize the MLAPI workspace

    Finalize();

  } 
  catch (const int e) {
    cout << "Caught integer exception, code = " << e << endl;
  } 
  catch (...) {
    cout << "problems here..." << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);

}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("This MLAPI example requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("\t--enable-galeri");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif // #if defined(HAVE_ML_MLAPI)
