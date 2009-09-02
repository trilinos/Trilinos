
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_common.h"

#if defined(HAVE_ML_MLAPI)

#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_SerialMatrix.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelSA.h"

using namespace MLAPI;

void Check(const double ActualNorm, const double ExpectedNorm)
{
  if (GetMyPID() == 0)
  {
    cout << "NormOne, actual = " << ActualNorm;
    cout << ", expected = " << ExpectedNorm;

    if (ActualNorm != ExpectedNorm)
      cout << ",  FAILED!" << endl;
    else
      cout << ",  OK!" << endl;
  }

  if (ActualNorm != ExpectedNorm)
    exit(EXIT_FAILURE);
}

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif
      
  try
  {
    Init();

    if (GetNumProcs() != 1)
      throw("Only with one processor");

    int NumGlobalElements = 10;
    Space S(NumGlobalElements);

    SerialMatrix B; // will hold a copy of current operator
    Operator C; /// will hold a composition of current operator
    MultiVector V(S), W;

    // ============================================================== //
    // Define the matrix inside, then reassign the DistributedMatrix  //
    // to another Operator, let the DistributedMatrix disappear, then //
    // work with the other Operator.                                  //
    // ============================================================== //
    
    if (true)
    {
      SerialMatrix A(S, S);
      for (int i = 0 ; i < S.GetNumGlobalElements() ; ++i)
      {
        A(i, i) = 1.0;
      }

      V = 1.0;
      W = A * V;

      Check(W.NormOne(), NumGlobalElements);

      // perform an operator using A, then change its values. 
      // NOTE: I can only assign B to a SerialMatrix object, not to a generic
      // Operator. Therefore B = A * A will not work, since A * A returns an
      // Operator, not a SerialMatrix.
      B = A;

      for (int i = 0 ; i < S.GetNumGlobalElements() ; ++i)
      {
        A(i, i) = 2.0;
      }

      B = A;
      // A is destroyed at this point.
    }

    // change the values of B
    for (int i = 0 ; i < S.GetNumGlobalElements() ; ++i)
    {
      B(i, i) = 2.0 * B(i, i);
    }

    // reassign B to C, then work with C
    C = B;
    C = C * B;

    // check the norm
    V = 1.0;
    W = C * V;

    Check(W.NormOne(), 16.0 * NumGlobalElements);
  }
  catch (const int e) 
  {
    cout << "Caught integer exception, code = " << e << endl;
  }
  catch (...) 
  {
    cout << "problems here..." << endl;
  } 

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-ifpack --enable-aztecoo");
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#endif /* #if defined(HAVE_ML_MLAPI) */
