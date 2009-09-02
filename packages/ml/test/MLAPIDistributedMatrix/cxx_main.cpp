
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
#include "MLAPI_DistributedMatrix.h"
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
    int NumGlobalElements = 10;
    Space S(NumGlobalElements);

    Operator B;
    MultiVector V(S), W;

    // ============================================================== //
    // Define the matrix inside, then reassign the DistributedMatrix  //
    // to another Operator, let the DistributedMatrix disappear, then //
    // work with the other Operator.                                  //
    // ============================================================== //
    
    if (true) 
    {
      DistributedMatrix A(S, S);

      // Set the elements, all on processor 0
      if (GetMyPID() == 0)
      {
        for (int i = 0 ; i < S.GetNumGlobalElements() ; ++i)
        {
          A(i, i) = 0.0;
        }
      }
      A.FillComplete();

      // Replace the elements, only local at this point
      for (int i = 0 ; i < S.GetNumMyElements() ; ++i)
      {
        int GlobalRow = S(i);
        A.ReplaceElement(GlobalRow, GlobalRow, 1.0 + GlobalRow);
      }

      B = A + A;
    }

    B = B / 2;

    V = 1;
    W = B * V;

    double ActualNorm = W.NormOne();
    double ExpectedNorm = NumGlobalElements * (NumGlobalElements + 1) / 2;

    Check(ActualNorm, ExpectedNorm);
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
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#endif /* #if defined(HAVE_ML_MLAPI) */
