
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
#include "MLAPI_Gallery.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelSA.h"

using namespace MLAPI;

void check(const double expected, const double actual)
{
  if (GetMyPID() == 0)
    cout << "expected norm = " << expected << ", actual = " << actual;
  if (expected != actual)
  {
    cout << ", FAILED!" << endl;
    exit(EXIT_FAILURE);
  }

  if (GetMyPID() == 0)
    cout << ",  OK!" << endl;
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
    Operator I = Gallery("Diag", S);

    MultiVector one(S);   one   = 1.0;
    MultiVector two(S);   two   = 2.0;
    MultiVector three(S); three = 3.0;
    MultiVector four(S);  four  = 4.0;
    MultiVector five(S);  five  = 5.0;
    MultiVector res;

    if (GetMyPID() == 0)
      cout << "Testing `1 + 1 + 1 + 1'... ";
    res = one + one + one + one;
    check(4.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + 2'... ";
    res = one + two;
    check(3.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + I * 2'... ";
    res = one + I * two;
    check(3.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + I * (2 + 2)'... ";
    res = one + I * (two + two);
    check(5.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + 3.5 * I * 2'... ";
    res = one + 3.5 * I * two;
    check(8.0, res.NormInf());

#if 0
    if (GetMyPID() == 0)
      cout << "Testing `1 + I * 2 * 3.5'... ";
    res = one + I * two * 3.5;
    check(8.0, res.NormInf());
#endif

    if (GetMyPID() == 0)
      cout << "Testing `1 + 3.5 * I * 2'... ";
    res = one + 3.5 * I * two;
    check(8.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `5.0 * I * 2 + 3'... ";
    res = 5.0 * I * two + three;
    check(13.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + 3.5 * I * 2 + 7.0 * I * (I * (2 + 3))'... ";
    res = one + 3.5 * I * two + 7.0 * I * (I * (two + three));
    check(43.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `5 - 3'... ";
    res = five - three;
    check(2.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `4.0 * 1 + 7.0 * 2'... ";
    res = 4.0 * one + 7.0 * two;
    check(18.0, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + 2 + 3'...";
    res = one + two + three;
    check(6, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + 2.0 * 2 + 3.0 * 3'...";
    res = one + 2.0 * two + 3.0 * three;
    check(14, res.NormInf());

    if (GetMyPID() == 0)
      cout << "Testing `1 + 2.0 * 2 + 3.0 * 3'...";
    res = one + 2.0 * two + 3.0 * three;
    check(14, res.NormInf());

    if (GetMyPID() == 0)
      cout << endl << "All test passed!" << endl;
  }
  catch (const int e) 
  {
    cout << "Caught integer exception, code = " << e << endl;
  }
  catch (...) 
  {
    cout << "problems here..." << endl;
    exit(EXIT_FAILURE);
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
