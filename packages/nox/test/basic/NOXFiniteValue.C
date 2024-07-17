// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests NOX's finite number test.

#include "NOX.H"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "NOX_StatusTest_FiniteValue.H"
#include "Teuchos_ScalarTraits.hpp"

using namespace std;

int main(int argc, char *argv[])
{

  int status = 0;

  int myPID = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPID);
#endif

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;


  NOX::StatusTest::FiniteValue fv_test;

  double finite = 1.0e-100;
  double nan = sqrt(-1.0);
  double infinity = 1.0 / Teuchos::ScalarTraits<double>::zero();

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  int result = fv_test.finiteNumberTest(finite);

  if (result != 0)
    status = 1;  // Nonzero is failure

  if ( verbose && (myPID == 0) ) {
    if (result == 0)
      std::cout << "\nFinite value = " << finite
       << ", test correctly identified value as finite!" << std::endl;
    else
      std::cout << "Finite value = " << finite
       << ", test failed to identify value as finite!" << std::endl;
  }

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  result = fv_test.finiteNumberTest(nan);

  if (result != -1)
    status = 1;  // Nonzero is failure

  if ( verbose && (myPID == 0) ) {
    if (result == -1)
      std::cout << "NaN value = " << nan
       << ", test correctly identified value as nan!" << std::endl;
    else
      std::cout << "NaN value = " << nan
       << ", test failed to identify value as nan!" << std::endl;
  }

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  result = fv_test.finiteNumberTest(infinity);

  if (result != -2)
    status = 1;  // Nonzero is failure

  if ( verbose && (myPID == 0) ) {
    if (result == -2)
      std::cout << "Inf value = " << infinity
       << ", test correctly identified value as inf!" << std::endl;
    else
      std::cout << "Inf value = " << infinity
       << ", test failed to identify value as inf!" << std::endl;
  }



#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (myPID == 0) {
    if (status == 0)
      std::cout << "\nTest passed!" << std::endl;
    else
      std::cout << "\nTest Failed!" << std::endl;
  }

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
