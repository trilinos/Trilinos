// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX.H"  // NOX headers
#include "NOX_LAPACK.H" // NOX LAPACK Interface headers
#include "NOX_TestError.H" // common file for testing
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#else
#endif

int main(int argc, char *argv[]) {

  bool success = false;
  bool verbose = false;
  try {
    // Set up the printing utilities
    Teuchos::ParameterList noxParams;
    Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
    printParams.set("Output Precision", 5);
    if (argc > 1) {
      if (argv[1][0]=='-' && argv[1][1]=='v')
         printParams.set("Output Information",
              NOX::Utils::OuterIteration +
              NOX::Utils::OuterIterationStatusTest +
              NOX::Utils::InnerIteration +
              NOX::Utils::Parameters +
              NOX::Utils::Details +
              NOX::Utils::Warning +
              NOX::Utils::TestDetails);
      else
         printParams.set("Output Information", NOX::Utils::Error);
    }
    NOX::Utils printing(printParams);

    if (printing.isPrintType(NOX::Utils::TestDetails)) {
      std::cout << "Starting lapack/NOX_Group/NOX_Group.exe" << std::endl;
    }

    int status = 0;

    success = (status == 0);

    // Begin real testing here!
    if (success)
      std::cout << "Test passed!" << std::endl;
    else
      std::cout << "Test failed!" << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

/*
  end of file main.cc
*/
