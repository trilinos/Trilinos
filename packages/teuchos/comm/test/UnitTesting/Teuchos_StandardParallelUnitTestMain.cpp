// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Teuchos_StandardUnitTestMain.cpp

\brief Standard Parallel Unit testing main program.

This file is ment to be used as a standard main program for a unit test
executable for parallel unit tests.

NOTE: This file should *not* be built and included as part of the Teuchos
library.  It is instead to be directly included in the build files for
specific unit test suites.

*/


#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"


int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
