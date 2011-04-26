/*! \file main.cpp

\brief MueLu unit testing main program.

This file is the main for the unit test executable.

NOTE: This file should *not* be built and included as part of the MueLu
library.  It is instead to be directly included in the build files for
specific unit test suites.

*/

#include <Teuchos_UnitTestRepository.hpp>
#include "MueLu_TestHelpers.hpp"

int main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Note: the command line parameter --linAlgebra= is take into accounts. 
  // Cthulhu parameters are added to the Teuchos::CommandLineProcessor of Teuchos::UnitTestRepository in MueLu_TestHelpers.cpp

  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
