/*! \file main.cpp

\brief MueLu unit testing main program.

This file is the main for the unit test executable.

NOTE: This file should *not* be built and included as part of the MueLu
library.  It is instead to be directly included in the build files for
specific unit test suites.

*/


#include <Teuchos_UnitTestRepository.hpp>
//#include <Teuchos_GlobalMPISession.hpp>

int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  return ret;
}
