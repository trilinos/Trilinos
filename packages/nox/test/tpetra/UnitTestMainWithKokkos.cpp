/*! \file UnitTestMainWithKokkos.cpp

\brief Copy of Teuchos_StandardUnitTestMain.cpp with Kokkos initialization/finalization.

*/

#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Kokkos_Core.hpp"

int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc, argv);
  int success = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  Kokkos::finalize();
  return success;
}
