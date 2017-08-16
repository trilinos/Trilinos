#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Kokkos_Core.hpp"

int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc,argv);
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
