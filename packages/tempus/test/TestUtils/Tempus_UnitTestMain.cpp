#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tempus_UnitTestMainUtils.hpp"

int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Tempus_Test::enableFPE(true);
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
