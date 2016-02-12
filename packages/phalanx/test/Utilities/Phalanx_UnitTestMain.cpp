#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Phalanx_config.hpp"
#include "Phalanx_KokkosUtilities.hpp"

int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc,argv);
  {
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setOutputToRootOnly(0);
    PHX::Device::execution_space::print_configuration(out);
  }
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
