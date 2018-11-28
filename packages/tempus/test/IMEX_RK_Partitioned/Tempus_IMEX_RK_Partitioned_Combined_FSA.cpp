// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_IMEX_RK_Partitioned_FSA.hpp"

#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

std::string method_name;

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(IMEX_RK_Partitioned, VanDerPol_Combined_FSA)
{
  test_vdp_fsa(method_name, true, false, out, success);
}

TEUCHOS_UNIT_TEST(IMEX_RK_Partitioned, VanDerPol_Combined_FSA_Tangent)
{
  test_vdp_fsa(method_name, true, true, out, success);
}

} // namespace Tempus_Test

int main( int argc, char* argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Add "--method" command line argument
  Teuchos::CommandLineProcessor& CLP = Teuchos::UnitTestRepository::getCLP();
  method_name = "";
  CLP.setOption("method", &method_name, "Stepper method");

  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
