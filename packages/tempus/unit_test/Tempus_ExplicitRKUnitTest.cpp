// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_UtilsUnitTest.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;


std::vector<std::string> getRKMethods()
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("Bogacki-Shampine 3(2) Pair");
  RKMethods.push_back("Merson 4(5) Pair");
  RKMethods.push_back("General ERK");
  RKMethods.push_back("RK Forward Euler");
  RKMethods.push_back("RK Explicit 4 Stage");
  RKMethods.push_back("RK Explicit 3/8 Rule");
  RKMethods.push_back("RK Explicit 4 Stage 3rd order by Runge");
  RKMethods.push_back("RK Explicit 5 Stage 3rd order by Kinnmark and Gray");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order TVD");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order by Heun");
  RKMethods.push_back("RK Explicit 2 Stage 2nd order by Runge");
  RKMethods.push_back("RK Explicit Trapezoidal");

  return RKMethods;
}

// Comment out any of the following tests to exclude from build/run.
#define INITIALIZE


#ifdef INITIALIZE
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ExplicitRKUnitTest, initialize)
{
  RCP<Tempus::StepperFactory<double> > sf =
    Teuchos::rcp(new Tempus::StepperFactory<double>());
  auto RKMethods = getRKMethods();

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    // Default construction.
    std::string stepperType = RKMethods[m];
    auto stepper = rcp_dynamic_cast<Tempus::StepperExplicitRK<double> >(
      sf->createStepper(stepperType));
    TEST_ASSERT(!(stepper->isInitialized()));

    auto model = rcp(new Tempus_Test::SinCosModel<double>());
    auto obs   = rcp(new Tempus::StepperExplicitRKObserver<double>());

    StepperInitializeBasic(model, stepper, obs, out, success);

    // Tableau
    stepper->setTableau(stepperType);  TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

  }
}
#endif // INITIALIZE

} // namespace Tempus_Test
