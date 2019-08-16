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
  RKMethods.push_back("RK Backward Euler");
  RKMethods.push_back("IRK 1 Stage Theta Method");
  RKMethods.push_back("Implicit Midpoint");
  RKMethods.push_back("SDIRK 1 Stage 1st order");
  RKMethods.push_back("SDIRK 2 Stage 2nd order");
  RKMethods.push_back("SDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage Theta Method");
  RKMethods.push_back("SDIRK 3 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 5th order");
  RKMethods.push_back("SDIRK 2(1) Pair");

  return RKMethods;
}

// Comment out any of the following tests to exclude from build/run.
#define INITIALIZE


#ifdef INITIALIZE
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRKUnitTest, initialize)
{
  RCP<Tempus::StepperFactory<double> > sf =
    Teuchos::rcp(new Tempus::StepperFactory<double>());
  auto RKMethods = getRKMethods();

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    // Default construction.
    std::string stepperType = RKMethods[m];
    auto stepper = rcp_dynamic_cast<Tempus::StepperDIRK<double> >(
      sf->createStepper(stepperType));
    TEST_ASSERT(!(stepper->isInitialized()));

    auto model = rcp(new Tempus_Test::SinCosModel<double>());
    auto obs   = rcp(new Tempus::StepperDIRKObserver<double>());

    StepperInitializeBasic(model, stepper, obs, out, success);

    // Tableau
    stepper->setTableau(stepperType);  TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

  }
}
#endif // INITIALIZE

} // namespace Tempus_Test
