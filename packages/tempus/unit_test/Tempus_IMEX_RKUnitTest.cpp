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
  RKMethods.push_back("IMEX RK 1st order");
  RKMethods.push_back("IMEX RK SSP2"     );
  RKMethods.push_back("IMEX RK ARS 233"  );
  RKMethods.push_back("General IMEX RK"  );

  return RKMethods;
}

// Comment out any of the following tests to exclude from build/run.
#define INITIALIZE


#ifdef INITIALIZE
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RKUnitTest, initialize)
{
  RCP<Tempus::StepperFactory<double> > sf =
    Teuchos::rcp(new Tempus::StepperFactory<double>());
  auto RKMethods = getRKMethods();

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    // Default construction.
    std::string stepperType = RKMethods[m];
    auto stepper = rcp_dynamic_cast<Tempus::StepperIMEX_RK<double> >(
      sf->createStepper(stepperType));
    TEST_ASSERT(!(stepper->isInitialized()));

    // Setup the explicit and implicit VanDerPol ModelEvaluator
    auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
    auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
    // Setup the IMEX Pair ModelEvaluator
    auto model = rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
                                               explicitModel, implicitModel));
    auto obs   = rcp(new Tempus::StepperIMEX_RKObserver<double>());

    StepperInitializeBasic(model, stepper, obs, out, success);

    // ModelPair
    stepper->setModelPair(model);      TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

    stepper->setModelPair(explicitModel, implicitModel);
                                       TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

    // Tableau
    stepper->setTableaus(Teuchos::null);
                                       TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

    stepper->setExplicitTableau("RK Explicit Trapezoidal", Teuchos::null);
                                       TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

    stepper->setImplicitTableau("SDIRK 2 Stage 3rd order", Teuchos::null);
                                       TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

  }
}
#endif // INITIALIZE

} // namespace Tempus_Test
