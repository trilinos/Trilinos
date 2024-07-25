//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_2Stage3rdOrder, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperSDIRK_2Stage3rdOrder<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 3);
  std::string gammaType = "3rd Order A-stable";
  TEUCHOS_ASSERT(stepper->getGammaType() == gammaType);
  double gamma = 0.7886751345948128;
  TEUCHOS_ASSERT(stepper->getGamma() == gamma);
  stepper->setGammaType(gammaType);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setGamma(gamma);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_2Stage3rdOrder, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("SDIRK 2 Stage 3rd order", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_2Stage3rdOrder, AppAction)
{
  auto stepper = rcp(new Tempus::StepperSDIRK_2Stage3rdOrder<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
