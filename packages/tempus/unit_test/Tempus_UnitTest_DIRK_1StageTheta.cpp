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
TEUCHOS_UNIT_TEST(DIRK_1StageTheta, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_1StageTheta<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
  double theta = 0.5;
  TEUCHOS_ASSERT(stepper->getTheta() == theta);
  stepper->setTheta(theta);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1StageTheta, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("DIRK 1 Stage Theta Method", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1StageTheta, AppAction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_1StageTheta<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
