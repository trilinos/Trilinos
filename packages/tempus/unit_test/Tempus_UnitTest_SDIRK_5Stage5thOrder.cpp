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
TEUCHOS_UNIT_TEST(SDIRK_5Stage5thOrder, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperSDIRK_5Stage5thOrder<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 5);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_5Stage5thOrder, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("SDIRK 5 Stage 5th order", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_5Stage5thOrder, StageNumber)
{
  auto stepper    = rcp(new Tempus::StepperSDIRK_5Stage5thOrder<double>());
  int stageNumber = 4;
  stepper->setStageNumber(stageNumber);
  int s = stepper->getStageNumber();
  TEST_COMPARE(s, ==, stageNumber);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_5Stage5thOrder, AppAction)
{
  auto stepper = rcp(new Tempus::StepperSDIRK_5Stage5thOrder<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
