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

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_3Stage3rdOrderTVD, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>());
  testExplicitRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 3);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_3Stage3rdOrderTVD, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Explicit 3 Stage 3rd order TVD", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_3Stage3rdOrderTVD, AppAction)
{
  auto stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
