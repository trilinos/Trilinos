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
TEUCHOS_UNIT_TEST(DIRK_2Stage2ndOrderLobattoIIIB, Default_Construction)
{
  auto stepper =
      rcp(new Tempus::StepperDIRK_2Stage2ndOrderLobattoIIIB<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_2Stage2ndOrderLobattoIIIB, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Implicit 2 Stage 2nd order Lobatto IIIB", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_2Stage2ndOrderLobattoIIIB, AppAction)
{
  auto stepper =
      rcp(new Tempus::StepperDIRK_2Stage2ndOrderLobattoIIIB<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
