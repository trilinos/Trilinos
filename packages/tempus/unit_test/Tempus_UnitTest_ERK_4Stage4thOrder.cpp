// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Tempus_UnitTest_Utils.hpp"

#include "../TestModels/SinCosModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_4Stage4thOrder, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperERK_4Stage4thOrder<double>());
  testExplicitRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 4);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_4Stage4thOrder, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Explicit 4 Stage", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_4Stage4thOrder, AppAction)
{
  auto stepper = rcp(new Tempus::StepperERK_4Stage4thOrder<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
