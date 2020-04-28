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
TEUCHOS_UNIT_TEST(ERK_5Stage3rdOrderKandG, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>());
  testExplicitRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 3);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_5Stage3rdOrderKandG, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Explicit 5 Stage 3rd order by Kinnmark and Gray", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_5Stage3rdOrderKandG, AppAction)
{
  auto stepper = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
