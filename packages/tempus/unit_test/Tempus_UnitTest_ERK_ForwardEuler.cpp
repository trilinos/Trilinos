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
TEUCHOS_UNIT_TEST(ERK_ForwardEuler, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>());
  testExplicitRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
  const auto rk_fe = stepper->getTableau();

  TEUCHOS_ASSERT( rk_fe->isTVD() );
  TEUCHOS_ASSERT( rk_fe->getTVDCoeff() == 1 );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_ForwardEuler, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Forward Euler", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_ForwardEuler, AppAction)
{
  auto stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
