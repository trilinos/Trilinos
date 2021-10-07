// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"


namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;
using Teuchos::sublist;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_2Stage3rdOrder, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperEDIRK_2Stage3rdOrder<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 3);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_2Stage3rdOrder, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("EDIRK 2 Stage 3rd order", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_2Stage3rdOrder, AppAction)
{
  auto stepper = rcp(new Tempus::StepperEDIRK_2Stage3rdOrder<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
