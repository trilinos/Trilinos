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
  testExplicitRKAccessorsFullConstruction("RK Explicit 4 Stage");
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
  testExplicitRKAppAction("RK Explicit 4 Stage", out, success);
}


} // namespace Tempus_Test
