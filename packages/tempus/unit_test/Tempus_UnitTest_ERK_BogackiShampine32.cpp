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
TEUCHOS_UNIT_TEST(ERK_BogackiShampine32, Default_Construction)
{
  testExplicitRKAccessorsFullConstruction("Bogacki-Shampine 3(2) Pair");
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_BogackiShampine32, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("Bogacki-Shampine 3(2) Pair", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_BogackiShampine32, AppAction)
{
  testExplicitRKAppAction("Bogacki-Shampine 3(2) Pair", out, success);
}


} // namespace Tempus_Test
