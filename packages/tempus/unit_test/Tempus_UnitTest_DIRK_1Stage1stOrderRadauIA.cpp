// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_UnitTest_Utils.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::StepperFactory;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1Stage1stOrderRadauIA, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_1Stage1stOrderRadauIA<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1Stage1stOrderRadauIA, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Implicit 1 Stage 1st order Radau IA", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1Stage1stOrderRadauIA, AppAction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_1Stage1stOrderRadauIA<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
