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
TEUCHOS_UNIT_TEST(SDIRK_2Stage2ndOrder, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperSDIRK_2Stage2ndOrder<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
  TEUCHOS_ASSERT(stepper->getGamma() == 0.2928932188134524);
  stepper->setGamma(0.5); stepper->initialize(); TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_2Stage2ndOrder, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("SDIRK 2 Stage 2nd order", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SDIRK_2Stage2ndOrder, AppAction)
{
  auto stepper = rcp(new Tempus::StepperSDIRK_2Stage2ndOrder<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
