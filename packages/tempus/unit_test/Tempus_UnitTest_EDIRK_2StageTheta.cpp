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
TEUCHOS_UNIT_TEST(EDIRK_2StageTheta, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperEDIRK_2StageTheta<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
  double theta = 0.5;
  TEUCHOS_ASSERT(stepper->getTheta() == theta);
  stepper->setTheta(theta); stepper->initialize(); TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_2StageTheta, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("EDIRK 2 Stage Theta Method", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_2StageTheta, AppAction)
{
  auto stepper = rcp(new Tempus::StepperEDIRK_2StageTheta<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
