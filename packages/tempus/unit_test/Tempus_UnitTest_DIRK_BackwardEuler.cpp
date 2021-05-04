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
#include "Tempus_IntegratorBasic.hpp"

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


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_BackwardEuler, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_BackwardEuler<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_BackwardEuler, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Backward Euler", model);
}


// ************************************************************
//* Test: construct the integrator from PL and make sure that
//* the solver PL is the same as the provided solver PL
//* and not the default solver PL
// ************************************************************

TEUCHOS_UNIT_TEST(DIRK_BackwardEuler, App_PL)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // read the params from xml file
  auto pList = getParametersFromXmlFile("Tempus_DIRK_VanDerPol.xml");
  auto pl = sublist(pList, "Tempus", true);
  auto appSolverPL = pl->sublist("App Stepper").sublist("App Solver");


  // setup the Integrator
  auto integrator = Tempus::createIntegratorBasic<double>(pl, model);
  auto stepperSolverPL = Teuchos::ParameterList();
  stepperSolverPL.set("NOX", *(integrator->getStepper()->getSolver()->getParameterList()));

  // make sure the app Solver PL is being used
  TEUCHOS_ASSERT( Teuchos::haveSameValues(appSolverPL, stepperSolverPL) );

}



// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_BackwardEuler, AppAction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_BackwardEuler<double>());
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
