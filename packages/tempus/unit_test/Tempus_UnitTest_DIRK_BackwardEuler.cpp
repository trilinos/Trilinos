//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

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
  auto pList       = Teuchos::getParametersFromXmlFile("Tempus_DIRK_VanDerPol.xml");
  auto pl          = sublist(pList, "Tempus", true);
  auto appSolverPL = pl->sublist("App Stepper").sublist("App Solver");

  // setup the Integrator
  auto integrator      = Tempus::createIntegratorBasic<double>(pl, model);
  auto stepperSolverPL = Teuchos::ParameterList();
  stepperSolverPL.set(
      "NOX", *(integrator->getStepper()->getSolver()->getParameterList()));

  // make sure the app Solver PL is being used
  TEUCHOS_ASSERT(Teuchos::haveSameValues(appSolverPL, stepperSolverPL));
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_BackwardEuler, AppAction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_BackwardEuler<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
