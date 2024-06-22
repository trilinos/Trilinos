//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategyConstant.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyConstant, Default_Construction)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyConstant<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  // Test the get functions (i.e., defaults).
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStrategyType() != "Constant");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Constant");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getName() != "Constant");
  TEST_FLOATING_EQUALITY(tscs->getConstantTimeStep(), 0.0, 1.0e-14);

  // Test the set functions.
  tscs->setConstantTimeStep(0.989);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  TEST_FLOATING_EQUALITY(tscs->getConstantTimeStep(), 0.989, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyConstant, Full_Construction)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyConstant<double>(
      0.123, "Full_Construction_Test"));
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStrategyType() != "Constant");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Constant");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getName() != "Full_Construction_Test");
  TEST_FLOATING_EQUALITY(tscs->getConstantTimeStep(), 0.123, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyConstant, Create_Construction)
{
  auto pl = Tempus::getTimeStepControlStrategyConstantPL<double>();

  // Set strategy parameters.
  pl->set<double>("Time Step", 0.02);

  auto tscsc = Tempus::createTimeStepControlStrategyConstant<double>(pl);
  TEUCHOS_TEST_FOR_EXCEPT(!tscsc->isInitialized());

  TEST_FLOATING_EQUALITY(tscsc->getConstantTimeStep(), 0.02, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyConstant, setNextTimeStep)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyConstant<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  double initTime = 1.0;
  int initIndex   = -100;

  // Setup the SolutionHistory --------------------------------
  auto model    = rcp(new Tempus_Test::SinCosModel<double>());
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX<double>(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);
  solutionHistory->getCurrentState()->setTimeStep(0.9);
  solutionHistory->getCurrentState()->setTime(initTime);
  solutionHistory->getCurrentState()->setIndex(initIndex);

  // Setup the TimeStepControl --------------------------------
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  tsc->setTimeStepControlStrategy(tscs);
  tsc->setInitTime(initTime);
  tsc->setFinalTime(100.0);
  tsc->setMinTimeStep(0.01);
  tsc->setInitTimeStep(0.02);
  tsc->setMaxTimeStep(0.05);
  tsc->setInitIndex(initIndex);
  tsc->setFinalIndex(100);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  Tempus::Status status = Tempus::Status::WORKING;

  // ** Mock Timestep ** //
  solutionHistory->initWorkingState();

  tsc->setNextTimeStep(solutionHistory, status);
  // ** Mock Timestep ** //

  auto workingState = solutionHistory->getWorkingState();
  TEST_FLOATING_EQUALITY(workingState->getTimeStep(), 0.02, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getTime(), 1.02, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyConstant, getValidParameters)
{
  auto tscsc = rcp(new Tempus::TimeStepControlStrategyConstant<double>());

  auto pl = tscsc->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Strategy Type"), ==, "Constant");
  TEST_FLOATING_EQUALITY(pl->get<double>("Time Step"), 0.0, 1.0e-14);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

}  // namespace Tempus_Unit_Test
