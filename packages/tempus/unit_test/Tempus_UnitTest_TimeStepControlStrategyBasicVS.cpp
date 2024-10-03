//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"

#include "../TestModels/DahlquistTestModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyBasicVS, Default_Construction)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  // Test the get functions (i.e., defaults).
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStrategyType() != "Basic VS");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getAmplFactor() != 1.75);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getReductFactor() != 0.5);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMinEta() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMaxEta() != 1.0e+16);

  // Test the set functions.
  tscs->setAmplFactor(1.33);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setReductFactor(0.75);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setMinEta(0.01);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setMaxEta(0.05);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  TEUCHOS_TEST_FOR_EXCEPT(tscs->getAmplFactor() != 1.33);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getReductFactor() != 0.75);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMinEta() != 0.01);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMaxEta() != 0.05);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyBasicVS, Full_Construction)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>(
      1.33, 0.75, 0.01, 0.05));
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStrategyType() != "Basic VS");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getAmplFactor() != 1.33);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getReductFactor() != 0.75);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMinEta() != 0.01);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMaxEta() != 0.05);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyBasicVS, Create_Construction)
{
  auto pl = Tempus::getTimeStepControlStrategyBasicVS_PL<double>();

  pl->set<double>("Amplification Factor", 1.33);
  pl->set<double>("Reduction Factor", 0.75);
  pl->set<double>("Minimum Value Monitoring Function", 0.01);
  pl->set<double>("Maximum Value Monitoring Function", 0.05);

  auto tscs = Tempus::createTimeStepControlStrategyBasicVS<double>(pl);

  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStrategyType() != "Basic VS");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getAmplFactor() != 1.33);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getReductFactor() != 0.75);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMinEta() != 0.01);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getMaxEta() != 0.05);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyBasicVS, setNextTimeStep)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  tscs->setAmplFactor(1.1);
  tscs->setReductFactor(0.5);
  tscs->setMinEta(0.01);
  tscs->setMaxEta(0.05);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  // Setup the TimeStepControl --------------------------------
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  tsc->setTimeStepControlStrategy(tscs);
  tsc->setInitTime(0.0);
  tsc->setFinalTime(10.0);
  tsc->setMinTimeStep(0.01);
  tsc->setInitTimeStep(0.1);
  tsc->setMaxTimeStep(1.0);
  tsc->setFinalIndex(100);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  Tempus::Status status = Tempus::Status::WORKING;

  // Setup the SolutionHistory --------------------------------
  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>());
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX<double>(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());

  // Test Reducing timestep
  {
    solutionHistory->addState(icState);
    solutionHistory->getCurrentState()->setTimeStep(0.5);
    solutionHistory->getCurrentState()->setTime(0.0);
    solutionHistory->getCurrentState()->setIndex(0);

    // Set up solution history with two time steps.
    for (int i = 0; i < 2; i++) {
      solutionHistory->initWorkingState();

      tsc->setNextTimeStep(solutionHistory, status);

      {  // Mock takeStep
        auto currentState = solutionHistory->getCurrentState();
        auto workingState = solutionHistory->getWorkingState();
        auto xN           = workingState->getX();
        Thyra::Vp_S(xN.ptr(), 1.0);
        workingState->setSolutionStatus(Tempus::Status::PASSED);
        workingState->computeNorms(currentState);
      }

      solutionHistory->promoteWorkingState();
      // out << "  x = " << get_ele(*(x), 0) << std::endl;
    }

    auto currentState = solutionHistory->getCurrentState();
    TEST_FLOATING_EQUALITY(currentState->getTimeStep(), 0.25, 1.0e-14);
    TEST_FLOATING_EQUALITY(currentState->getTime(), 0.75, 1.0e-14);
  }

  // Test increasing timestep
  {
    solutionHistory->clear();
    solutionHistory->addState(icState);
    solutionHistory->getCurrentState()->setTimeStep(0.5);
    solutionHistory->getCurrentState()->setTime(0.0);
    solutionHistory->getCurrentState()->setIndex(0);

    // Set up solution history with two time steps.
    for (int i = 0; i < 2; i++) {
      solutionHistory->initWorkingState();

      tsc->setNextTimeStep(solutionHistory, status);

      {  // Mock takeStep
        auto currentState = solutionHistory->getCurrentState();
        auto workingState = solutionHistory->getWorkingState();
        auto xN           = workingState->getX();
        Thyra::Vp_S(xN.ptr(), 0.0);
        workingState->setSolutionStatus(Tempus::Status::PASSED);
        workingState->computeNorms(currentState);
      }

      solutionHistory->promoteWorkingState();
      // auto x = solutionHistory->getCurrentState()->getX();
      // out << "  x = " << get_ele(*(x), 0) << std::endl;
    }

    auto currentState = solutionHistory->getCurrentState();
    TEST_FLOATING_EQUALITY(currentState->getTimeStep(), 0.55, 1.0e-14);
    TEST_FLOATING_EQUALITY(currentState->getTime(), 1.05, 1.0e-14);
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyBasicVS, getValidParameters)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());

  auto pl = tscs->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Strategy Type"), ==, "Basic VS");
  TEST_FLOATING_EQUALITY(pl->get<double>("Amplification Factor"), 1.75,
                         1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Reduction Factor"), 0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Minimum Value Monitoring Function"),
                         0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Maximum Value Monitoring Function"),
                         1.0e+16, 1.0e-14);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

}  // namespace Tempus_Unit_Test
