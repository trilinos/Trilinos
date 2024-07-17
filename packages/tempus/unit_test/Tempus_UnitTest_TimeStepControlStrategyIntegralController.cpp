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
#include "Tempus_TimeStepControlStrategyIntegralController.hpp"

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
TEUCHOS_UNIT_TEST(TimeStepControlStrategyIntegralController,
                  Default_Construction)
{
  auto tscs =
      rcp(new Tempus::TimeStepControlStrategyIntegralController<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  // Test the get functions (i.e., defaults).
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getController() != "PID");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKI() != 0.58);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKP() != 0.21);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKD() != 0.10);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactor() != 0.90);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactorAfterReject() != 0.90);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMax() != 5.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMin() != 0.5);

  // Test the set functions.
  tscs->setController("I");
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setKI(0.6);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setKP(0.0);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setKD(0.0);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setSafetyFactor(0.8);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setSafetyFactorAfterReject(0.8);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setFacMax(4.0);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());
  tscs->setFacMin(0.4);
  tscs->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getController() != "I");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKI() != 0.6);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKP() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKD() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactor() != 0.8);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactorAfterReject() != 0.8);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMax() != 4.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMin() != 0.4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyIntegralController, Full_Construction)
{
  auto tscs = rcp(new Tempus::TimeStepControlStrategyIntegralController<double>(
      "I", 0.6, 0.0, 0.0, 0.8, 0.8, 4.0, 0.4));
  TEUCHOS_TEST_FOR_EXCEPT(!tscs->isInitialized());

  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getController() != "I");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKI() != 0.6);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKP() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKD() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactor() != 0.8);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactorAfterReject() != 0.8);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMax() != 4.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMin() != 0.4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyIntegralController,
                  Create_Construction)
{
  auto pl = Tempus::getTimeStepControlStrategyIntegralControllerPL<double>();

  pl->set<std::string>("Controller Type", "I");
  pl->set<double>("KI", 0.6);
  pl->set<double>("KP", 0.0);
  pl->set<double>("KD", 0.0);
  pl->set<double>("Safety Factor", 0.8);
  pl->set<double>("Safety Factor After Step Rejection", 0.8);
  pl->set<double>("Maximum Safety Factor", 4.0);
  pl->set<double>("Minimum Safety Factor", 0.4);

  auto tscs =
      Tempus::createTimeStepControlStrategyIntegralController<double>(pl);

  TEUCHOS_TEST_FOR_EXCEPT(tscs->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getController() != "I");
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKI() != 0.6);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKP() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getKD() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactor() != 0.8);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getSafetyFactorAfterReject() != 0.8);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMax() != 4.0);
  TEUCHOS_TEST_FOR_EXCEPT(tscs->getFacMin() != 0.4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyIntegralController, setNextTimeStep)
{
  double KI                      = 0.5;
  double KP                      = 0.25;
  double KD                      = 0.15;
  double safetyFactor            = 0.9;
  double safetyFactorAfterReject = 0.9;
  double facMax                  = 5.0;
  double facMin                  = 0.5;

  auto tscs = rcp(new Tempus::TimeStepControlStrategyIntegralController<double>(
      "PID", KI, KP, KD, safetyFactor, safetyFactorAfterReject, facMax,
      facMin));

  // Setup the TimeStepControl --------------------------------
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  tsc->setTimeStepControlStrategy(tscs);
  tsc->setInitTime(0.0);
  tsc->setFinalTime(10.0);
  tsc->setMinTimeStep(0.01);
  tsc->setInitTimeStep(1.0);
  tsc->setMaxTimeStep(10.0);
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

  double order = 2.0;
  solutionHistory->addState(icState);
  solutionHistory->getCurrentState()->setTimeStep(1.0);
  solutionHistory->getCurrentState()->setTime(0.0);
  solutionHistory->getCurrentState()->setIndex(0);
  solutionHistory->getCurrentState()->setOrder(order);

  // Mock Integrator

  // -- First Time Step
  solutionHistory->initWorkingState();
  auto currentState = solutionHistory->getCurrentState();
  auto workingState = solutionHistory->getWorkingState();

  TEST_FLOATING_EQUALITY(workingState->getErrorRel(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm1(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm2(), 0.0, 1.0e-14);

  tsc->setNextTimeStep(solutionHistory, status);

  // First time step should cause no change to dt because
  // internal relative errors = 1.
  TEST_FLOATING_EQUALITY(workingState->getTimeStep(), 1.0, 1.0e-14);

  // Mock takeStep
  double errN = 0.1;
  workingState->setErrorRel(errN);
  workingState->setSolutionStatus(Tempus::Status::PASSED);

  // -- Second Time Step
  solutionHistory->initWorkingState();
  currentState = solutionHistory->getCurrentState();
  workingState = solutionHistory->getWorkingState();
  double dt    = workingState->getTimeStep();

  TEST_FLOATING_EQUALITY(workingState->getErrorRel(), 0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm1(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm2(), 0.0, 1.0e-14);

  tsc->setNextTimeStep(solutionHistory, status);

  double p     = order - 1.0;
  double dtNew = dt * safetyFactor * std::pow(errN, -KI / p);
  TEST_FLOATING_EQUALITY(workingState->getTimeStep(), dtNew, 1.0e-14);

  // Mock takeStep
  errN          = 0.2;
  double errNm1 = 0.1;
  workingState->setErrorRel(errN);
  workingState->setSolutionStatus(Tempus::Status::PASSED);

  // -- Third Time Step
  solutionHistory->initWorkingState();
  currentState = solutionHistory->getCurrentState();
  workingState = solutionHistory->getWorkingState();
  dt           = workingState->getTimeStep();

  TEST_FLOATING_EQUALITY(workingState->getErrorRel(), 0.2, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm1(), 0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm2(), 0.0, 1.0e-14);

  tsc->setNextTimeStep(solutionHistory, status);

  dtNew =
      dt * safetyFactor * std::pow(errN, -KI / p) * std::pow(errNm1, KP / p);
  TEST_FLOATING_EQUALITY(workingState->getTimeStep(), dtNew, 1.0e-14);

  // Mock takeStep
  errN          = 0.3;
  errNm1        = 0.2;
  double errNm2 = 0.1;
  workingState->setErrorRel(errN);
  workingState->setSolutionStatus(Tempus::Status::PASSED);

  // -- Fourth Time Step
  solutionHistory->initWorkingState();
  currentState = solutionHistory->getCurrentState();
  workingState = solutionHistory->getWorkingState();
  dt           = workingState->getTimeStep();

  TEST_FLOATING_EQUALITY(workingState->getErrorRel(), 0.3, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm1(), 0.2, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getErrorRelNm2(), 0.1, 1.0e-14);

  tsc->setNextTimeStep(solutionHistory, status);

  dtNew = dt * safetyFactor * std::pow(errN, -KI / p) *
          std::pow(errNm1, KP / p) * std::pow(errNm2, -KD / p);
  TEST_FLOATING_EQUALITY(workingState->getTimeStep(), dtNew, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyIntegralController, getValidParameters)
{
  auto tscs =
      rcp(new Tempus::TimeStepControlStrategyIntegralController<double>());

  auto pl = tscs->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Strategy Type"), ==,
               "Integral Controller");
  TEST_COMPARE(pl->get<std::string>("Controller Type"), ==, "PID");
  TEST_FLOATING_EQUALITY(pl->get<double>("KI"), 0.58, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("KP"), 0.21, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("KD"), 0.10, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Safety Factor"), 0.9, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Safety Factor After Step Rejection"),
                         0.9, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Maximum Safety Factor"), 5.0,
                         1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Minimum Safety Factor"), 0.5,
                         1.0e-14);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

}  // namespace Tempus_Unit_Test
