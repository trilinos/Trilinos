//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"

#include "../TestModels/DahlquistTestModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_ForwardEuler, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>());
  testExplicitRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
  const auto rk_fe = stepper->getTableau();

  TEUCHOS_ASSERT(rk_fe->isTVD());
  TEUCHOS_ASSERT(rk_fe->getTVDCoeff() == 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_ForwardEuler, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Forward Euler", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_ForwardEuler, AppAction)
{
  auto stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_ForwardEuler, FSAL)
{
  auto stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>());
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, true));

  stepper->setModel(model);
  stepper->setICConsistency("Consistent");
  stepper->setUseFSAL(true);
  stepper->initialize();

  // Create a SolutionHistory.
  auto solutionHistory = Tempus::createSolutionHistoryME(model);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 1.0;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  solutionHistory->getWorkingState()->setTime(dt);
  stepper->takeStep(solutionHistory);

  // Test solution.
  const double relTol = 1.0e-14;

  // ICs
  auto currentState   = solutionHistory->getCurrentState();
  const double x_0    = get_ele(*(currentState->getX()), 0);
  const double xDot_0 = get_ele(*(currentState->getXDot()), 0);
  TEST_FLOATING_EQUALITY(x_0, 1.0, relTol);
  TEST_FLOATING_EQUALITY(xDot_0, -1.0, relTol);
  TEST_ASSERT(std::abs(currentState->getTime()) < relTol);

  // After one step.
  auto workingState   = solutionHistory->getWorkingState();
  const double x_1    = get_ele(*(workingState->getX()), 0);
  const double xDot_1 = get_ele(*(workingState->getXDot()), 0);
  // out << "xDot_1 = " << xDot_1 << std::endl;
  TEST_ASSERT(std::abs(x_1) < relTol);
  TEST_ASSERT(std::abs(xDot_1) < relTol);
  TEST_FLOATING_EQUALITY(workingState->getTime(), 1.0, relTol);
}

}  // namespace Tempus_Unit_Test
