//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_General, Default_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperDIRK_General<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier  = rcp(new Tempus::StepperRKModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperRKModifierXDefault<double>());
  auto observer  = rcp(new Tempus::StepperRKObserverDefault<double>());
  auto solver    = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();
  bool useEmbedded          = stepper->getUseEmbedded();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();

  using Teuchos::as;
  int NumStages = 2;
  Teuchos::SerialDenseMatrix<int, double> A(NumStages, NumStages);
  Teuchos::SerialDenseVector<int, double> b(NumStages);
  Teuchos::SerialDenseVector<int, double> c(NumStages);
  Teuchos::SerialDenseVector<int, double> bstar(0);

  // Fill A:
  A(0, 0) = 0.2928932188134524;
  A(0, 1) = 0.0;
  A(1, 0) = 0.7071067811865476;
  A(1, 1) = 0.2928932188134524;

  // Fill b:
  b(0) = 0.7071067811865476;
  b(1) = 0.2928932188134524;

  // Fill c:
  c(0) = 0.2928932188134524;
  c(1) = 1.0;

  int order = 2;

  // Test the set functions.
  stepper->setAppAction(modifier);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(modifierX);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(observer);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setSolver(solver);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseEmbedded(useEmbedded);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setZeroInitialGuess(zeroInitialGuess);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  stepper->setTableau(A, b, c, order, order, order);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Full argument list construction.
  stepper = rcp(new Tempus::StepperDIRK_General<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded,
      zeroInitialGuess, modifier, A, b, c, order, order, order, bstar));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_General, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("General DIRK", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_General, AppAction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_General<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
