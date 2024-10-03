//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"
#include "Tempus_StepperRKButcherTableau.hpp"

#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperForwardEulerModifierBase.hpp"
#include "Tempus_StepperForwardEulerModifierXBase.hpp"
#include "Tempus_StepperForwardEulerObserverBase.hpp"
#include "Tempus_StepperForwardEulerModifierDefault.hpp"
#include "Tempus_StepperForwardEulerModifierXDefault.hpp"
#include "Tempus_StepperForwardEulerObserverDefault.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

using Tempus::StepperExplicitRK;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, Default_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto modifier = rcp(new Tempus::StepperForwardEulerModifierDefault<double>());
  auto modifierX =
      rcp(new Tempus::StepperForwardEulerModifierXDefault<double>());
  auto observer = rcp(new Tempus::StepperForwardEulerObserverDefault<double>());
  auto stepper  = rcp(new Tempus::StepperForwardEuler<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  // Default values for construction.
  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();

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
  stepper->setUseFSAL(useFSAL);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  stepper = rcp(new Tempus::StepperForwardEuler<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("Forward Euler", model);
}

// ************************************************************
// ************************************************************
class StepperForwardEulerModifierTest
  : virtual public Tempus::StepperForwardEulerModifierBase<double> {
 public:
  /// Constructor
  StepperForwardEulerModifierTest()
    : testBEGIN_STEP(false),
      testBEFORE_EXPLICIT_EVAL(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testName("")
  {
  }

  /// Destructor
  virtual ~StepperForwardEulerModifierTest() {}

  /// Observe ForwardEuler Stepper at end of takeStep.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperForwardEuler<double> > stepper,
      const typename Tempus::StepperForwardEulerAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperForwardEulerAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperForwardEulerAppAction<double>::BEFORE_EXPLICIT_EVAL: {
        testBEFORE_EXPLICIT_EVAL = true;
        testDt                   = sh->getWorkingState()->getTimeStep() / 10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        testName = "Forward Euler - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperForwardEulerAppAction<double>::END_STEP: {
        testEND_STEP     = true;
        auto x           = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }
  }
  bool testBEGIN_STEP;
  bool testBEFORE_EXPLICIT_EVAL;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(ForwardEuler, AppAction_Modifier)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperForwardEuler<double>());
  stepper->setModel(model);
  auto modifier = rcp(new StepperForwardEulerModifierTest());
  stepper->setAppAction(modifier);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(0.0);
  icState->setIndex(0);
  icState->setTimeStep(0.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 0.1;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-14);

  TEST_COMPARE(modifier->testName, ==, "Forward Euler - Modifier");
}

// ************************************************************
// ************************************************************
class StepperForwardEulerObserverTest
  : virtual public Tempus::StepperForwardEulerObserverBase<double> {
 public:
  /// Constructor
  StepperForwardEulerObserverTest()
    : testBEGIN_STEP(false),
      testBEFORE_EXPLICIT_EVAL(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testName("")
  {
  }

  /// Destructor
  virtual ~StepperForwardEulerObserverTest() {}

  /// Observe ForwardEuler Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<const Tempus::StepperForwardEuler<double> > stepper,
      const typename Tempus::StepperForwardEulerAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperForwardEulerAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperForwardEulerAppAction<double>::BEFORE_EXPLICIT_EVAL: {
        testBEFORE_EXPLICIT_EVAL = true;
        testDt                   = sh->getWorkingState()->getTimeStep();
        testName                 = stepper->getStepperName();
        break;
      }
      case StepperForwardEulerAppAction<double>::END_STEP: {
        testEND_STEP     = true;
        auto x           = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }
  }

  bool testBEGIN_STEP;
  bool testBEFORE_EXPLICIT_EVAL;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(ForwardEuler, AppAction_Observer)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperForwardEuler<double>());
  stepper->setModel(model);
  auto observer = rcp(new StepperForwardEulerObserverTest());
  stepper->setAppAction(observer);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(0.0);
  icState->setIndex(0);
  icState->setTimeStep(0.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 0.1;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
  TEST_COMPARE(observer->testBEFORE_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(observer->testEND_STEP, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  TEST_FLOATING_EQUALITY(observer->testDt, dt, 1.0e-14);

  TEST_COMPARE(observer->testName, ==, "Forward Euler");
}

// ************************************************************
// ************************************************************
class StepperForwardEulerModifierXTest
  : virtual public Tempus::StepperForwardEulerModifierXBase<double> {
 public:
  /// Constructor
  StepperForwardEulerModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEFORE_EXPLICIT_EVAL(false),
      testXDOT_END_STEP(false),
      testX(-0.99),
      testXDot(-0.99),
      testDt(-1.5),
      testTime(-1.5)
  {
  }

  /// Destructor
  virtual ~StepperForwardEulerModifierXTest() {}

  /// Observe BackwardEuler Stepper at end of takeStep.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<double> > x,
                      const double time, const double dt,
                      const typename Tempus::StepperForwardEulerModifierXBase<
                          double>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperForwardEulerModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testX            = get_ele(*(x), 0);
        break;
      }
      case StepperForwardEulerModifierXBase<double>::X_BEFORE_EXPLICIT_EVAL: {
        testX_BEFORE_EXPLICIT_EVAL = true;
        testDt                     = dt;
        testTime                   = time;
        break;
      }
      case StepperForwardEulerModifierXBase<double>::XDOT_END_STEP: {
        testXDOT_END_STEP = true;
        testXDot          = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }
  }
  bool testX_BEGIN_STEP;
  bool testX_BEFORE_EXPLICIT_EVAL;
  bool testXDOT_END_STEP;
  double testX;
  double testXDot;
  double testDt;
  double testTime;
};

TEUCHOS_UNIT_TEST(ForwardEuler, AppAction_ModifierX)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperForwardEuler<double>());
  stepper->setModel(model);
  auto modifierX = rcp(new StepperForwardEulerModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(0.0);
  icState->setIndex(0);
  icState->setTimeStep(0.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 0.1;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifierX->testXDOT_END_STEP, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(x), 0), 1.0e-14);
  // Temporary memory for xDot is not guarranteed to exist outside the Stepper.
  auto xDot = solutionHistory->getWorkingState()->getXDot();
  if (xDot == Teuchos::null) xDot = stepper->getStepperXDot();

  TEST_FLOATING_EQUALITY(modifierX->testXDot, get_ele(*(xDot), 0), 1.0e-14);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-14);

  auto time = solutionHistory->getWorkingState()->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-14);
}

}  // namespace Tempus_Unit_Test
