//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"

#include "Tempus_StepperSubcycling.hpp"
#include "Tempus_StepperSubcyclingModifierBase.hpp"
#include "Tempus_StepperSubcyclingModifierXBase.hpp"
#include "Tempus_StepperSubcyclingObserverBase.hpp"
#include "Tempus_StepperSubcyclingModifierDefault.hpp"
#include "Tempus_StepperSubcyclingModifierXDefault.hpp"
#include "Tempus_StepperSubcyclingObserverDefault.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

using Tempus::StepperFactory;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, Default_Construction)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double> >(model);

  // Setup SolutionHistory ------------------------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);
  solutionHistory->initWorkingState();

  // Default construction.
  auto stepper   = rcp(new Tempus::StepperSubcycling<double>());
  auto stepperBE = Tempus::createStepperBackwardEuler(modelME, Teuchos::null);
  stepper->setSubcyclingStepper(stepperBE);
  stepper->setInitialConditions(solutionHistory);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier  = rcp(new Tempus::StepperSubcyclingModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperSubcyclingModifierXDefault<double>());
  auto observer  = rcp(new Tempus::StepperSubcyclingObserverDefault<double>());
  auto solver    = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();

  // Test the set functions
  stepper->setSolver(solver);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
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

  // Full argument list construction.
  auto scIntegrator = Teuchos::rcp(new Tempus::IntegratorBasic<double>());
  auto stepperFE    = Tempus::createStepperForwardEuler(modelME, Teuchos::null);
  scIntegrator->setStepper(stepperFE);
  scIntegrator->setSolutionHistory(solutionHistory);
  scIntegrator->initialize();

  stepper = rcp(new Tempus::StepperSubcycling<double>(
      model, scIntegrator, useFSAL, ICConsistency, ICConsistencyCheck,
      modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, MaxTimeStepDoesNotChangeDuring_takeStep)
{
  // Setup the stepper ----------------------------------------
  auto model     = rcp(new Tempus_Test::SinCosModel<double>());
  auto modelME   = rcp_dynamic_cast<const Thyra::ModelEvaluator<double> >(model);
  auto stepper   = rcp(new Tempus::StepperSubcycling<double>());
  auto stepperBE = Tempus::createStepperBackwardEuler(modelME, Teuchos::null);
  stepper->setSubcyclingStepper(stepperBE);

  // Setup SolutionHistory ------------------------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);
  solutionHistory->initWorkingState();

  // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
  stepper->setInitialConditions(solutionHistory);
  stepper->initialize();

  // Test
  stepper->setSubcyclingInitTimeStep(0.25);
  stepper->setSubcyclingMaxTimeStep(0.5);
  double maxTimeStep_Set = stepper->getSubcyclingMaxTimeStep();
  stepper->takeStep(solutionHistory);
  double maxTimeStep_After = stepper->getSubcyclingMaxTimeStep();

  TEST_FLOATING_EQUALITY(maxTimeStep_Set, maxTimeStep_After, 1.0e-14);
}

// ************************************************************
// ************************************************************
class StepperSubcyclingModifierTest
  : virtual public Tempus::StepperSubcyclingModifierBase<double> {
 public:
  /// Constructor
  StepperSubcyclingModifierTest()
    : testBEGIN_STEP(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(1.5),
      testName("")
  {
  }
  /// Destructor
  virtual ~StepperSubcyclingModifierTest() {}

  /// Modify Subcycling Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperSubcycling<double> > stepper,
      const typename Tempus::StepperSubcyclingAppAction<double>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperSubcyclingAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        testName         = "Subcycling - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperSubcyclingAppAction<double>::END_STEP: {
        testEND_STEP     = true;
        auto x           = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        testDt           = sh->getWorkingState()->getTimeStep() / 10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }
  }

  bool testBEGIN_STEP;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(Subcycling, AppAction_Modifier)
{
  // Setup the SinCosModel ------------------------------------
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double> >(model);

  // Setup Stepper for field solve ----------------------------
  auto stepper   = rcp(new Tempus::StepperSubcycling<double>());
  auto stepperFE = Tempus::createStepperForwardEuler(modelME, Teuchos::null);
  auto modifier  = rcp(new StepperSubcyclingModifierTest());
  stepper->setAppAction(modifier);
  stepper->setSubcyclingStepper(stepperFE);

  stepper->setSubcyclingMinTimeStep(15);
  stepper->setSubcyclingInitTimeStep(15.0);
  stepper->setSubcyclingMaxTimeStep(15.0);
  stepper->setSubcyclingMaxFailures(10);
  stepper->setSubcyclingMaxConsecFailures(5);
  stepper->setSubcyclingScreenOutputIndexInterval(1);
  stepper->setSubcyclingPrintDtChanges(true);

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime(0.0);
  timeStepControl->setFinalTime(1.0);
  timeStepControl->setInitTimeStep(15.0);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(timeStepControl->getInitTime());
  icState->setIndex(timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);                           // dt for ICs are indicated by zero.
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
  stepper->setInitialConditions(solutionHistory);
  stepper->initialize();

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  solutionHistory->getWorkingState()->setTimeStep(15.0);
  stepper->takeStep(solutionHistory);

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-14);

  TEST_COMPARE(modifier->testName, ==, "Subcycling - Modifier");
}

// ************************************************************
// ************************************************************
class StepperSubcyclingObserverTest
  : virtual public Tempus::StepperSubcyclingObserverBase<double> {
 public:
  /// Constructor
  StepperSubcyclingObserverTest()
    : testBEGIN_STEP(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(15.0),
      testName("Subcyling")
  {
  }

  /// Destructor
  virtual ~StepperSubcyclingObserverTest() {}

  /// Observe Subcycling Stepper at action location.
  virtual void observe(
      Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<const Tempus::StepperSubcycling<double> > stepper,
      const typename Tempus::StepperSubcyclingAppAction<double>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperSubcyclingAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperSubcyclingAppAction<double>::END_STEP: {
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
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(Subcycling, AppAction_Observer)
{
  // Setup the SinCosModel ------------------------------------
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double> >(model);

  // Setup Stepper for field solve ----------------------------
  auto stepper   = rcp(new Tempus::StepperSubcycling<double>());
  auto stepperFE = Tempus::createStepperForwardEuler(modelME, Teuchos::null);
  auto observer  = rcp(new StepperSubcyclingObserverTest());
  stepper->setAppAction(observer);
  stepper->setSubcyclingStepper(stepperFE);

  stepper->setSubcyclingMinTimeStep(15);
  stepper->setSubcyclingInitTimeStep(15.0);
  stepper->setSubcyclingMaxTimeStep(15.0);
  stepper->setSubcyclingMaxFailures(10);
  stepper->setSubcyclingMaxConsecFailures(5);
  stepper->setSubcyclingScreenOutputIndexInterval(1);
  stepper->setSubcyclingPrintDtChanges(true);

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime(0.0);
  timeStepControl->setFinalTime(1.0);
  timeStepControl->setInitTimeStep(15.0);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(timeStepControl->getInitTime());
  icState->setIndex(timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);                           // dt for ICs are indicated by zero.
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
  stepper->setInitialConditions(solutionHistory);
  stepper->initialize();

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  solutionHistory->getWorkingState()->setTimeStep(15.0);
  stepper->takeStep(solutionHistory);

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
  TEST_COMPARE(observer->testEND_STEP, ==, true);

  // Testing that values can be observed through the observer.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  TEST_FLOATING_EQUALITY(observer->testDt, 15.0, 1.0e-14);

  TEST_COMPARE(observer->testName, ==, "Subcyling");
}

// ************************************************************
// ************************************************************
class StepperSubcyclingModifierXTest
  : virtual public Tempus::StepperSubcyclingModifierXBase<double> {
 public:
  /// Constructor
  StepperSubcyclingModifierXTest()
    : testX_BEGIN_STEP(false),
      testXDOT_END_STEP(false),
      testX(-0.99),
      testXDot(-0.99),
      testDt(1.5),
      testTime(1.5)
  {
  }

  /// Destructor
  virtual ~StepperSubcyclingModifierXTest() {}

  /// Modify Subcycling Stepper at action location.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<double> > x,
                      const double time, const double dt,
                      const typename Tempus::StepperSubcyclingModifierXBase<
                          double>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperSubcyclingModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testX            = get_ele(*(x), 0);
        testDt           = dt;
        break;
      }
      case StepperSubcyclingModifierXBase<double>::XDOT_END_STEP: {
        testXDOT_END_STEP = true;
        testXDot          = get_ele(*(x), 0);
        testTime          = time;
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }
  }

  bool testX_BEGIN_STEP;
  bool testXDOT_END_STEP;
  double testX;
  double testXDot;
  double testDt;
  double testTime;
};

TEUCHOS_UNIT_TEST(Subcycling, AppAction_ModifierX)
{
  // Setup the SinCosModel ------------------------------------
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double> >(model);

  // Setup Stepper for field solve ----------------------------
  auto stepper   = rcp(new Tempus::StepperSubcycling<double>());
  auto stepperFE = Tempus::createStepperForwardEuler(modelME, Teuchos::null);
  auto modifierX = rcp(new StepperSubcyclingModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->setSubcyclingStepper(stepperFE);

  stepper->setSubcyclingMinTimeStep(15);
  stepper->setSubcyclingInitTimeStep(15.0);
  stepper->setSubcyclingMaxTimeStep(15.0);
  stepper->setSubcyclingMaxFailures(10);
  stepper->setSubcyclingMaxConsecFailures(5);
  stepper->setSubcyclingScreenOutputIndexInterval(1);
  stepper->setSubcyclingPrintDtChanges(true);

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime(0.0);
  timeStepControl->setFinalTime(1.0);
  timeStepControl->setInitTimeStep(15.0);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icSolutionDot =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x_dot());
  auto icState = Tempus::createSolutionStateX(icSolution, icSolutionDot);
  icState->setTime(timeStepControl->getInitTime());
  icState->setIndex(timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);                           // dt for ICs are indicated by zero.
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
  stepper->setInitialConditions(solutionHistory);
  stepper->initialize();

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  solutionHistory->getWorkingState()->setTimeStep(15.0);
  stepper->takeStep(solutionHistory);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  solutionHistory->getWorkingState()->setTimeStep(15.0);
  stepper->takeStep(solutionHistory);

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testXDOT_END_STEP, ==, true);

  // Testing that values can be set through the Modifier.
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
