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

#include "Tempus_StepperLeapfrogModifierBase.hpp"
#include "Tempus_StepperLeapfrogObserverBase.hpp"
#include "Tempus_StepperLeapfrogModifierXBase.hpp"
#include "Tempus_StepperLeapfrogModifierDefault.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestModels/HarmonicOscillatorModel.hpp"
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
TEUCHOS_UNIT_TEST(Leapfrog, Default_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperLeapfrog<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Default values for construction.
  auto modifier = rcp(new Tempus::StepperLeapfrogModifierDefault<double>());
  stepper->setAppAction(modifier);
  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();


  // Test the set functions.
  stepper->setAppAction(modifier);                      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Full argument list construction.
  stepper = rcp(new Tempus::StepperLeapfrog<double>(
    model, useFSAL, ICConsistency, ICConsistencyCheck,modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Leapfrog, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());
  testFactoryConstruction("Leapfrog", model);
}

// ************************************************************
// ************************************************************
class StepperLeapfrogModifierTest
  : virtual public Tempus::StepperLeapfrogModifierBase<double>
{
public:

  /// Constructor
  StepperLeapfrogModifierTest()
    : testBEGIN_STEP(false), testBEFORE_X_UPDATE(false),testBEFORE_EXPLICIT_EVAL(false),
      testBEFORE_XDOT_UPDATE(false), testCurrentValue(-0.99), testWorkingValue(-0.99),
      testDt(-1.5), testType("")
  {}

  /// Destructor
  virtual ~StepperLeapfrogModifierTest(){}

  /// Observe Leapfrog Stepper at end of takeStep.
  virtual void modify(
                      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
                      Teuchos::RCP<Tempus::StepperLeapfrog<double> > stepper,
                      const typename Tempus::StepperLeapfrogAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
    case StepperLeapfrogAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
    case StepperLeapfrogAppAction<double>::BEFORE_EXPLICIT_EVAL:
      {
        testBEFORE_EXPLICIT_EVAL = true;
        testDt = sh->getWorkingState()->getTimeStep()/10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
    case StepperLeapfrogAppAction<double>::BEFORE_X_UPDATE:
      {
        testBEFORE_X_UPDATE = true;
        testType = "Leapfrog - Modifier";
        stepper->setStepperType(testType);
        break;
      }
    case StepperLeapfrogAppAction<double>::BEFORE_XDOT_UPDATE:
      {
        testBEFORE_XDOT_UPDATE = true;
        auto x = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        break;
      }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "Error - unknown action location.\n");
    }
  }
  bool testBEGIN_STEP;
  bool testBEFORE_X_UPDATE;
  bool testBEFORE_EXPLICIT_EVAL;
  bool testBEFORE_XDOT_UPDATE;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};

TEUCHOS_UNIT_TEST(Leapfrog, AppAction_Modifier)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperLeapfrog<double>());
  stepper->setModel(model);
  auto modifier = rcp(new StepperLeapfrogModifierTest());
  stepper->setAppAction(modifier);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitTimeStep(15.0);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  auto icX = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icXDot = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
  auto icXDotDot = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
  auto icState = Tempus::createSolutionStateX<double>(icX, icXDot, icXDotDot);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(15.0);
  icState->setOrder   (stepper->getOrder());
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
  solutionHistory->getWorkingState()->setTimeStep(15.0);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifier->testBEFORE_X_UPDATE, ==, true);
  TEST_COMPARE(modifier->testBEFORE_XDOT_UPDATE, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-15);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(x), 0), 1.0e-15);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-15);

  TEST_COMPARE(modifier->testType, ==, "Leapfrog - Modifier");
}
// ************************************************************
// ************************************************************
class StepperLeapfrogModifierXTest
  : virtual public Tempus::StepperLeapfrogModifierXBase<double>
{
public:

  /// Constructor
  StepperLeapfrogModifierXTest()
    : testX_BEGIN_STEP(false), testX_BEFORE_EXPLICIT_EVAL(false),
      testX_BEFORE_X_UPDATE(false), testX_BEFORE_XDOT_UPDATE(false),
      testX(0.0), testDt(-1.25), testTime(-1.25),testType("")
  {}

  /// Destructor
  virtual ~StepperLeapfrogModifierXTest(){}

  /// Observe Leapfrog Stepper at end of takeStep.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<double> > x,
    const double time, const double dt,
    const typename Tempus::StepperLeapfrogModifierXBase<double>::MODIFIER_TYPE modType)
  {
    switch(modType) {
    case StepperLeapfrogModifierXBase<double>::X_BEGIN_STEP:
      {
        testX_BEGIN_STEP = true;
        testX = get_ele(*(x), 0);
        break;
      }
    case StepperLeapfrogModifierXBase<double>::X_BEFORE_EXPLICIT_EVAL:
      {
        testX_BEFORE_EXPLICIT_EVAL = true;
        testDt = dt;
        break;
      }
    case StepperLeapfrogModifierXBase<double>::X_BEFORE_X_UPDATE:
      {
        testX_BEFORE_X_UPDATE = true;
        testTime = time;
        break;
      }
    case StepperLeapfrogModifierXBase<double>::X_BEFORE_XDOT_UPDATE:
      {
        testX_BEFORE_XDOT_UPDATE = true;
	testType = "Leapfrog - ModifierX";
        break;
      }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "Error - unknown action location.\n");
    }
  }
  bool testX_BEGIN_STEP;
  bool testX_BEFORE_EXPLICIT_EVAL;
  bool testX_BEFORE_X_UPDATE;
  bool testX_BEFORE_XDOT_UPDATE;
  double testX;
  double testDt;
  double testTime;
  std::string testType;
};

TEUCHOS_UNIT_TEST(LeapFrog, AppAction_ModifierX)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperLeapfrog<double>());
  stepper->setModel(model);
  auto modifierX = rcp(new StepperLeapfrogModifierXTest());
  stepper->setAppAction(modifierX);
   stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitTimeStep(15.0);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  auto icX = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icXDot = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
  auto icXDotDot = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
  auto icState = Tempus::createSolutionStateX<double>(icX, icXDot, icXDotDot);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(15.0);
  icState->setOrder   (stepper->getOrder());
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
  solutionHistory->getWorkingState()->setTimeStep(15.0);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_XDOT_UPDATE, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_X_UPDATE, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(x), 0), 1.0e-15);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-15);
  auto time = solutionHistory->getWorkingState()->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-15);
  TEST_COMPARE(modifierX->testType, ==, "Leapfrog - ModifierX");

}

// ************************************************************
// ************************************************************
class StepperLeapfrogObserverTest
  : virtual public Tempus::StepperLeapfrogObserverBase<double>
{
public:
  /// Constructor
  StepperLeapfrogObserverTest()
    : testBEGIN_STEP(false), testBEFORE_EXPLICIT_EVAL(false),
      testBEFORE_X_UPDATE(false), testBEFORE_XDOT_UPDATE(false),
      testCurrentValue(-0.99), testWorkingValue(-0.99),
      testDt(-1.5), testType("")
  {}
  /// Destructor
  virtual ~StepperLeapfrogObserverTest(){}

  /// Observe Leapfrog Stepper at action location.
  virtual void observe(
		       Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
		       Teuchos::RCP<const Tempus::StepperLeapfrog<double> > stepper,
		       const typename Tempus::StepperLeapfrogAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
    case StepperLeapfrogAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
    case StepperLeapfrogAppAction<double>::BEFORE_EXPLICIT_EVAL:
      {
        testBEFORE_EXPLICIT_EVAL = true;
        testDt = sh->getWorkingState()->getTimeStep();
        break;
      }
    case StepperLeapfrogAppAction<double>::BEFORE_X_UPDATE:
      {
        testBEFORE_X_UPDATE = true;
        testType = stepper->getStepperType();
        break;
      }
    case StepperLeapfrogAppAction<double>::BEFORE_XDOT_UPDATE:
      {
        testBEFORE_XDOT_UPDATE = true;
        auto x = sh->getWorkingState()->getX();
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
  bool testBEFORE_X_UPDATE;
  bool testBEFORE_XDOT_UPDATE;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};

TEUCHOS_UNIT_TEST(Leapfrog, AppAction_Observer)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperLeapfrog<double>());
  stepper->setModel(model);
  auto observer = rcp(new StepperLeapfrogObserverTest());
  stepper->setAppAction(observer);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  double dt = 0.173;
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  auto icX = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icXDot = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
  auto icXDotDot = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
  auto icState = Tempus::createSolutionStateX<double>(icX, icXDot, icXDotDot);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(dt);
  icState->setOrder   (stepper->getOrder());
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
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);

  // Testing that values can be observed through the observer.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(x), 0), 1.0e-15);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(x), 0), 1.0e-15);
  TEST_FLOATING_EQUALITY(observer->testDt, dt, 1.0e-15);

  TEST_COMPARE(observer->testType, ==, "Leapfrog");
}


} // namespace Tempus_Test


