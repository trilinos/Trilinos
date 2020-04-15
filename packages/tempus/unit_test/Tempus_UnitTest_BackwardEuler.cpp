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
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBackwardEulerModifierBase.hpp"
#include "Tempus_StepperBackwardEulerObserverBase.hpp"
#include "Tempus_StepperBackwardEulerModifierXBase.hpp"
#include "Tempus_StepperBackwardEulerModifierDefault.hpp"
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

using Tempus::StepperFactory;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, Construction)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Default values for construction.
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  auto predictorStepper = rcp(new Tempus::StepperForwardEuler<double>());
  predictorStepper->setModel(model);  // Can use the same model since both steppers are implicit ODEs.
  predictorStepper->initialize();

  auto defaultStepper = rcp(new Tempus::StepperBackwardEuler<double>());
  bool useFSAL              = defaultStepper->getUseFSALDefault();
  std::string ICConsistency = defaultStepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = defaultStepper->getICConsistencyCheckDefault();
  bool zeroInitialGuess     = defaultStepper->getZeroInitialGuess();

  // Test the set functions.
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  auto obs    = rcp(new Tempus::StepperBackwardEulerObserver<double>());
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
  auto modifier = rcp(new Tempus::StepperBackwardEulerModifierDefault<double>());
  stepper->setAppAction(modifier);
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setPredictor(predictorStepper);             stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  // Full argument list construction.
  stepper = rcp(new Tempus::StepperBackwardEuler<double>(
    model, obs, solver, predictorStepper, useFSAL,
    ICConsistency, ICConsistencyCheck, zeroInitialGuess));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif

  // Full argument list construction.
  stepper = rcp(new Tempus::StepperBackwardEuler<double>(
    model, solver, predictorStepper, useFSAL,
    ICConsistency, ICConsistencyCheck, zeroInitialGuess, modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("Backward Euler", model);
}


// ************************************************************
// ************************************************************
class StepperBackwardEulerModifierTest
  : virtual public Tempus::StepperBackwardEulerModifierBase<double>
{
public:

  /// Constructor
  StepperBackwardEulerModifierTest()
    : testBEGIN_STEP(false), testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false), testEND_STEP(false),
      testCurrentValue(-0.99), testWorkingValue(-0.99),
      testDt(-1.5), testType("")
  {}

  /// Destructor
  virtual ~StepperBackwardEulerModifierTest(){}

  /// Modify BackwardEuler Stepper at action location.
  virtual void modify(
    Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
    Teuchos::RCP<Tempus::StepperBackwardEuler<double> > stepper,
    const typename Tempus::StepperBackwardEulerAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperBackwardEulerAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperBackwardEulerAppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testDt = sh->getWorkingState()->getTimeStep()/10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperBackwardEulerAppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        testType = "Backward Euler - Modifier";
        stepper->setStepperType(testType);
        break;
      }
      case StepperBackwardEulerAppAction<double>::END_STEP:
      {
        testEND_STEP = true;
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
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};

TEUCHOS_UNIT_TEST(BackwardEuler, AppAction_Modifier)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> >
    model = rcp(new Tempus_Test::SinCosModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
  stepper->setModel(model);
  auto modifier = rcp(new StepperBackwardEulerModifierTest());
  stepper->setAppAction(modifier);
  stepper->initialize();

  // Create a SolutionHistory.
  auto solutionHistory = Tempus::createSolutionHistoryME(model);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 0.1;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifier->testAFTER_SOLVE, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-15);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(x), 0), 1.0e-15);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-15);

  TEST_COMPARE(modifier->testType, ==, "Backward Euler - Modifier");

}


// ************************************************************
// ************************************************************
class StepperBackwardEulerObserverTest
  : virtual public Tempus::StepperBackwardEulerObserverBase<double>
{
public:

  /// Constructor
  StepperBackwardEulerObserverTest()
    : testBEGIN_STEP(false), testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false), testEND_STEP(false),
      testCurrentValue(-0.99), testWorkingValue(-0.99),
      testDt(-1.5), testType("")
  {}

  /// Destructor
  virtual ~StepperBackwardEulerObserverTest(){}

  /// Observe BackwardEuler Stepper at end of takeStep.
  virtual void observe(
    Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
    Teuchos::RCP<const Tempus::StepperBackwardEuler<double> > stepper,
    const typename Tempus::StepperBackwardEulerAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperBackwardEulerAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperBackwardEulerAppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testDt = sh->getWorkingState()->getTimeStep();
        break;
      }
      case StepperBackwardEulerAppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        testType = stepper->getStepperType();
        break;
      }
      case StepperBackwardEulerAppAction<double>::END_STEP:
      {
        testEND_STEP = true;
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
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};

TEUCHOS_UNIT_TEST(BackwardEuler, AppAction_Observer)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> >
    model = rcp(new Tempus_Test::SinCosModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
  stepper->setModel(model);
  auto observer = rcp(new StepperBackwardEulerObserverTest());
  stepper->setAppAction(observer);
  stepper->initialize();

  // Setup a SolutionHistory.
  auto solutionHistory = Tempus::createSolutionHistoryME(model);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 0.1;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
  TEST_COMPARE(observer->testBEFORE_SOLVE, ==, true);
  TEST_COMPARE(observer->testAFTER_SOLVE, ==, true);
  TEST_COMPARE(observer->testEND_STEP, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(x), 0), 1.0e-15);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(x), 0), 1.0e-15);
  TEST_FLOATING_EQUALITY(observer->testDt, dt, 1.0e-15);

  TEST_COMPARE(observer->testType, ==, "Backward Euler");
}


// ************************************************************
// ************************************************************
class StepperBackwardEulerModifierXTest
  : virtual public Tempus::StepperBackwardEulerModifierXBase<double>
{
public:

  /// Constructor
  StepperBackwardEulerModifierXTest()
    : testX_BEGIN_STEP(false), testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false), testXDOT_END_STEP(false),
      testX(-0.99), testXDot(-0.99),
      testDt(-1.5), testTime(-1.5)
  {}

  /// Destructor
  virtual ~StepperBackwardEulerModifierXTest(){}

  /// Observe BackwardEuler Stepper at end of takeStep.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<double> > x,
    const double time, const double dt,
    const typename Tempus::StepperBackwardEulerModifierXBase<double>::MODIFIER_TYPE modType)
  {
    switch(modType) {
      case StepperBackwardEulerModifierXBase<double>::X_BEGIN_STEP:
      {
        testX_BEGIN_STEP = true;
        testX = get_ele(*(x), 0);
        break;
      }
      case StepperBackwardEulerModifierXBase<double>::X_BEFORE_SOLVE:
      {
        testX_BEFORE_SOLVE = true;
        testDt = dt;
        break;
      }
      case StepperBackwardEulerModifierXBase<double>::X_AFTER_SOLVE:
      {
        testX_AFTER_SOLVE = true;
        testTime = time;
        break;
      }
      case StepperBackwardEulerModifierXBase<double>::XDOT_END_STEP:
      {
        testXDOT_END_STEP = true;
        testXDot = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

  bool testX_BEGIN_STEP;
  bool testX_BEFORE_SOLVE;
  bool testX_AFTER_SOLVE;
  bool testXDOT_END_STEP;
  double testX;
  double testXDot;
  double testDt;
  double testTime;
};

TEUCHOS_UNIT_TEST(BackwardEuler, AppAction_ModifierX)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> >
    model = rcp(new Tempus_Test::SinCosModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
  stepper->setModel(model);
  auto modifierX = rcp(new StepperBackwardEulerModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->initialize();

  // Setup a SolutionHistory.
  auto solutionHistory = Tempus::createSolutionHistoryME(model);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 0.1;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testXDOT_END_STEP, ==, true);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(x), 0), 1.0e-15);
  // Temporary memory for xDot is not guarranteed to exist outside the Stepper.
  auto xDot = stepper->getStepperXDot(solutionHistory->getWorkingState());
  TEST_FLOATING_EQUALITY(modifierX->testXDot, get_ele(*(xDot), 0),1.0e-15);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-15);

  auto time = solutionHistory->getWorkingState()->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-15);
}


} // namespace Tempus_Test
