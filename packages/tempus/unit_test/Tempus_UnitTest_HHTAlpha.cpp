//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_StepperHHTAlpha.hpp"
#include "Tempus_StepperHHTAlphaModifierBase.hpp"
#include "Tempus_StepperHHTAlphaModifierXBase.hpp"
#include "Tempus_StepperHHTAlphaObserverBase.hpp"
#include "Tempus_StepperHHTAlphaModifierDefault.hpp"
#include "Tempus_StepperHHTAlphaModifierXDefault.hpp"
#include "Tempus_StepperHHTAlphaObserverDefault.hpp"

#include "../TestModels/HarmonicOscillatorModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, Default_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperHHTAlpha<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();
  std::string schemeName    = "Newmark Beta User Defined";
  double beta               = 0.25;
  double gamma              = 0.5;
  double alpha_f            = 0.0;
  double alpha_m            = 0.0;

  auto modifier  = rcp(new Tempus::StepperHHTAlphaModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperHHTAlphaModifierXDefault<double>());
  auto observer  = rcp(new Tempus::StepperHHTAlphaObserverDefault<double>());

  // Test the set functions.
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
  stepper->setZeroInitialGuess(zeroInitialGuess);
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

  stepper->setSchemeName(schemeName);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setBeta(beta);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setGamma(gamma);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAlphaF(alpha_f);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAlphaM(alpha_m);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Full argument list construction.
  stepper = rcp(new Tempus::StepperHHTAlpha<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      zeroInitialGuess, schemeName, beta, gamma, alpha_f, alpha_m, modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());
  testFactoryConstruction("HHT-Alpha", model);
}

// ************************************************************
// ************************************************************
class StepperHHTAlphaModifierTest
  : virtual public Tempus::StepperHHTAlphaModifierBase<double> {
 public:
  /// Constructor
  StepperHHTAlphaModifierTest()
    : testBEGIN_STEP(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testName("")
  {
  }

  /// Destructor
  virtual ~StepperHHTAlphaModifierTest() {}

  /// Modify HHT Alpha Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperHHTAlpha<double> > stepper,
      const typename Tempus::StepperHHTAlphaAppAction<double>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperHHTAlphaAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperHHTAlphaAppAction<double>::BEFORE_SOLVE: {
        testBEFORE_SOLVE = true;
        testDt           = sh->getWorkingState()->getTimeStep() / 10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperHHTAlphaAppAction<double>::AFTER_SOLVE: {
        testAFTER_SOLVE = true;
        testName        = "HHT Alpha - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperHHTAlphaAppAction<double>::END_STEP: {
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
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(HHTAlpha, AppAction_Modifier)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperHHTAlpha<double>());
  stepper->setModel(model);
  auto modifier = rcp(new StepperHHTAlphaModifierTest());
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

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifier->testAFTER_SOLVE, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-14);
  TEST_COMPARE(modifier->testName, ==, "HHT Alpha - Modifier");
}

// ************************************************************
// ************************************************************
class StepperHHTAlphaObserverTest
  : virtual public Tempus::StepperHHTAlphaObserverBase<double> {
 public:
  /// Constructor
  StepperHHTAlphaObserverTest()
    : testBEGIN_STEP(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testName("")
  {
  }

  /// Destructor
  virtual ~StepperHHTAlphaObserverTest() {}

  /// Observe HHTAlpha Stepper at action location.
  virtual void observe(
      Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<const Tempus::StepperHHTAlpha<double> > stepper,
      const typename Tempus::StepperHHTAlphaAppAction<double>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperHHTAlphaAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperHHTAlphaAppAction<double>::BEFORE_SOLVE: {
        testBEFORE_SOLVE = true;
        testDt           = sh->getWorkingState()->getTimeStep();
        break;
      }
      case StepperHHTAlphaAppAction<double>::AFTER_SOLVE: {
        testAFTER_SOLVE = true;
        testName        = stepper->getStepperType();
        break;
      }
      case StepperHHTAlphaAppAction<double>::END_STEP: {
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
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(HHTAlpha, AppAction_Observer)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperHHTAlpha<double>());
  stepper->setModel(model);
  auto observer = rcp(new StepperHHTAlphaObserverTest());
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

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
  TEST_COMPARE(observer->testBEFORE_SOLVE, ==, true);
  TEST_COMPARE(observer->testAFTER_SOLVE, ==, true);
  TEST_COMPARE(observer->testEND_STEP, ==, true);
  // Testing that values can be observed through the observer.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  TEST_FLOATING_EQUALITY(observer->testDt, dt, 1.0e-14);

  TEST_COMPARE(observer->testName, ==, "HHT-Alpha");
}

// ************************************************************
// ************************************************************
class StepperHHTAlphaModifierXTest
  : virtual public Tempus::StepperHHTAlphaModifierXBase<double> {
 public:
  /// Constructor
  StepperHHTAlphaModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false),
      testX_END_STEP(false),
      testXbegin(-0.99),
      testXend(-0.99),
      testDt(-1.5),
      testTime(-1.5)
  {
  }

  /// Destructor
  virtual ~StepperHHTAlphaModifierXTest() {}

  /// Modify HHTAlpha Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<double> > x, const double time,
      const double dt,
      const typename Tempus::StepperHHTAlphaModifierXBase<double>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperHHTAlphaModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testXbegin       = get_ele(*(x), 0);
        break;
      }
      case StepperHHTAlphaModifierXBase<double>::X_BEFORE_SOLVE: {
        testX_BEFORE_SOLVE = true;
        testDt             = dt;
        break;
      }
      case StepperHHTAlphaModifierXBase<double>::X_AFTER_SOLVE: {
        testX_AFTER_SOLVE = true;
        testTime          = time;
        break;
      }
      case StepperHHTAlphaModifierXBase<double>::X_END_STEP: {
        testX_END_STEP = true;
        testXend       = get_ele(*(x), 0);
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
  bool testX_END_STEP;
  double testXbegin;
  double testXend;
  double testDt;
  double testTime;
};

TEUCHOS_UNIT_TEST(HHTAlpha, AppAction_ModifierX)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperHHTAlpha<double>());
  stepper->setModel(model);
  auto modifierX = rcp(new StepperHHTAlphaModifierXTest());
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

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_END_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto xbegin = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testXbegin, get_ele(*(xbegin), 0), 1.0e-14);
  auto xend = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testXend, get_ele(*(xend), 0), 1.0e-14);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-14);
  auto time = solutionHistory->getWorkingState()->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-14);
}

}  // namespace Tempus_Unit_Test
