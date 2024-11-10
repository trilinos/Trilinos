//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tempus_StepperRKButcherTableau.hpp"

#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"

#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_StepperOperatorSplitModifierBase.hpp"
#include "Tempus_StepperOperatorSplitModifierXBase.hpp"
#include "Tempus_StepperOperatorSplitObserverBase.hpp"
#include "Tempus_StepperOperatorSplitModifierDefault.hpp"
#include "Tempus_StepperOperatorSplitModifierXDefault.hpp"
#include "Tempus_StepperOperatorSplitObserverDefault.hpp"

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"

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
TEUCHOS_UNIT_TEST(OperatorSplit, Default_Construction)
{
  RCP<const Thyra::ModelEvaluator<double> > explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  RCP<const Thyra::ModelEvaluator<double> > implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperOperatorSplit<double>());
  auto subStepper1 =
      Tempus::createStepperForwardEuler(explicitModel, Teuchos::null);
  auto subStepper2 =
      Tempus::createStepperBackwardEuler(implicitModel, Teuchos::null);
  stepper->addStepper(subStepper1);
  stepper->addStepper(subStepper2);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier =
      rcp(new Tempus::StepperOperatorSplitModifierDefault<double>());
  auto modifierX =
      rcp(new Tempus::StepperOperatorSplitModifierXDefault<double>());
  auto observer =
      rcp(new Tempus::StepperOperatorSplitObserverDefault<double>());
  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();
  int order                 = 1;

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
  stepper->setOrder(order);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setOrderMin(order);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setOrderMax(order);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Full argument list construction.
  std::vector<RCP<const Thyra::ModelEvaluator<double> > > models;
  models.push_back(explicitModel);
  models.push_back(implicitModel);

  std::vector<Teuchos::RCP<Tempus::Stepper<double> > > subStepperList;
  subStepperList.push_back(subStepper1);
  subStepperList.push_back(subStepper2);

  stepper = rcp(new Tempus::StepperOperatorSplit<double>(
      models, subStepperList, useFSAL, ICConsistency, ICConsistencyCheck, order,
      order, order, modifier));

  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(OperatorSplit, StepperFactory_Construction)
{
  // Read params from .xml file
  auto pList = Teuchos::getParametersFromXmlFile(
      "../test/OperatorSplit/Tempus_OperatorSplit_VanDerPol.xml");
  auto tempusPL  = sublist(pList, "Tempus", true);
  auto stepperPL = sublist(tempusPL, "Demo Stepper", true);

  auto explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  std::vector<RCP<const Thyra::ModelEvaluator<double> > > models;
  models.push_back(explicitModel);
  models.push_back(implicitModel);

  auto sf = Teuchos::rcp(new Tempus::StepperFactory<double>());

  // Test using ParameterList.
  // Passing in model.
  auto stepper = sf->createStepper(stepperPL, models);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}

// ************************************************************
// ************************************************************
class StepperOperatorSplitModifierTest
  : virtual public Tempus::StepperOperatorSplitModifierBase<double> {
 public:
  /// Constructor
  StepperOperatorSplitModifierTest()
    : testBEGIN_STEP(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testName("")
  {
  }

  /// Destructor
  virtual ~StepperOperatorSplitModifierTest() {}

  /// Modify OperatorSplit Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperOperatorSplit<double> > stepper,
      const typename Tempus::StepperOperatorSplitAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperOperatorSplitAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperOperatorSplitAppAction<double>::BEFORE_STEPPER: {
        testBEFORE_STEPPER = true;
        testDt             = sh->getWorkingState()->getTimeStep() / 10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperOperatorSplitAppAction<double>::AFTER_STEPPER: {
        testAFTER_STEPPER = true;
        testName          = "OperatorSplit - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperOperatorSplitAppAction<double>::END_STEP: {
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
  bool testBEFORE_STEPPER;
  bool testAFTER_STEPPER;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(OperatorSplit, AppAction_Modifier)
{
  RCP<const Thyra::ModelEvaluator<double> > explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  RCP<const Thyra::ModelEvaluator<double> > implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  // Default construction.
  auto stepper = rcp(new Tempus::StepperOperatorSplit<double>());
  auto subStepper1 =
      Tempus::createStepperForwardEuler(explicitModel, Teuchos::null);
  auto subStepper2 =
      Tempus::createStepperBackwardEuler(implicitModel, Teuchos::null);
  auto modifier = rcp(new StepperOperatorSplitModifierTest());
  stepper->setAppAction(modifier);
  stepper->addStepper(subStepper1);
  stepper->addStepper(subStepper2);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = stepper->getModel()->getNominalValues();
  auto icX      = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icXDot =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x_dot());
  auto icState = Tempus::createSolutionStateX(icX, icXDot);
  icState->setTime(0.0);
  icState->setIndex(1);
  icState->setTimeStep(-15.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Create a SolutionHistory.
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  solutionHistory->getWorkingState()->setTimeStep(-15.0);
  stepper->takeStep(solutionHistory);

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_STEPPER, ==, true);
  TEST_COMPARE(modifier->testAFTER_STEPPER, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-14);

  TEST_COMPARE(modifier->testName, ==, "OperatorSplit - Modifier");
}

// ************************************************************
// ************************************************************
class StepperOperatorSplitObserverTest
  : virtual public Tempus::StepperOperatorSplitObserverBase<double> {
 public:
  /// Constructor
  StepperOperatorSplitObserverTest()
    : testBEGIN_STEP(false),
      testBEFORE_STEPPER(false),
      testAFTER_STEPPER(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testName("Operator Split")
  {
  }

  /// Destructor
  virtual ~StepperOperatorSplitObserverTest() {}

  /// Observe OperatorSplit Stepper at action location.
  virtual void observe(
      Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<const Tempus::StepperOperatorSplit<double> > stepper,
      const typename Tempus::StepperOperatorSplitAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperOperatorSplitAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperOperatorSplitAppAction<double>::BEFORE_STEPPER: {
        testBEFORE_STEPPER = true;
        testDt             = sh->getWorkingState()->getTimeStep();
        break;
      }
      case StepperOperatorSplitAppAction<double>::AFTER_STEPPER: {
        testAFTER_STEPPER = true;
        testName          = stepper->getStepperType();
        break;
      }
      case StepperOperatorSplitAppAction<double>::END_STEP: {
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
  bool testBEFORE_STEPPER;
  bool testAFTER_STEPPER;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testName;
};

TEUCHOS_UNIT_TEST(OperatorSplit, AppAction_Observer)
{
  RCP<const Thyra::ModelEvaluator<double> > explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  RCP<const Thyra::ModelEvaluator<double> > implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  // Default construction.
  auto stepper = rcp(new Tempus::StepperOperatorSplit<double>());
  auto subStepper1 =
      Tempus::createStepperForwardEuler(explicitModel, Teuchos::null);
  auto subStepper2 =
      Tempus::createStepperBackwardEuler(implicitModel, Teuchos::null);
  auto observer = rcp(new StepperOperatorSplitObserverTest());
  stepper->setAppAction(observer);
  stepper->addStepper(subStepper1);
  stepper->addStepper(subStepper2);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = stepper->getModel()->getNominalValues();
  auto icX      = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icXDot =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x_dot());
  auto icState = Tempus::createSolutionStateX(icX, icXDot);
  icState->setTime(0.0);
  icState->setIndex(1);
  icState->setTimeStep(-1.5);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Create a SolutionHistory.
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  solutionHistory->getWorkingState()->setTimeStep(-1.5);
  stepper->takeStep(solutionHistory);

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
  TEST_COMPARE(observer->testBEFORE_STEPPER, ==, true);
  TEST_COMPARE(observer->testAFTER_STEPPER, ==, true);
  TEST_COMPARE(observer->testEND_STEP, ==, true);

  // Testing that values can be observed through the observer.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(x), 0), 1.0e-14);
  TEST_FLOATING_EQUALITY(observer->testDt, -1.5, 1.0e-14);

  TEST_COMPARE(observer->testName, ==, "Operator Split");
}

// ************************************************************
// ************************************************************
class StepperOperatorSplitModifierXTest
  : virtual public Tempus::StepperOperatorSplitModifierXBase<double> {
 public:
  /// Constructor
  StepperOperatorSplitModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEFORE_STEPPER(false),
      testX_AFTER_STEPPER(false),
      testXDOT_END_STEP(false),
      testX(-0.99),
      testXDot(-0.99),
      testDt(-1.5),
      testTime(-1.5)
  {
  }

  /// Destructor
  virtual ~StepperOperatorSplitModifierXTest() {}

  /// Modify OperatorSplit Stepper at action location.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<double> > x,
                      const double time, const double dt,
                      const typename Tempus::StepperOperatorSplitModifierXBase<
                          double>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperOperatorSplitModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testX            = get_ele(*(x), 0);
        break;
      }
      case StepperOperatorSplitModifierXBase<double>::X_BEFORE_STEPPER: {
        testX_BEFORE_STEPPER = true;
        testDt               = dt;
        break;
      }
      case StepperOperatorSplitModifierXBase<double>::X_AFTER_STEPPER: {
        testX_AFTER_STEPPER = true;
        testTime            = time;
        break;
      }
      case StepperOperatorSplitModifierXBase<double>::XDOT_END_STEP: {
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
  bool testX_BEFORE_STEPPER;
  bool testX_AFTER_STEPPER;
  bool testXDOT_END_STEP;
  double testX;
  double testXDot;
  double testDt;
  double testTime;
};

TEUCHOS_UNIT_TEST(OperatorSplit, AppAction_ModifierX)
{
  RCP<const Thyra::ModelEvaluator<double> > explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  RCP<const Thyra::ModelEvaluator<double> > implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  // Default construction.
  auto stepper = rcp(new Tempus::StepperOperatorSplit<double>());
  auto subStepper1 =
      Tempus::createStepperForwardEuler(explicitModel, Teuchos::null);
  auto subStepper2 =
      Tempus::createStepperBackwardEuler(implicitModel, Teuchos::null);
  auto modifierX = rcp(new StepperOperatorSplitModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->addStepper(subStepper1);
  stepper->addStepper(subStepper2);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = stepper->getModel()->getNominalValues();
  auto icX      = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icXDot =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x_dot());
  auto icState = Tempus::createSolutionStateX(icX, icXDot);
  icState->setTime(0.0);
  icState->setIndex(1);
  icState->setTimeStep(-1.5);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Create a SolutionHistory.
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  solutionHistory->getWorkingState()->setTimeStep(-1.5);
  stepper->takeStep(solutionHistory);

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_STEPPER, ==, true);
  TEST_COMPARE(modifierX->testX_AFTER_STEPPER, ==, true);
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
