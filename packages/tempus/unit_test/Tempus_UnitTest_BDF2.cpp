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
#include "Tempus_StepperHHTAlphaModifierBase.hpp"
#include "Tempus_StepperHHTAlphaModifierXBase.hpp"
#include "Tempus_StepperHHTAlphaObserverBase.hpp"
#include "Tempus_StepperHHTAlphaModifierDefault.hpp"
#include "Tempus_StepperHHTAlphaModifierXDefault.hpp"
#include "Tempus_StepperHHTAlphaObserverDefault.hpp"

#include "Tempus_StepperBDF2ModifierBase.hpp"
#include "Tempus_StepperBDF2ObserverBase.hpp"
#include "Tempus_StepperBDF2ModifierXBase.hpp"
#include "Tempus_StepperBDF2ModifierDefault.hpp"

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

  // Comment out any of the following tests to exclude from build/run.


  // ************************************************************
  // ************************************************************
  TEUCHOS_UNIT_TEST(BDF2, Default_Construction)
  {
    auto model = rcp(new Tempus_Test::SinCosModel<double>());

    // Default construction.
    auto stepper = rcp(new Tempus::StepperBDF2<double>());
    stepper->setModel(model);
    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


    // Default values for construction.
    auto modifier = rcp(new Tempus::StepperBDF2ModifierDefault<double>());
    auto solver = rcp(new Thyra::NOXNonlinearSolver());
    solver->setParameterList(Tempus::defaultSolverParameters());

    auto startUpStepper = rcp(new Tempus::StepperDIRK_1StageTheta<double>());
    startUpStepper->setModel(model);  // Can use the same model since both steppers are implicit ODEs.
    startUpStepper->initialize();

    auto defaultStepper = rcp(new Tempus::StepperBDF2<double>());
    bool useFSAL              = defaultStepper->getUseFSAL();
    std::string ICConsistency = defaultStepper->getICConsistency();
    bool ICConsistencyCheck   = defaultStepper->getICConsistencyCheck();
    bool zeroInitialGuess     = defaultStepper->getZeroInitialGuess();

    // Test the set functions.
    stepper->setAppAction(modifier);                     stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    stepper->setStartUpStepper(startUpStepper);          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

    stepper = rcp(new Tempus::StepperBDF2<double>(model, solver, startUpStepper, useFSAL,
						    ICConsistency, ICConsistencyCheck, zeroInitialGuess,modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    // Test stepper properties.
    TEUCHOS_ASSERT(stepper->getOrder() == 2);
  }


  // ************************************************************
  // ************************************************************
  TEUCHOS_UNIT_TEST(BDF2, StepperFactory_Construction)
  {
    auto model = rcp(new Tempus_Test::SinCosModel<double>());
    testFactoryConstruction("BDF2", model);
  }


  // ************************************************************
  // ************************************************************
class StepperBDF2ModifierTest
  : virtual public Tempus::StepperBDF2ModifierBase<double>
{
public:

  /// Constructor
  StepperBDF2ModifierTest()
    : testBEGIN_STEP(false),testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),testEND_STEP(false),
      testCurrentValue(-0.99), testWorkingValue(-0.99),
      testDt(.99), testType("")
  {}

  /// Destructor
  virtual ~StepperBDF2ModifierTest(){}

  /// Modify BDF2 Stepper at end of takeStep.
  virtual void modify(Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
		      Teuchos::RCP<Tempus::StepperBDF2<double> > stepper,
		      const typename Tempus::StepperBDF2AppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
    case StepperBDF2AppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getWorkingState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
    case StepperBDF2AppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testType = "BDF2 - Modifier";
        break;
      }
    case StepperBDF2AppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        testDt = sh->getCurrentState()->getTimeStep()/10.0;
        sh->getCurrentState()->setTimeStep(testDt);
        break;
      }
    case StepperBDF2AppAction<double>::END_STEP:
      {
        testEND_STEP = true;
	auto x  = sh->getWorkingState()->getX();
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

  TEUCHOS_UNIT_TEST(BDF2, AppAction_Modifier)
  {
    auto model = rcp(new Tempus_Test::SinCosModel<double>());

    // Setup Stepper for field solve ----------------------------
    auto stepper = rcp(new Tempus::StepperBDF2<double>());
    stepper->setModel(model);
    auto modifier = rcp(new StepperBDF2ModifierTest());
    stepper->setAppAction(modifier);
    stepper->initialize();

    // Setup initial condition SolutionState --------------------
    Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
      stepper->getModel()->getNominalValues();
  auto icSolution =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime    (0.0);
  icState->setIndex   (0);
  icState->setTimeStep(1.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setStepType ("Constant");
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime (0.0);
  timeStepControl->setFinalTime(2.0);
  timeStepControl->setInitTimeStep(1.0);
  timeStepControl->initialize();

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(3);
  solutionHistory->addState(icState);

  // Take two time steps (the first step will not test BDF2's modifier)
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 1.0;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);
  solutionHistory->promoteWorkingState();
  solutionHistory->initWorkingState();
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifier->testAFTER_SOLVE, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  auto Dt = solutionHistory->getCurrentState()->getTimeStep();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-15);
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue,  get_ele(*(x), 0), 1.0e-15);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue,  get_ele(*(x), 0), 1.0e-15);
  TEST_COMPARE(modifier->testType, ==, "BDF2 - Modifier");
  }

  // ************************************************************
  // ************************************************************
class StepperBDF2ObserverTest
  : virtual public Tempus::StepperBDF2ObserverBase<double>
{
public:

  /// Constructor
  StepperBDF2ObserverTest()
    : testBEGIN_STEP(false),testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),testEND_STEP(false),
      testCurrentValue(0.99), testWorkingValue(0.99),
      testDt(.99), testType("")
  {}

  /// Destructor
  virtual ~StepperBDF2ObserverTest(){}

  /// Modify BDF2 Stepper at end of takeStep.
  virtual void modify(Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
		      Teuchos::RCP<Tempus::StepperBDF2<double> > stepper,
		      const typename Tempus::StepperBDF2AppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
    case StepperBDF2AppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
    case StepperBDF2AppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
	testType = stepper->getStepperType();
        break;
      }
    case StepperBDF2AppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        testDt = sh->getCurrentState()->getTimeStep();
        break;
      }
    case StepperBDF2AppAction<double>::END_STEP:
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

TEUCHOS_UNIT_TEST(BDF2, AppAction_Observer)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperBDF2<double>());
  stepper->setModel(model);
  auto observer = rcp(new StepperBDF2ModifierTest());
  stepper->setAppAction(observer);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  auto icSolution =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime    (0.0);
  icState->setIndex   (0);
  icState->setTimeStep(1.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setStepType ("Constant");
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime (0.0);
  timeStepControl->setFinalTime(2.0);
  timeStepControl->setInitTimeStep(1.0);
  timeStepControl->initialize();

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(3);
  solutionHistory->addState(icState);

  // Take two time steps (the first step will not test BDF2's modifier)
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 1.0;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);
  solutionHistory->promoteWorkingState();
  solutionHistory->initWorkingState();
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
  TEST_COMPARE(observer->testBEFORE_SOLVE, ==, true);
  TEST_COMPARE(observer->testAFTER_SOLVE, ==, true);
  TEST_COMPARE(observer->testEND_STEP, ==, true);

  auto Dt = solutionHistory->getCurrentState()->getTimeStep();
  TEST_FLOATING_EQUALITY(observer->testDt, Dt, 1.0e-15);

  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(observer->testCurrentValue,  get_ele(*(x), 0), 1.0e-15);
  x      = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue,  get_ele(*(x), 0), 1.0e-15);
  TEST_COMPARE(observer->testType, ==, "BDF2 - Modifier");
}


class StepperBDF2ModifierXTest
  : virtual public Tempus::StepperBDF2ModifierXBase<double>
{
public:

  /// Constructor
  StepperBDF2ModifierXTest()
    : testX_BEGIN_STEP(false),testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false),testX_END_STEP(false),
      testXbegin(-.99),testXend(-.99),testTime(0.0),testDt(0.0)
  {}

  /// Destructor
  virtual ~StepperBDF2ModifierXTest(){}

  /// Modify BDF2 Stepper at end of takeStep.
  virtual void modify(
		      Teuchos::RCP<Thyra::VectorBase<double> > x,
		      const double time, const double dt,
		      const typename Tempus::StepperBDF2ModifierXBase<double>::MODIFIER_TYPE modType)
  {
    switch(modType) {
    case StepperBDF2ModifierXBase<double>::X_BEGIN_STEP:
      {
        testX_BEGIN_STEP = true;
        testXbegin = get_ele(*(x), 0);
        break;
      }
    case StepperBDF2ModifierXBase<double>::X_BEFORE_SOLVE:
      {
        testX_BEFORE_SOLVE = true;
        testDt = dt;
        break;
      }
    case StepperBDF2ModifierXBase<double>::X_AFTER_SOLVE:
      {
        testX_AFTER_SOLVE = true;
        testTime = time;
        break;
      }
    case StepperBDF2ModifierXBase<double>::X_END_STEP:
      {
        testX_END_STEP = true;
        testXend = get_ele(*(x), 0);
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
  double testTime;
  double testDt;
};

TEUCHOS_UNIT_TEST(BDF2, AppAction_ModifierX)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperBDF2<double>());
  stepper->setModel(model);
  auto modifierX = rcp(new StepperBDF2ModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->initialize();

  // Setup initial condition SolutionState --------------------
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  auto icSolution =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime    (0.0);
  icState->setIndex   (0);
  icState->setTimeStep(1.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setStepType ("Constant");
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime (0.0);
  timeStepControl->setFinalTime(2.0);
  timeStepControl->setInitTimeStep(1.0);
  timeStepControl->initialize();

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(3);
  solutionHistory->addState(icState);


  // Take two time steps (the first step will not test BDF2's modifier)
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 1.0;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  stepper->takeStep(solutionHistory);
  solutionHistory->promoteWorkingState();
  solutionHistory->initWorkingState();
  stepper->takeStep(solutionHistory);

  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_END_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto x = solutionHistory->getCurrentState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testXbegin, get_ele(*(x), 0), 1.0e-15);
  auto Dt = solutionHistory->getWorkingState()->getTimeStep();
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testXend, get_ele(*(x), 0), 1.0e-15);
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-15);
  auto time = solutionHistory->getWorkingState()->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-15);
  }

} // namespace Tempus_Test
