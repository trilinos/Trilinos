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

#include "Tempus_StepperBackwardEulerModifierBase.hpp"
#include "Tempus_StepperBackwardEulerModifierXBase.hpp"
#include "Tempus_StepperBackwardEulerObserverBase.hpp"
#include "Tempus_StepperBackwardEulerModifierDefault.hpp"
#include "Tempus_StepperBackwardEulerModifierXDefault.hpp"
#include "Tempus_StepperBackwardEulerObserverDefault.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, Default_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier =
      rcp(new Tempus::StepperBackwardEulerModifierDefault<double>());
  auto modifierX =
      rcp(new Tempus::StepperBackwardEulerModifierXDefault<double>());
  auto observer =
      rcp(new Tempus::StepperBackwardEulerObserverDefault<double>());
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  auto predictorStepper = rcp(new Tempus::StepperForwardEuler<double>());
  predictorStepper->setModel(
      model);  // Can use the same model since both steppers are implicit ODEs.
  predictorStepper->initialize();

  auto defaultStepper       = rcp(new Tempus::StepperBackwardEuler<double>());
  bool useFSAL              = defaultStepper->getUseFSAL();
  std::string ICConsistency = defaultStepper->getICConsistency();
  bool ICConsistencyCheck   = defaultStepper->getICConsistencyCheck();
  bool zeroInitialGuess     = defaultStepper->getZeroInitialGuess();

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
  stepper->setPredictor(predictorStepper);
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

  // Full argument list construction.
  stepper = rcp(new Tempus::StepperBackwardEuler<double>(
      model, solver, predictorStepper, useFSAL, ICConsistency,
      ICConsistencyCheck, zeroInitialGuess, modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
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
  : virtual public Tempus::StepperBackwardEulerModifierBase<double> {
 public:
  /// Constructor
  StepperBackwardEulerModifierTest()
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
  virtual ~StepperBackwardEulerModifierTest() {}

  /// Modify BackwardEuler Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperBackwardEuler<double> > stepper,
      const typename Tempus::StepperBackwardEulerAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperBackwardEulerAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperBackwardEulerAppAction<double>::BEFORE_SOLVE: {
        testBEFORE_SOLVE = true;
        testDt           = sh->getWorkingState()->getTimeStep() / 10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperBackwardEulerAppAction<double>::AFTER_SOLVE: {
        testAFTER_SOLVE = true;
        testName        = "Backward Euler - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperBackwardEulerAppAction<double>::END_STEP: {
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

TEUCHOS_UNIT_TEST(BackwardEuler, AppAction_Modifier)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::SinCosModel<double>());

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

  TEST_COMPARE(modifier->testName, ==, "Backward Euler - Modifier");
}

// ************************************************************
// ************************************************************
class StepperBackwardEulerObserverTest
  : virtual public Tempus::StepperBackwardEulerObserverBase<double> {
 public:
  /// Constructor
  StepperBackwardEulerObserverTest()
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
  virtual ~StepperBackwardEulerObserverTest() {}

  /// Observe BackwardEuler Stepper at action location.
  virtual void observe(
      Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<const Tempus::StepperBackwardEuler<double> > stepper,
      const typename Tempus::StepperBackwardEulerAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperBackwardEulerAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperBackwardEulerAppAction<double>::BEFORE_SOLVE: {
        testBEFORE_SOLVE = true;
        testDt           = sh->getWorkingState()->getTimeStep();
        break;
      }
      case StepperBackwardEulerAppAction<double>::AFTER_SOLVE: {
        testAFTER_SOLVE = true;
        testName        = stepper->getStepperType();
        break;
      }
      case StepperBackwardEulerAppAction<double>::END_STEP: {
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

TEUCHOS_UNIT_TEST(BackwardEuler, AppAction_Observer)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::SinCosModel<double>());

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

  TEST_COMPARE(observer->testName, ==, "Backward Euler");
}

// ************************************************************
// ************************************************************
class StepperBackwardEulerModifierXTest
  : virtual public Tempus::StepperBackwardEulerModifierXBase<double> {
 public:
  /// Constructor
  StepperBackwardEulerModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false),
      testXDOT_END_STEP(false),
      testX(-0.99),
      testXDot(-0.99),
      testDt(-1.5),
      testTime(-1.5)
  {
  }

  /// Destructor
  virtual ~StepperBackwardEulerModifierXTest() {}

  /// Modify BackwardEuler Stepper at action location.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<double> > x,
                      const double time, const double dt,
                      const typename Tempus::StepperBackwardEulerModifierXBase<
                          double>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperBackwardEulerModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testX            = get_ele(*(x), 0);
        break;
      }
      case StepperBackwardEulerModifierXBase<double>::X_BEFORE_SOLVE: {
        testX_BEFORE_SOLVE = true;
        testDt             = dt;
        break;
      }
      case StepperBackwardEulerModifierXBase<double>::X_AFTER_SOLVE: {
        testX_AFTER_SOLVE = true;
        testTime          = time;
        break;
      }
      case StepperBackwardEulerModifierXBase<double>::XDOT_END_STEP: {
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
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::SinCosModel<double>());

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

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
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

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      Teuchos::getParametersFromXmlFile("Tempus_BackwardEuler_SinCos.xml");

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  auto model                = rcp(new Tempus_Test::SinCosModel<double>(scm_pl));

  RCP<ParameterList> tempusPL = sublist(pList, "Tempus", true);

  // Test constructor IntegratorBasic(tempusPL, model)
  {
    RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::createIntegratorBasic<double>(tempusPL, model);

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
    RCP<const ParameterList> defaultPL =
        integrator->getStepper()->getValidParameters();

    bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
    if (!pass) {
      out << std::endl;
      out << "stepperPL -------------- \n"
          << *stepperPL << std::endl;
      out << "defaultPL -------------- \n"
          << *defaultPL << std::endl;
    }
    TEST_ASSERT(pass)
  }

  // Test constructor IntegratorBasic(model, stepperType)
  {
    RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::createIntegratorBasic<double>(model,
                                              std::string("Backward Euler"));

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
    // Match Predictor for comparison
    stepperPL->set("Predictor Stepper Type", "None");
    RCP<const ParameterList> defaultPL =
        integrator->getStepper()->getValidParameters();

    bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
    if (!pass) {
      out << std::endl;
      out << "stepperPL -------------- \n"
          << *stepperPL << std::endl;
      out << "defaultPL -------------- \n"
          << *defaultPL << std::endl;
    }
    TEST_ASSERT(pass)
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ConstructingFromDefaults)
{
  double dt = 0.1;
  std::vector<std::string> options;
  options.push_back("Default Parameters");
  options.push_back("ICConsistency and Check");

  for (const auto& option : options) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        Teuchos::getParametersFromXmlFile("Tempus_BackwardEuler_SinCos.xml");
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    // RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    auto model = rcp(new Tempus_Test::SinCosModel<double>(scm_pl));

    // Setup Stepper for field solve ----------------------------
    auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
    stepper->setModel(model);
    if (option == "ICConsistency and Check") {
      stepper->setICConsistency("Consistent");
      stepper->setICConsistencyCheck(true);
    }
    stepper->initialize();

    // Setup TimeStepControl ------------------------------------
    auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
    ParameterList tscPL =
        pl->sublist("Default Integrator").sublist("Time Step Control");
    timeStepControl->setInitIndex(tscPL.get<int>("Initial Time Index"));
    timeStepControl->setInitTime(tscPL.get<double>("Initial Time"));
    timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
    timeStepControl->setInitTimeStep(dt);
    timeStepControl->initialize();

    // Setup initial condition SolutionState --------------------
    auto inArgsIC = model->getNominalValues();
    auto icSolution =
        rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
    auto icState = Tempus::createSolutionStateX(icSolution);
    icState->setTime(timeStepControl->getInitTime());
    icState->setIndex(timeStepControl->getInitIndex());
    icState->setTimeStep(0.0);
    icState->setOrder(stepper->getOrder());
    icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

    // Setup SolutionHistory ------------------------------------
    auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
    solutionHistory->setName("Forward States");
    solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
    solutionHistory->setStorageLimit(2);
    solutionHistory->addState(icState);

    // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
    stepper->setInitialConditions(solutionHistory);

    // Setup Integrator -----------------------------------------
    RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::createIntegratorBasic<double>();
    integrator->setStepper(stepper);
    integrator->setTimeStepControl(timeStepControl);
    integrator->setSolutionHistory(solutionHistory);
    // integrator->setObserver(...);
    integrator->initialize();

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time      = integrator->getTime();
    double timeFinal = pl->sublist("Default Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
        model->getExactSolution(time).get_x();

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

    // Check the order and intercept
    out << "  Stepper = " << stepper->description() << " with " << option
        << std::endl;
    out << "  =========================" << std::endl;
    out << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
        << get_ele(*(x_exact), 1) << std::endl;
    out << "  Computed solution: " << get_ele(*(x), 0) << "   "
        << get_ele(*(x), 1) << std::endl;
    out << "  Difference       : " << get_ele(*(xdiff), 0) << "   "
        << get_ele(*(xdiff), 1) << std::endl;
    out << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.798923, 1.0e-4);
    TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.516729, 1.0e-4);
  }
}

}  // namespace Tempus_Unit_Test
