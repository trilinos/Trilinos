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

#include "Tempus_StepperBDF2.hpp"
#include "Tempus_StepperBDF2ModifierBase.hpp"
#include "Tempus_StepperBDF2ObserverBase.hpp"
#include "Tempus_StepperBDF2ModifierXBase.hpp"
#include "Tempus_StepperBDF2ModifierDefault.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

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
  auto solver   = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  auto startUpStepper = rcp(new Tempus::StepperDIRK_1StageTheta<double>());
  startUpStepper->setModel(
      model);  // Can use the same model since both steppers are implicit ODEs.
  startUpStepper->initialize();

  auto defaultStepper       = rcp(new Tempus::StepperBDF2<double>());
  bool useFSAL              = defaultStepper->getUseFSAL();
  std::string ICConsistency = defaultStepper->getICConsistency();
  bool ICConsistencyCheck   = defaultStepper->getICConsistencyCheck();
  bool zeroInitialGuess     = defaultStepper->getZeroInitialGuess();

  // Test the set functions.
  stepper->setAppAction(modifier);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setSolver(solver);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setStartUpStepper(startUpStepper);
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

  stepper = rcp(new Tempus::StepperBDF2<double>(
      model, solver, startUpStepper, useFSAL, ICConsistency, ICConsistencyCheck,
      zeroInitialGuess, modifier));
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
TEUCHOS_UNIT_TEST(BDF2, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      Teuchos::getParametersFromXmlFile("Tempus_BDF2_SinCos.xml");

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
        Tempus::createIntegratorBasic<double>(model, std::string("BDF2"));

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
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, ConstructingFromDefaults)
{
  double dt = 0.1;
  std::vector<std::string> options;
  options.push_back("Default Parameters");
  options.push_back("ICConsistency and Check");

  for (const auto& option : options) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        Teuchos::getParametersFromXmlFile("Tempus_BDF2_SinCos.xml");
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    // RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    auto model = rcp(new Tempus_Test::SinCosModel<double>(scm_pl));

    // Setup Stepper for field solve ----------------------------
    auto stepper = rcp(new Tempus::StepperBDF2<double>());
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
    solutionHistory->setStorageLimit(3);
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
    out << "  Stepper = " << stepper->description() << "\n            with "
        << option << std::endl;
    out << "  =========================" << std::endl;
    out << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
        << get_ele(*(x_exact), 1) << std::endl;
    out << "  Computed solution: " << get_ele(*(x), 0) << "   "
        << get_ele(*(x), 1) << std::endl;
    out << "  Difference       : " << get_ele(*(xdiff), 0) << "   "
        << get_ele(*(xdiff), 1) << std::endl;
    out << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.839732, 1.0e-4);
    TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.542663, 1.0e-4);
  }
}

// ************************************************************
// ************************************************************
class StepperBDF2ModifierTest
  : virtual public Tempus::StepperBDF2ModifierBase<double> {
 public:
  /// Constructor
  StepperBDF2ModifierTest()
    : testBEGIN_STEP(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(.99),
      testType("")
  {
  }

  /// Destructor
  virtual ~StepperBDF2ModifierTest() {}

  /// Modify BDF2 Stepper at end of takeStep.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperBDF2<double> > stepper,
      const typename Tempus::StepperBDF2AppAction<double>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperBDF2AppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getWorkingState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperBDF2AppAction<double>::BEFORE_SOLVE: {
        testBEFORE_SOLVE = true;
        testType         = "BDF2 - Modifier";
        break;
      }
      case StepperBDF2AppAction<double>::AFTER_SOLVE: {
        testAFTER_SOLVE = true;
        testDt          = sh->getCurrentState()->getTimeStep() / 10.0;
        sh->getCurrentState()->setTimeStep(testDt);
        break;
      }
      case StepperBDF2AppAction<double>::END_STEP: {
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
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(0.0);
  icState->setIndex(0);
  icState->setTimeStep(1.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime(0.0);
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
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-15);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(x), 0), 1.0e-15);
  TEST_COMPARE(modifier->testType, ==, "BDF2 - Modifier");
}

// ************************************************************
// ************************************************************
class StepperBDF2ObserverTest
  : virtual public Tempus::StepperBDF2ObserverBase<double> {
 public:
  /// Constructor
  StepperBDF2ObserverTest()
    : testBEGIN_STEP(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testEND_STEP(false),
      testCurrentValue(0.99),
      testWorkingValue(0.99),
      testDt(.99),
      testType("")
  {
  }

  /// Destructor
  virtual ~StepperBDF2ObserverTest() {}

  /// Modify BDF2 Stepper at end of takeStep.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperBDF2<double> > stepper,
      const typename Tempus::StepperBDF2AppAction<double>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperBDF2AppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP   = true;
        auto x           = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperBDF2AppAction<double>::BEFORE_SOLVE: {
        testBEFORE_SOLVE = true;
        testType         = stepper->getStepperType();
        break;
      }
      case StepperBDF2AppAction<double>::AFTER_SOLVE: {
        testAFTER_SOLVE = true;
        testDt          = sh->getCurrentState()->getTimeStep();
        break;
      }
      case StepperBDF2AppAction<double>::END_STEP: {
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
  std::string testType;
};

// ************************************************************
// ************************************************************
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
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(0.0);
  icState->setIndex(0);
  icState->setTimeStep(1.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime(0.0);
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
  TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(x), 0), 1.0e-15);
  x = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(x), 0), 1.0e-15);
  TEST_COMPARE(observer->testType, ==, "BDF2 - Modifier");
}

// ************************************************************
// ************************************************************
class StepperBDF2ModifierXTest
  : virtual public Tempus::StepperBDF2ModifierXBase<double> {
 public:
  /// Constructor
  StepperBDF2ModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false),
      testX_END_STEP(false),
      testXbegin(-.99),
      testXend(-.99),
      testTime(0.0),
      testDt(0.0)
  {
  }

  /// Destructor
  virtual ~StepperBDF2ModifierXTest() {}

  /// Modify BDF2 Stepper at end of takeStep.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<double> > x, const double time,
      const double dt,
      const typename Tempus::StepperBDF2ModifierXBase<double>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperBDF2ModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testXbegin       = get_ele(*(x), 0);
        break;
      }
      case StepperBDF2ModifierXBase<double>::X_BEFORE_SOLVE: {
        testX_BEFORE_SOLVE = true;
        testDt             = dt;
        break;
      }
      case StepperBDF2ModifierXBase<double>::X_AFTER_SOLVE: {
        testX_AFTER_SOLVE = true;
        testTime          = time;
        break;
      }
      case StepperBDF2ModifierXBase<double>::X_END_STEP: {
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
  double testTime;
  double testDt;
};

// ************************************************************
// ************************************************************
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
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime(0.0);
  icState->setIndex(0);
  icState->setTimeStep(1.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitIndex(0);
  timeStepControl->setInitTime(0.0);
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
  x       = solutionHistory->getWorkingState()->getX();
  TEST_FLOATING_EQUALITY(modifierX->testXend, get_ele(*(x), 0), 1.0e-15);
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-15);
  auto time = solutionHistory->getWorkingState()->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-15);
}

}  // namespace Tempus_Unit_Test
