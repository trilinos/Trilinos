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

#include "Tempus_TimeStepControl.hpp"

#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitAFormModifierBase.hpp"
#include "Tempus_StepperNewmarkImplicitAFormModifierXBase.hpp"
#include "Tempus_StepperNewmarkImplicitAFormModifierDefault.hpp"
#include "Tempus_StepperNewmarkImplicitAFormModifierXDefault.hpp"

#include "../TestModels/HarmonicOscillatorModel.hpp"

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
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, Default_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperNewmarkImplicitAForm<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  auto modifier =
      rcp(new Tempus::StepperNewmarkImplicitAFormModifierDefault<double>());
  auto modifierX =
      rcp(new Tempus::StepperNewmarkImplicitAFormModifierXDefault<double>());

  // Default values for construction.
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();
  std::string schemeName    = "Average Acceleration";
  double beta               = 0.25;
  double gamma              = 0.5;

  // Test the set functions.
  stepper->setAppAction(modifier);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(modifierX);
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
  stepper->setZeroInitialGuess(zeroInitialGuess);
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

  // Full argument list construction.
  stepper = rcp(new Tempus::StepperNewmarkImplicitAForm<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      zeroInitialGuess, schemeName, beta, gamma, modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());
  testFactoryConstruction("Newmark Implicit a-Form", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, ConstructingFromDefaults)
{
  double dt = 1.0;
  std::vector<std::string> options;
  options.push_back("Default Parameters");
  options.push_back("ICConsistency and Check");

  for (const auto& option : options) {
    // Read params from .xml file
    RCP<ParameterList> pList = Teuchos::getParametersFromXmlFile(
        "Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder."
        "xml");
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);

    // Setup the HarmonicOscillatorModel
    RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
    auto model =
        Teuchos::rcp(new Tempus_Test::HarmonicOscillatorModel<double>(hom_pl));
    auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double>>(model);

    // Setup Stepper for field solve ----------------------------
    RCP<Tempus::StepperNewmarkImplicitAForm<double>> stepper =
        Tempus::createStepperNewmarkImplicitAForm(modelME, Teuchos::null);
    if (option == "ICConsistency and Check") {
      stepper->setICConsistency("Consistent");
      stepper->setICConsistencyCheck(true);
    }
    stepper->initialize();

    // Setup TimeStepControl ------------------------------------
    RCP<Tempus::TimeStepControl<double>> timeStepControl =
        Teuchos::rcp(new Tempus::TimeStepControl<double>());
    ParameterList tscPL =
        pl->sublist("Default Integrator").sublist("Time Step Control");
    timeStepControl->setInitIndex(tscPL.get<int>("Initial Time Index"));
    timeStepControl->setInitTime(tscPL.get<double>("Initial Time"));
    timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
    timeStepControl->setInitTimeStep(dt);
    timeStepControl->initialize();

    // Setup initial condition SolutionState --------------------
    using Teuchos::rcp_const_cast;
    auto inArgsIC = model->getNominalValues();
    RCP<Thyra::VectorBase<double>> icX =
        rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
    RCP<Thyra::VectorBase<double>> icXDot =
        rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot());
    RCP<Thyra::VectorBase<double>> icXDotDot =
        rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot_dot());
    RCP<Tempus::SolutionState<double>> icState =
        Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
    icState->setTime(timeStepControl->getInitTime());
    icState->setIndex(timeStepControl->getInitIndex());
    icState->setTimeStep(0.0);
    icState->setOrder(stepper->getOrder());
    icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

    // Setup SolutionHistory ------------------------------------
    RCP<Tempus::SolutionHistory<double>> solutionHistory =
        Teuchos::rcp(new Tempus::SolutionHistory<double>());
    solutionHistory->setName("Forward States");
    solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
    solutionHistory->setStorageLimit(2);
    solutionHistory->addState(icState);

    // Ensure ICs are consistent.
    stepper->setInitialConditions(solutionHistory);

    // Setup Integrator -----------------------------------------
    RCP<Tempus::IntegratorBasic<double>> integrator =
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
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Calculate the error
    RCP<Thyra::VectorBase<double>> xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

    // Check the order and intercept
    out << "  Stepper = " << stepper->description() << "\n            with "
        << option << std::endl;
    out << "  =========================" << std::endl;
    out << "  Exact solution   : " << get_ele(*(x_exact), 0) << std::endl;
    out << "  Computed solution: " << get_ele(*(x), 0) << std::endl;
    out << "  Difference       : " << get_ele(*(xdiff), 0) << std::endl;
    out << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(get_ele(*(x), 0), -0.222222, 1.0e-4);
  }
}

// ************************************************************
// ************************************************************
class StepperNewmarkImplicitAFormModifierTest
  : virtual public Tempus::StepperNewmarkImplicitAFormModifierBase<double> {
 public:
  /// Constructor
  StepperNewmarkImplicitAFormModifierTest()
    : testBEGIN_STEP(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testDt(-1.5),
      testName("")
  {
  }

  /// Destructor
  virtual ~StepperNewmarkImplicitAFormModifierTest() {}

  /// Modify NewmarkImplicitAForm Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double>> sh,
      Teuchos::RCP<Tempus::StepperNewmarkImplicitAForm<double>> stepper,
      const typename Tempus::StepperNewmarkImplicitAFormAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperNewmarkImplicitAFormAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP = true;
        break;
      }
      case StepperNewmarkImplicitAFormAppAction<double>::BEFORE_SOLVE: {
        testBEFORE_SOLVE = true;
        testDt           = sh->getWorkingState()->getTimeStep();
        // sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperNewmarkImplicitAFormAppAction<double>::AFTER_SOLVE: {
        testAFTER_SOLVE = true;
        testName        = "Newmark Implicit A Form - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperNewmarkImplicitAFormAppAction<double>::END_STEP: {
        testEND_STEP     = true;
        auto x           = sh->getWorkingState()->getX();
        testCurrentValue = get_ele(*(x), 0);
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
  double testDt;
  std::string testName;
};

// ************************************************************
// ************************************************************
class StepperNewmarkImplicitAFormModifierXTest
  : virtual public Tempus::StepperNewmarkImplicitAFormModifierXBase<double> {
 public:
  /// Constructor
  StepperNewmarkImplicitAFormModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false),
      testX_END_STEP(false),
      testX(-0.99),
      testXDot(-0.99),
      testDt(-1.5),
      testTime(-1.5)
  {
  }

  /// Destructor
  virtual ~StepperNewmarkImplicitAFormModifierXTest() {}

  /// Modify NewmarkImplicitAForm Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<double>> x, const double time,
      const double dt,
      const typename Tempus::StepperNewmarkImplicitAFormModifierXBase<
          double>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperNewmarkImplicitAFormModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testX            = get_ele(*(x), 0);
        break;
      }
      case StepperNewmarkImplicitAFormModifierXBase<double>::X_BEFORE_SOLVE: {
        testX_BEFORE_SOLVE = true;
        testDt             = dt;
        break;
      }
      case StepperNewmarkImplicitAFormModifierXBase<double>::X_AFTER_SOLVE: {
        testX_AFTER_SOLVE = true;
        testTime          = time;
        testX             = get_ele(*(x), 0);
        break;
      }
      case StepperNewmarkImplicitAFormModifierXBase<double>::X_END_STEP: {
        testX_END_STEP = true;
        testX          = get_ele(*(x), 0);
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
  double testX;
  double testXDot;
  double testDt;
  double testTime;
};

TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, AppAction_Modifier)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::sublist;

  double dt = 1.0;

  // Read params from .xml file
  RCP<ParameterList> pList = Teuchos::getParametersFromXmlFile(
      "Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<const Thyra::ModelEvaluator<double>> model =
      Teuchos::rcp(new Tempus_Test::HarmonicOscillatorModel<double>(hom_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperNewmarkImplicitAForm<double>> stepper =
      Tempus::createStepperNewmarkImplicitAForm(model, Teuchos::null);

  auto modifier = rcp(new StepperNewmarkImplicitAFormModifierTest());
  stepper->setAppAction(modifier);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double>> timeStepControl =
      Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL =
      pl->sublist("Default Integrator").sublist("Time Step Control");
  timeStepControl->setInitIndex(tscPL.get<int>("Initial Time Index"));
  timeStepControl->setInitTime(tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(dt);
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  using Teuchos::rcp_const_cast;
  auto inArgsIC = model->getNominalValues();
  RCP<Thyra::VectorBase<double>> icX =
      rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  RCP<Thyra::VectorBase<double>> icXDot =
      rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot());
  RCP<Thyra::VectorBase<double>> icXDotDot =
      rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot_dot());
  RCP<Tempus::SolutionState<double>> icState =
      Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
  icState->setTime(timeStepControl->getInitTime());
  icState->setIndex(timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  RCP<Tempus::SolutionHistory<double>> solutionHistory =
      Teuchos::rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double>> integrator =
      Tempus::createIntegratorBasic<double>();
  integrator->setStepper(stepper);
  integrator->setTimeStepControl(timeStepControl);
  integrator->setSolutionHistory(solutionHistory);
  integrator->initialize();

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifier->testAFTER_SOLVE, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto x  = integrator->getX();
  auto Dt = integrator->getTime();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-14);
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  TEST_COMPARE(modifier->testName, ==, stepper->getStepperName());
}

TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, AppAction_ModifierX)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::sublist;

  double dt = 1.0;

  // Read params from .xml file
  RCP<ParameterList> pList = Teuchos::getParametersFromXmlFile(
      "Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<const Thyra::ModelEvaluator<double>> model =
      Teuchos::rcp(new Tempus_Test::HarmonicOscillatorModel<double>(hom_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperNewmarkImplicitAForm<double>> stepper =
      Tempus::createStepperNewmarkImplicitAForm(model, Teuchos::null);

  auto modifierX = rcp(new StepperNewmarkImplicitAFormModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double>> timeStepControl =
      Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL =
      pl->sublist("Default Integrator").sublist("Time Step Control");
  timeStepControl->setInitIndex(tscPL.get<int>("Initial Time Index"));
  timeStepControl->setInitTime(tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(dt);
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  using Teuchos::rcp_const_cast;
  auto inArgsIC = model->getNominalValues();
  RCP<Thyra::VectorBase<double>> icX =
      rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  RCP<Thyra::VectorBase<double>> icXDot =
      rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot());
  RCP<Thyra::VectorBase<double>> icXDotDot =
      rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot_dot());
  RCP<Tempus::SolutionState<double>> icState =
      Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
  icState->setTime(timeStepControl->getInitTime());
  icState->setIndex(timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  RCP<Tempus::SolutionHistory<double>> solutionHistory =
      Teuchos::rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double>> integrator =
      Tempus::createIntegratorBasic<double>();
  integrator->setStepper(stepper);
  integrator->setTimeStepControl(timeStepControl);
  integrator->setSolutionHistory(solutionHistory);
  integrator->initialize();

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
  TEST_COMPARE(modifierX->testX_END_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto Dt = integrator->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-14);

  const auto x = integrator->getX();
  TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(x), 0), 1.0e-14);
}

}  // namespace Tempus_Unit_Test
