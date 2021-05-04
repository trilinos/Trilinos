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

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_UnitTest_Utils.hpp"
#include "Tempus_TimeStepControl.hpp"

#include "Tempus_StepperNewmarkImplicitDForm.hpp"
#include "Tempus_StepperNewmarkImplicitDFormModifierBase.hpp"
#include "Tempus_StepperNewmarkImplicitDFormModifierXBase.hpp"
#include "Tempus_StepperNewmarkImplicitDFormModifierDefault.hpp"
#include "Tempus_StepperNewmarkImplicitDFormModifierXDefault.hpp"

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
class StepperNewmarkImplicitDFormModifierTest
  : virtual public Tempus::StepperNewmarkImplicitDFormModifierBase<double>
{
public:

  /// Constructor
  StepperNewmarkImplicitDFormModifierTest()
    : testBEGIN_STEP(false), testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false), testEND_STEP(false),
      testCurrentValue(-0.99),
      testDt(-1.5), testName("")
  {}

  /// Destructor
  virtual ~StepperNewmarkImplicitDFormModifierTest(){}

  /// Modify NewmarkImplicitDForm Stepper at action location.
  virtual void modify(
    Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
    Teuchos::RCP<Tempus::StepperNewmarkImplicitDForm<double> > stepper,
    const typename Tempus::StepperNewmarkImplicitDFormAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperNewmarkImplicitDFormAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        break;
      }
      case StepperNewmarkImplicitDFormAppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testDt = sh->getWorkingState()->getTimeStep();
        //sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperNewmarkImplicitDFormAppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        testName = "Newmark Implicit A Form - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperNewmarkImplicitDFormAppAction<double>::END_STEP:
      {
        testEND_STEP = true;
        auto x = sh->getWorkingState()->getX();
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
class StepperNewmarkImplicitDFormModifierXTest
  : virtual public Tempus::StepperNewmarkImplicitDFormModifierXBase<double>
{
public:

  /// Constructor
  StepperNewmarkImplicitDFormModifierXTest()
    : testX_BEGIN_STEP(false), testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false), testX_END_STEP(false),
      testX(-0.99), testXDot(-0.99),
      testDt(-1.5), testTime(-1.5)
  {}

  /// Destructor
  virtual ~StepperNewmarkImplicitDFormModifierXTest(){}

  /// Modify NewmarkImplicitDForm Stepper at action location.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<double> > x,
    const double time, const double dt,
    const typename Tempus::StepperNewmarkImplicitDFormModifierXBase<double>::MODIFIER_TYPE modType)
  {
    switch(modType) {
      case StepperNewmarkImplicitDFormModifierXBase<double>::X_BEGIN_STEP:
      {
        testX_BEGIN_STEP = true;
        testX = get_ele(*(x), 0);
        break;
      }
      case StepperNewmarkImplicitDFormModifierXBase<double>::X_BEFORE_SOLVE:
      {
        testX_BEFORE_SOLVE = true;
        testDt = dt;
        break;
      }
      case StepperNewmarkImplicitDFormModifierXBase<double>::X_AFTER_SOLVE:
      {
        testX_AFTER_SOLVE = true;
        testTime = time;
        testX = get_ele(*(x), 0);
        break;
      }
      case StepperNewmarkImplicitDFormModifierXBase<double>::X_END_STEP:
      {
        testX_END_STEP = true;
        testX = get_ele(*(x), 0);
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

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitDForm, Default_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperNewmarkImplicitDForm<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  auto modifier  = rcp(new Tempus::StepperNewmarkImplicitDFormModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperNewmarkImplicitDFormModifierXDefault<double>());


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
  stepper->setAppAction(modifier);                     stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(modifierX);                    stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  stepper->setSchemeName(schemeName);                  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setBeta(beta);                              stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setGamma(gamma);                            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  stepper = rcp(new Tempus::StepperNewmarkImplicitDForm<double>(
    model, solver, useFSAL,
    ICConsistency, ICConsistencyCheck, zeroInitialGuess,
    schemeName, beta, gamma, modifier));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitDForm, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());
  testFactoryConstruction("Newmark Implicit d-Form", model);
}

TEUCHOS_UNIT_TEST(NewmarkImplicitDForm, AppAction_Modifier)
{
  using Teuchos::RCP;
  using Teuchos::getParametersFromXmlFile;
  using Teuchos::sublist;
  using Teuchos::ParameterList;

  double dt = 1.0;

  // Read params from .xml file
  RCP<ParameterList> pList = getParametersFromXmlFile(
    "Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<const Thyra::ModelEvaluator<double> > model =
    Teuchos::rcp(new Tempus_Test::HarmonicOscillatorModel<double>(hom_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperNewmarkImplicitDForm<double> > stepper =
    Tempus::createStepperNewmarkImplicitDForm(model, Teuchos::null);

  auto modifier = rcp(new StepperNewmarkImplicitDFormModifierTest());
  stepper->setAppAction(modifier);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double> > timeStepControl =
    Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL = pl->sublist("Default Integrator")
                           .sublist("Time Step Control");
  timeStepControl->setInitIndex(tscPL.get<int>   ("Initial Time Index"));
  timeStepControl->setInitTime (tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(dt);
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  using Teuchos::rcp_const_cast;
  auto inArgsIC = model->getNominalValues();
  RCP<Thyra::VectorBase<double> > icX =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  RCP<Thyra::VectorBase<double> > icXDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
  RCP<Thyra::VectorBase<double> > icXDotDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
  RCP<Tempus::SolutionState<double> > icState =
    Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  RCP<Tempus::SolutionHistory<double> > solutionHistory =
    Teuchos::rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double> > integrator =
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
  auto x = integrator->getX();
  auto Dt = integrator->getTime();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-14);
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  TEST_COMPARE(modifier->testName, ==, stepper->getStepperName());
}

TEUCHOS_UNIT_TEST(NewmarkImplicitDForm, AppAction_ModifierX)
{
  using Teuchos::RCP;
  using Teuchos::getParametersFromXmlFile;
  using Teuchos::sublist;
  using Teuchos::ParameterList;

  double dt = 1.0;

  // Read params from .xml file
  RCP<ParameterList> pList = getParametersFromXmlFile(
    "Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<const Thyra::ModelEvaluator<double> > model =
    Teuchos::rcp(new Tempus_Test::HarmonicOscillatorModel<double>(hom_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperNewmarkImplicitDForm<double> > stepper =
    Tempus::createStepperNewmarkImplicitDForm(model, Teuchos::null);

  auto modifierX = rcp(new StepperNewmarkImplicitDFormModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double> > timeStepControl =
    Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL = pl->sublist("Default Integrator")
                           .sublist("Time Step Control");
  timeStepControl->setInitIndex(tscPL.get<int>   ("Initial Time Index"));
  timeStepControl->setInitTime (tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(dt);
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  using Teuchos::rcp_const_cast;
  auto inArgsIC = model->getNominalValues();
  RCP<Thyra::VectorBase<double> > icX =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  RCP<Thyra::VectorBase<double> > icXDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
  RCP<Thyra::VectorBase<double> > icXDotDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
  RCP<Tempus::SolutionState<double> > icState =
    Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  RCP<Tempus::SolutionHistory<double> > solutionHistory =
    Teuchos::rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double> > integrator =
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

} // namespace Tempus_Test
