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

#include "Tempus_StepperNewmarkExplicitAForm.hpp"
#include "Tempus_StepperNewmarkExplicitAFormModifierBase.hpp"
#include "Tempus_StepperNewmarkExplicitAFormModifierXBase.hpp"
#include "Tempus_StepperNewmarkExplicitAFormModifierDefault.hpp"
#include "Tempus_StepperNewmarkExplicitAFormModifierXDefault.hpp"

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
class StepperNewmarkExplicitAFormModifierTest
  : virtual public Tempus::StepperNewmarkExplicitAFormModifierBase<double> {
 public:
  /// Constructor
  StepperNewmarkExplicitAFormModifierTest()
    : testBEGIN_STEP(false),
      testBEFORE_EXPLICIT_EVAL(false),
      testAFTER_EXPLICIT_EVAL(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testDt(-1.5),
      testName("")
  {
  }

  /// Destructor
  virtual ~StepperNewmarkExplicitAFormModifierTest() {}

  /// Modify NewmarkExplicitAForm Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double>> sh,
      Teuchos::RCP<Tempus::StepperNewmarkExplicitAForm<double>> stepper,
      const typename Tempus::StepperNewmarkExplicitAFormAppAction<
          double>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperNewmarkExplicitAFormAppAction<double>::BEGIN_STEP: {
        testBEGIN_STEP = true;
        break;
      }
      case StepperNewmarkExplicitAFormAppAction<double>::BEFORE_EXPLICIT_EVAL: {
        testBEFORE_EXPLICIT_EVAL = true;
        testName                 = "Newmark Explicit A Form - Modifier";
        stepper->setStepperName(testName);
        break;
      }
      case StepperNewmarkExplicitAFormAppAction<double>::AFTER_EXPLICIT_EVAL: {
        testAFTER_EXPLICIT_EVAL = true;
        testDt                  = sh->getWorkingState()->getTimeStep();
        // sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperNewmarkExplicitAFormAppAction<double>::END_STEP: {
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
  bool testBEFORE_EXPLICIT_EVAL;
  bool testAFTER_EXPLICIT_EVAL;
  bool testEND_STEP;
  double testCurrentValue;
  double testDt;
  std::string testName;
};

// ************************************************************
// ************************************************************
class StepperNewmarkExplicitAFormModifierXTest
  : virtual public Tempus::StepperNewmarkExplicitAFormModifierXBase<double> {
 public:
  /// Constructor
  StepperNewmarkExplicitAFormModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEFORE_EXPLICIT_EVAL(false),
      testX_AFTER_EXPLICIT_EVAL(false),
      testX_END_STEP(false),
      testX(-0.99),
      testXDot(-0.99),
      testDt(-1.5),
      testTime(-1.5)
  {
  }

  /// Destructor
  virtual ~StepperNewmarkExplicitAFormModifierXTest() {}

  /// Modify NewmarkExplicitAForm Stepper at action location.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<double>> x, const double time,
      const double dt,
      const typename Tempus::StepperNewmarkExplicitAFormModifierXBase<
          double>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperNewmarkExplicitAFormModifierXBase<double>::X_BEGIN_STEP: {
        testX_BEGIN_STEP = true;
        testDt           = dt;
        testX            = get_ele(*(x), 0);
        break;
      }
      case StepperNewmarkExplicitAFormModifierXBase<
          double>::X_BEFORE_EXPLICIT_EVAL: {
        testX_BEFORE_EXPLICIT_EVAL = true;
        testTime                   = time;
        testX                      = get_ele(*(x), 0);
        break;
      }
      case StepperNewmarkExplicitAFormModifierXBase<
          double>::X_AFTER_EXPLICIT_EVAL: {
        testX_AFTER_EXPLICIT_EVAL = true;
        testX                     = get_ele(*(x), 0);
        break;
      }
      case StepperNewmarkExplicitAFormModifierXBase<double>::X_END_STEP: {
        testX_END_STEP = true;
        testTime       = time;
        testX          = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }
  }

  bool testX_BEGIN_STEP;
  bool testX_BEFORE_EXPLICIT_EVAL;
  bool testX_AFTER_EXPLICIT_EVAL;
  bool testX_END_STEP;
  double testX;
  double testXDot;
  double testDt;
  double testTime;
};

TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, AppAction_Modifier)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::sublist;

  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;

  // Read params from .xml file
  RCP<ParameterList> pList = Teuchos::getParametersFromXmlFile(
      "Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<Tempus_Test::HarmonicOscillatorModel<double>> model =
      Teuchos::rcp(new Tempus_Test::HarmonicOscillatorModel<double>(hom_pl));

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl        = sublist(pList, "Tempus", true);
  RCP<ParameterList> stepperPL = sublist(pl, "Default Stepper", true);
  stepperPL->remove("Zero Initial Guess");

  double dt = pl->sublist("Default Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  pl->sublist("Default Integrator")
      .sublist("Time Step Control")
      .set("Initial Time Step", dt);
  integrator = Tempus::createIntegratorBasic<double>(pl, model);

  RCP<Tempus::StepperNewmarkExplicitAForm<double>> stepper =
      Teuchos::rcp_dynamic_cast<Tempus::StepperNewmarkExplicitAForm<double>>(
          integrator->getStepper(), true);

  auto modifier = rcp(new StepperNewmarkExplicitAFormModifierTest());
  stepper->setAppAction(modifier);
  stepper->initialize();
  integrator->initialize();

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
  TEST_COMPARE(modifier->testBEFORE_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifier->testAFTER_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifier->testEND_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto x  = integrator->getX();
  auto Dt = integrator->getTime();
  TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-14);
  TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(x), 0), 1.0e-14);
  TEST_COMPARE(modifier->testName, ==, stepper->getStepperName());
}

TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, AppAction_ModifierX)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::sublist;

  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;

  // Read params from .xml file
  RCP<ParameterList> pList = Teuchos::getParametersFromXmlFile(
      "Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<Tempus_Test::HarmonicOscillatorModel<double>> model =
      Teuchos::rcp(new Tempus_Test::HarmonicOscillatorModel<double>(hom_pl));

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl        = sublist(pList, "Tempus", true);
  RCP<ParameterList> stepperPL = sublist(pl, "Default Stepper", true);
  stepperPL->remove("Zero Initial Guess");

  double dt = pl->sublist("Default Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  pl->sublist("Default Integrator")
      .sublist("Time Step Control")
      .set("Initial Time Step", dt);
  integrator = Tempus::createIntegratorBasic<double>(pl, model);

  RCP<Tempus::StepperNewmarkExplicitAForm<double>> stepper =
      Teuchos::rcp_dynamic_cast<Tempus::StepperNewmarkExplicitAForm<double>>(
          integrator->getStepper(), true);

  auto modifierX = rcp(new StepperNewmarkExplicitAFormModifierXTest());
  stepper->setAppAction(modifierX);
  stepper->initialize();
  integrator->initialize();

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Testing that each ACTION_LOCATION has been called.
  TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
  TEST_COMPARE(modifierX->testX_BEFORE_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifierX->testX_AFTER_EXPLICIT_EVAL, ==, true);
  TEST_COMPARE(modifierX->testX_END_STEP, ==, true);

  // Testing that values can be set through the Modifier.
  auto Dt = integrator->getTime();
  TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-14);

  const auto x = integrator->getX();
  TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(x), 0), 1.0e-14);
}

}  // namespace Tempus_Unit_Test
