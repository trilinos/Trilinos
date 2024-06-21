//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"

#include "Tempus_StepperRKModifierBase.hpp"

#include "../TestModels/DahlquistTestModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_TrapezoidalRule, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperEDIRK_TrapezoidalRule<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_TrapezoidalRule, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Trapezoidal Rule", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EDIRK_TrapezoidalRule, AppAction)
{
  auto stepper = rcp(new Tempus::StepperEDIRK_TrapezoidalRule<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

// ************************************************************
// ************************************************************

class StepperRKModifierEDIRK_TrapezoidaTest
  : virtual public Tempus::StepperRKModifierBase<double> {
 public:
  /// Constructor
  StepperRKModifierEDIRK_TrapezoidaTest(Teuchos::FancyOStream &Out,
                                        bool &Success)
    : out(Out), success(Success)
  {
  }

  // FSAL
  /// Destructor
  virtual ~StepperRKModifierEDIRK_TrapezoidaTest() {}

  /// Test the modify RK Stepper at the Action Locations.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
      Teuchos::RCP<Tempus::StepperRKBase<double> > stepper,
      const typename Tempus::StepperRKAppAction<double>::ACTION_LOCATION actLoc)
  {
    const double relTol                       = 1.0e-14;
    auto stageNumber                          = stepper->getStageNumber();
    Teuchos::SerialDenseVector<int, double> c = stepper->getTableau()->c();

    auto currentState = sh->getCurrentState();
    auto workingState = sh->getWorkingState();
    const double dt   = workingState->getTimeStep();
    double time       = currentState->getTime();
    if (stageNumber >= 0) time += c(stageNumber) * dt;

    auto x    = workingState->getX();
    auto xDot = workingState->getXDot();
    if (xDot == Teuchos::null) xDot = stepper->getStepperXDot();

    switch (actLoc) {
      case StepperRKAppAction<double>::BEGIN_STEP: {
        {
          auto DME = Teuchos::rcp_dynamic_cast<
              const Tempus_Test::DahlquistTestModel<double> >(
              stepper->getModel());
          TEST_FLOATING_EQUALITY(DME->getLambda(), -1.0, relTol);
        }
        TEST_FLOATING_EQUALITY(dt, 1.0, relTol);

        const double x_0    = get_ele(*(x), 0);
        const double xDot_0 = get_ele(*(xDot), 0);
        TEST_FLOATING_EQUALITY(x_0, 1.0, relTol);      // Should be x_0
        TEST_FLOATING_EQUALITY(xDot_0, -1.0, relTol);  // Should be xDot_0
        TEST_ASSERT(std::abs(time) < relTol);
        TEST_FLOATING_EQUALITY(dt, 1.0, relTol);
        TEST_COMPARE(stageNumber, ==, -1);
        break;
      }
      case StepperRKAppAction<double>::BEGIN_STAGE:
      case StepperRKAppAction<double>::BEFORE_SOLVE: {
        const double X_i = get_ele(*(x), 0);
        const double f_i = get_ele(*(xDot), 0);

        if (stageNumber == 0) {
          TEST_FLOATING_EQUALITY(X_i, 1.0, relTol);
          TEST_ASSERT(std::abs(f_i) < relTol);
          TEST_ASSERT(std::abs(time) < relTol);
        }
        else if (stageNumber == 1) {
          TEST_FLOATING_EQUALITY(X_i, 1.0, relTol);
          TEST_FLOATING_EQUALITY(f_i, -1.0, relTol);
          TEST_FLOATING_EQUALITY(time, 1.0, relTol);
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPT(!(-1 < stageNumber && stageNumber < 2));
        }

        break;
      }
      case StepperRKAppAction<double>::AFTER_SOLVE:
      case StepperRKAppAction<double>::BEFORE_EXPLICIT_EVAL:
      case StepperRKAppAction<double>::END_STAGE: {
        const double X_i = get_ele(*(x), 0);
        const double f_i = get_ele(*(xDot), 0);

        if (stageNumber == 0) {
          // X_i = 1, f_1 = 0
          TEST_FLOATING_EQUALITY(X_i, 1.0, relTol);
          TEST_FLOATING_EQUALITY(f_i, -1.0, relTol);
          TEST_ASSERT(std::abs(time) < relTol);
        }
        else if (stageNumber == 1) {
          // X_i = , f_i =
          TEST_FLOATING_EQUALITY(X_i, 1.0 / 3.0, relTol);
          TEST_FLOATING_EQUALITY(f_i, -1.0 / 3.0, relTol);
          TEST_FLOATING_EQUALITY(time, 1.0, relTol);
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPT(!(-1 < stageNumber && stageNumber < 2));
        }

        break;
      }
      case StepperRKAppAction<double>::END_STEP: {
        const double x_1 = get_ele(*(x), 0);
        time             = workingState->getTime();
        TEST_FLOATING_EQUALITY(x_1, 1.0 / 3.0, relTol);  // Should be x_1
        TEST_FLOATING_EQUALITY(time, 1.0, relTol);
        TEST_FLOATING_EQUALITY(dt, 1.0, relTol);
        TEST_COMPARE(stageNumber, ==, -1);

        if (stepper->getUseEmbedded() == true) {
          TEST_FLOATING_EQUALITY(workingState->getTolRel(), 1.0, relTol);
          TEST_ASSERT(std::abs(workingState->getTolAbs()) < relTol);
          // e = 0 from doxygen above.
          TEST_ASSERT(std::abs(workingState->getErrorRel()) < relTol);
        }
      }
    }
  }

 private:
  Teuchos::FancyOStream &out;
  bool &success;
};

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK, EDIRK_Trapezoida_Modifier)
{
  auto stepper = rcp(new Tempus::StepperEDIRK_TrapezoidalRule<double>());
  Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::DahlquistTestModel<double>());

  auto modifier = rcp(new StepperRKModifierEDIRK_TrapezoidaTest(out, success));

  stepper->setModel(model);
  stepper->setAppAction(modifier);
  stepper->setICConsistency("Consistent");
  stepper->setUseFSAL(false);
  stepper->initialize();

  // Create a SolutionHistory.
  auto solutionHistory = Tempus::createSolutionHistoryME(model);

  // Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 1.0;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  solutionHistory->getWorkingState()->setTime(dt);
  stepper->takeStep(solutionHistory);  // Primary testing occurs in
                                       // modifierX during takeStep().
  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}
}  // namespace Tempus_Unit_Test
