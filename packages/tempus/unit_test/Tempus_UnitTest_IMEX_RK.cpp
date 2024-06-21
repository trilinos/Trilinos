//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"

#include "Tempus_StepperRKButcherTableau.hpp"

#include "../TestModels/DahlquistTestModel.hpp"
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
using Tempus::StepperFactory;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK, Default_Construction)
{
  // Setup the IMEX Pair ModelEvaluator
  auto explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
      explicitModel, implicitModel));

  // Default construction.
  auto stepper = rcp(new Tempus::StepperIMEX_RK<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier  = rcp(new Tempus::StepperRKModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperRKModifierXDefault<double>());
  auto observer  = rcp(new Tempus::StepperRKObserverDefault<double>());
  auto solver    = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();
  std::string stepperType   = "IMEX RK SSP2";
  auto stepperERK           = Teuchos::rcp(new Tempus::StepperERK_Trapezoidal<double>());
  auto explicitTableau      = stepperERK->getTableau();
  auto stepperSDIRK =
      Teuchos::rcp(new Tempus::StepperSDIRK_2Stage3rdOrder<double>());
  auto implicitTableau = stepperSDIRK->getTableau();
  int order            = 2;

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

  stepper->setStepperName(stepperType);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setExplicitTableau(explicitTableau);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setImplicitTableau(implicitTableau);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setOrder(order);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  TEUCHOS_TEST_FOR_EXCEPT(explicitTableau != stepper->getTableau());
  TEUCHOS_TEST_FOR_EXCEPT(explicitTableau != stepper->getExplicitTableau());
  TEUCHOS_TEST_FOR_EXCEPT(implicitTableau != stepper->getImplicitTableau());

  // Full argument list construction.
  stepper = rcp(new Tempus::StepperIMEX_RK<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      zeroInitialGuess, modifier, stepperType, explicitTableau, implicitTableau,
      order));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK, StepperFactory_Construction)
{
  // Setup the IMEX Pair ModelEvaluator
  auto explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
      explicitModel, implicitModel));

  testFactoryConstruction("IMEX RK SSP2", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK, StepperFactory_Construction_General_wo_Parameterlist)
{
  // Setup the IMEX Pair ModelEvaluator
  auto explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
      explicitModel, implicitModel));

  RCP<StepperFactory<double>> sf = Teuchos::rcp(new StepperFactory<double>());

  auto stepper = sf->createStepper("General IMEX RK", model);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK,
                  StepperFactory_Construction_General_wo_Parameterlist_Model)
{
  // Setup the IMEX Pair ModelEvaluator
  auto explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
      explicitModel, implicitModel));

  RCP<StepperFactory<double>> sf = Teuchos::rcp(new StepperFactory<double>());

  auto stepper = sf->createStepper("General IMEX RK");
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK, AppAction)
{
  auto stepper = rcp(new Tempus::StepperIMEX_RK<double>());
  // Setup the IMEX Pair ModelEvaluator
  auto explicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel =
      rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
      explicitModel, implicitModel));

  testRKAppAction(stepper, model, out, success);
}

class StepperRKModifierIMEX_TrapezoidaTest
  : virtual public Tempus::StepperRKModifierBase<double> {
 public:
  /// Constructor
  StepperRKModifierIMEX_TrapezoidaTest(Teuchos::FancyOStream &Out,
                                       bool &Success)
    : out(Out), success(Success)
  {
  }

  // FSAL
  /// Destructor
  virtual ~StepperRKModifierIMEX_TrapezoidaTest() {}

  /// Test the modify RK Stepper at the Action Locations.
  virtual void modify(
      Teuchos::RCP<Tempus::SolutionHistory<double>> sh,
      Teuchos::RCP<Tempus::StepperRKBase<double>> stepper,
      const typename Tempus::StepperRKAppAction<double>::ACTION_LOCATION actLoc)
  {
    const double relTol = 1.0e-14;
    auto stepper_imex =
        Teuchos::rcp_dynamic_cast<const Tempus::StepperIMEX_RK<double>>(stepper,
                                                                        true);
    auto stageNumber                          = stepper->getStageNumber();
    Teuchos::SerialDenseVector<int, double> c = stepper_imex->getTableau()->c();
    Teuchos::SerialDenseVector<int, double> chat =
        stepper_imex->getImplicitTableau()->c();

    auto currentState = sh->getCurrentState();
    auto workingState = sh->getWorkingState();
    const double dt   = workingState->getTimeStep();
    double time       = currentState->getTime();
    double imp_time   = time;
    if (stageNumber >= 0) {
      time += c(stageNumber) * dt;
      imp_time += chat(stageNumber) * dt;
    }

    auto x    = workingState->getX();
    auto xDot = workingState->getXDot();
    if (xDot == Teuchos::null) xDot = stepper->getStepperXDot();

    switch (actLoc) {
      case StepperRKAppAction<double>::BEGIN_STEP: {
        {
          auto imex_me = Teuchos::rcp_dynamic_cast<
              const Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>>(
              stepper->getModel(), true);
          auto explicitModel = Teuchos::rcp_dynamic_cast<
              const Tempus_Test::DahlquistTestModel<double>>(
              imex_me->getExplicitModel(), true);
          auto implicitModel = Teuchos::rcp_dynamic_cast<
              const Tempus_Test::DahlquistTestModel<double>>(
              imex_me->getImplicitModel(), true);

          TEST_FLOATING_EQUALITY(explicitModel->getLambda(), 1.0, relTol);
          TEST_FLOATING_EQUALITY(implicitModel->getLambda(), 2.0, relTol);
        }
        TEST_FLOATING_EQUALITY(dt, 1.0, relTol);

        const double x_0 = get_ele(*(x), 0);
        TEST_FLOATING_EQUALITY(x_0, 1.0, relTol);
        TEST_ASSERT(std::abs(time) < relTol);
        TEST_ASSERT(std::abs(imp_time) < relTol);
        TEST_FLOATING_EQUALITY(dt, 1.0, relTol);
        TEST_COMPARE(stageNumber, ==, -1);
        break;
      }
      case StepperRKAppAction<double>::BEGIN_STAGE:
      case StepperRKAppAction<double>::BEFORE_SOLVE: {
        const double X_i = get_ele(*(x), 0);
        const double f_i = get_ele(*(xDot), 0);

        if (stageNumber == 0) {
          TEST_FLOATING_EQUALITY(X_i, 1.0, relTol);                       // 1
          TEST_FLOATING_EQUALITY(f_i, 1.0, relTol);                       // 1
          TEST_FLOATING_EQUALITY(imp_time, 0.78867513459481275, relTol);  // 1
          TEST_ASSERT(std::abs(time) < relTol);
        }
        else if (stageNumber == 1) {
          TEST_FLOATING_EQUALITY(X_i, -std::sqrt(3), relTol);  // -sqrt(3)
          TEST_FLOATING_EQUALITY(f_i, 1.0, relTol);            // 1
          TEST_FLOATING_EQUALITY(time, 1.0, relTol);
          TEST_FLOATING_EQUALITY(imp_time, 0.21132486540518725, relTol);  // 1
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
          TEST_FLOATING_EQUALITY(X_i, -std::sqrt(3), relTol);             // -sqrt(3)
          TEST_FLOATING_EQUALITY(f_i, 1.0, relTol);                       // 1
          TEST_FLOATING_EQUALITY(imp_time, 0.78867513459481275, relTol);  // 1
          TEST_ASSERT(std::abs(time) < relTol);
        }
        else if (stageNumber == 1) {
          // X_i = , f_i =
          TEST_FLOATING_EQUALITY(X_i, 3.0 - 3.0 * std::sqrt(3.0),
                                 relTol);            // -3sqrt(3) - 3
          TEST_FLOATING_EQUALITY(f_i, 1.0, relTol);  // 1
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
        TEST_FLOATING_EQUALITY(x_1, -(6.0 * std::sqrt(3) - 11.0 / 2.0),
                               relTol);  // -( 6sqrt(3) - 11/2)
        TEST_FLOATING_EQUALITY(time, 1.0, relTol);
        TEST_FLOATING_EQUALITY(dt, 1.0, relTol);
        TEST_COMPARE(stageNumber, ==, -1);
      }
    }
  }

 private:
  Teuchos::FancyOStream &out;
  bool &success;
};

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK, IMEX_RK_Modifier)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<double>> explicitModel =
      rcp(new Tempus_Test::DahlquistTestModel<double>(1.0, true));
  Teuchos::RCP<const Thyra::ModelEvaluator<double>> implicitModel =
      rcp(new Tempus_Test::DahlquistTestModel<double>(2.0, true));

  auto model = rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
      explicitModel, implicitModel));

  auto modifier = rcp(new StepperRKModifierIMEX_TrapezoidaTest(out, success));

  // Default construction.
  auto stepper = rcp(new Tempus::StepperIMEX_RK<double>());
  stepper->setModel(model);

  // Default values for construction.
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSAL();
  std::string ICConsistency = stepper->getICConsistency();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheck();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();
  std::string stepperType   = "IMEX RK SSP2";
  auto stepperERK           = Teuchos::rcp(new Tempus::StepperERK_Trapezoidal<double>());
  auto explicitTableau      = stepperERK->getTableau();
  auto stepperSDIRK =
      Teuchos::rcp(new Tempus::StepperSDIRK_2Stage3rdOrder<double>());
  auto implicitTableau = stepperSDIRK->getTableau();
  int order            = 2;

  stepper->setStepperName(stepperType);
  stepper->setExplicitTableau(explicitTableau);
  stepper->setImplicitTableau(implicitTableau);
  stepper->setOrder(order);
  stepper->setSolver(solver);
  stepper->setUseFSAL(useFSAL);
  stepper->setICConsistency(ICConsistency);
  stepper->setICConsistencyCheck(ICConsistencyCheck);
  stepper->setZeroInitialGuess(zeroInitialGuess);

  stepper->setModel(model);
  stepper->setAppAction(modifier);
  stepper->setUseFSAL(false);
  stepper->initialize();

  // Create a SolutionHistory.
  auto solutionHistory = Tempus::createSolutionHistoryME(explicitModel);

  //// Take one time step.
  stepper->setInitialConditions(solutionHistory);
  solutionHistory->initWorkingState();
  double dt = 1.0;
  solutionHistory->getWorkingState()->setTimeStep(dt);
  solutionHistory->getWorkingState()->setTime(dt);
  stepper->takeStep(solutionHistory);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}

}  // namespace Tempus_Unit_Test
