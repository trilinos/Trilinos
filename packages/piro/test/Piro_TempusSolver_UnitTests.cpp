// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"

#include "Tpetra_Core.hpp"

#ifdef HAVE_PIRO_TEMPUS
#include "Piro_TempusSolver.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Piro_ObserverToTempusIntegrationObserverAdapter.hpp"
#include "Piro_Test_MockTempusObserver.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_NOXSolver.hpp"
#endif /* HAVE_PIRO_NOX */

#include "Piro_Test_ThyraSupport.hpp"
#include "Piro_Test_WeakenedModelEvaluator.hpp"
#include "Piro_Test_MockObserver.hpp"
#include "Piro_TempusIntegrator.hpp"

#include "MockModelEval_A_Tpetra.hpp"
#include "Piro_StratimikosUtils.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"

#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include <stdexcept>


using namespace Teuchos;
using namespace Piro;
using namespace Piro::Test;

namespace Thyra {
  typedef ModelEvaluatorBase MEB;
} // namespace Thyra

// Setup support

RCP<Thyra::ModelEvaluator<double> > defaultModelNew()
{
  auto comm = Tpetra::getDefaultComm();
  auto model = rcp(new MockModelEval_A_Tpetra(comm));

  std::string inputFile = "input_tempus_be_nox_solver.yaml";
  auto piroParams = rcp(new Teuchos::ParameterList("Piro Parameters"));
  Teuchos::updateParametersFromYamlFile(inputFile, piroParams.ptr());

  auto stratParams = Piro::extractStratimikosParams(piroParams);

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  linearSolverBuilder.setParameterList(stratParams);

  auto lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

  return rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(model, lowsFactory));
}

const RCP<TempusSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluator<double> > &thyraModel,
    double finalTime, 
    const std::string sens_method_string)
{
  const RCP<ParameterList> tempusPL(new ParameterList("Tempus"));
  tempusPL->set("Integrator Name", "Demo Integrator");
  tempusPL->sublist("Demo Integrator").set("Integrator Type", "Integrator Basic");
  tempusPL->sublist("Demo Integrator").set("Stepper Name", "Demo Stepper");
  tempusPL->sublist("Demo Integrator").sublist("Solution History").set("Storage Type", "Unlimited");
  tempusPL->sublist("Demo Integrator").sublist("Solution History").set("Storage Limit", 20);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Initial Time", 0.0);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Final Time", finalTime);
  tempusPL->sublist("Demo Stepper").set("Stepper Type", "Backward Euler");
  tempusPL->sublist("Demo Stepper").set("Zero Initial Guess", false);
  tempusPL->sublist("Demo Stepper").set("Solver Name", "Demo Solver");
  tempusPL->sublist("Demo Stepper").sublist("Demo Solver").sublist("NOX").sublist("Direction").set("Method","Newton");
  SENS_METHOD sens_method; 
  if (sens_method_string == "None") sens_method = Piro::NONE; 
  else if (sens_method_string == "Forward") sens_method = Piro::FORWARD; 
  else if (sens_method_string == "Adjoint") sens_method = Piro::ADJOINT;
  else TEUCHOS_ASSERT(false);
  Teuchos::RCP<Piro::TempusIntegrator<double> > integrator 
      = Teuchos::rcp(new Piro::TempusIntegrator<double>(tempusPL, thyraModel, sens_method));
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver = Teuchos::null;

  RCP<ParameterList> stepperPL = Teuchos::rcp(&(tempusPL->sublist("Demo Stepper")), false);

  RCP<Tempus::StepperFactory<double> > sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  const RCP<Tempus::Stepper<double> > stepper = sf->createStepper(stepperPL, thyraModel);
  return rcp(new TempusSolver<double>(integrator, stepper, stepSolver, thyraModel, finalTime, sens_method_string));
}

const RCP<TempusSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluator<double> > &thyraModel,
    double initialTime,
    double finalTime,
    const RCP<Piro::ObserverBase<double> > &observer,
    const std::string sens_method_string)
{
  const RCP<ParameterList> tempusPL(new ParameterList("Tempus"));
  tempusPL->set("Integrator Name", "Demo Integrator");
  tempusPL->sublist("Demo Integrator").set("Integrator Type", "Integrator Basic");
  tempusPL->sublist("Demo Integrator").set("Stepper Name", "Demo Stepper");
  tempusPL->sublist("Demo Integrator").sublist("Solution History").set("Storage Type", "Unlimited");
  tempusPL->sublist("Demo Integrator").sublist("Solution History").set("Storage Limit", 20);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Initial Time", initialTime);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Final Time", finalTime);
  tempusPL->sublist("Demo Stepper").set("Stepper Type", "Backward Euler");
  tempusPL->sublist("Demo Stepper").set("Zero Initial Guess", false);
  tempusPL->sublist("Demo Stepper").set("Solver Name", "Demo Solver");
  tempusPL->sublist("Demo Stepper").sublist("Demo Solver").sublist("NOX").sublist("Direction").set("Method","Newton");
  SENS_METHOD sens_method; 
  if (sens_method_string == "None") sens_method = Piro::NONE; 
  else if (sens_method_string == "Forward") sens_method = Piro::FORWARD; 
  else if (sens_method_string == "Adjoint") sens_method = Piro::ADJOINT;
  else TEUCHOS_ASSERT(false);
  Teuchos::RCP<Piro::TempusIntegrator<double> > integrator 
      = Teuchos::rcp(new Piro::TempusIntegrator<double>(tempusPL, thyraModel, sens_method));
  const RCP<const Tempus::SolutionHistory<double> > solutionHistory = integrator->getSolutionHistory();
  const RCP<const Tempus::TimeStepControl<double> > timeStepControl = integrator->getTimeStepControl();

  const Teuchos::RCP<Tempus::IntegratorObserver<double> > tempusObserver = 
	  Teuchos::rcp(new ObserverToTempusIntegrationObserverAdapter<double>(solutionHistory, timeStepControl, observer, 
				  false, false, sens_method));
  integrator->setObserver(tempusObserver);
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver = Teuchos::null;
  RCP<ParameterList> stepperPL = Teuchos::rcp(&(tempusPL->sublist("Demo Stepper")), false);
  RCP<Tempus::StepperFactory<double> > sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  const RCP<Tempus::Stepper<double> > stepper = sf->createStepper(stepperPL, thyraModel);

  return rcp(new TempusSolver<double>(integrator, stepper, stepSolver, thyraModel, initialTime, finalTime, sens_method_string));
}

const RCP<TempusSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluator<double> > &thyraModel,
    double initialTime,
    double finalTime,
    double fixedTimeStep,
    const RCP<Piro::ObserverBase<double> > &observer,
    const std::string sens_method_string)
{
  const RCP<ParameterList> tempusPL(new ParameterList("Tempus"));
  tempusPL->set("Integrator Name", "Demo Integrator");
  tempusPL->sublist("Demo Integrator").set("Integrator Type", "Integrator Basic");
  tempusPL->sublist("Demo Integrator").set("Stepper Name", "Demo Stepper");
  tempusPL->sublist("Demo Integrator").sublist("Solution History").set("Storage Type", "Unlimited");
  tempusPL->sublist("Demo Integrator").sublist("Solution History").set("Storage Limit", 20);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Initial Time", initialTime);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Final Time", finalTime);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Minimum Time Step", fixedTimeStep);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Initial Time Step", fixedTimeStep);
  tempusPL->sublist("Demo Integrator").sublist("Time Step Control").set("Maximum Time Step", fixedTimeStep);
  tempusPL->sublist("Demo Stepper").set("Stepper Type", "Backward Euler");
  tempusPL->sublist("Demo Stepper").set("Zero Initial Guess", false);
  tempusPL->sublist("Demo Stepper").set("Solver Name", "Demo Solver");
  tempusPL->sublist("Demo Stepper").sublist("Demo Solver").sublist("NOX").sublist("Direction").set("Method","Newton");

  Teuchos::RCP<Piro::TempusIntegrator<double> > integrator = Teuchos::rcp(new Piro::TempusIntegrator<double>(tempusPL, thyraModel));
  const RCP<const Tempus::SolutionHistory<double> > solutionHistory = integrator->getSolutionHistory();
  const RCP<const Tempus::TimeStepControl<double> > timeStepControl = integrator->getTimeStepControl();
  SENS_METHOD sens_method; 
  if (sens_method_string == "None") sens_method = Piro::NONE; 
  else if (sens_method_string == "Forward") sens_method = Piro::FORWARD; 
  else if (sens_method_string == "Adjoint") sens_method = Piro::ADJOINT;
  else TEUCHOS_ASSERT(false);
  const Teuchos::RCP<Tempus::IntegratorObserver<double> > tempusObserver 
      = Teuchos::rcp(new ObserverToTempusIntegrationObserverAdapter<double>(solutionHistory, 
			      timeStepControl, observer, false, false, sens_method));
  integrator->setObserver(tempusObserver);
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver = Teuchos::null;
  RCP<ParameterList> stepperPL = Teuchos::rcp(&(tempusPL->sublist("Demo Stepper")), false);
  RCP<Tempus::StepperFactory<double> > sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  const RCP<Tempus::Stepper<double> > stepper = sf->createStepper(stepperPL, thyraModel);

  return rcp(new TempusSolver<double>(integrator, stepper, stepSolver, thyraModel, initialTime, finalTime, sens_method_string));
}


Thyra::ModelEvaluatorBase::InArgs<double> getStaticNominalValues(const Thyra::ModelEvaluator<double> &model)
{
  Thyra::ModelEvaluatorBase::InArgs<double> result = model.getNominalValues();
  if (result.supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot)) {
    result.set_x_dot(Teuchos::null);
  }
  return result;
}


// Floating point tolerance
const double tol = 1.0e-8;

TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_Solution)
{
  const auto model = defaultModelNew();
  const double finalTime = 0.0;

  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "None");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const RCP<Thyra::VectorBase<double> > solution =
    Thyra::createMember(solver->get_g_space(solutionResponseIndex));
  outArgs.set_g(solutionResponseIndex, solution);

  solver->evalModel(inArgs, outArgs);

  const RCP<const Thyra::VectorBase<double> > initialCondition =
    model->getNominalValues().get_x();

  TEST_COMPARE_FLOATING_ARRAYS(
      arrayFromVector(*solution),
      arrayFromVector(*initialCondition),
      tol);
}


TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_Response)
{
  const auto model = defaultModelNew();

  const int responseIndex = 0;

  const RCP<Thyra::VectorBase<double> > expectedResponse =
    Thyra::createMember(model->get_g_space(responseIndex));
  {
    const Thyra::MEB::InArgs<double> modelInArgs = getStaticNominalValues(*model);
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_g(responseIndex, expectedResponse);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "None");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const RCP<Thyra::VectorBase<double> > response =
    Thyra::createMember(solver->get_g_space(responseIndex));
  outArgs.set_g(responseIndex, response);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = arrayFromVector(*expectedResponse);
  const Array<double> actual = arrayFromVector(*response);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_NoDfDpMv_NoSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model(
      new WeakenedModelEvaluator_NoDfDpMv(defaultModelNew()));

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "None");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  /*const int responseIndex = 0;
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;

  TEST_ASSERT(outArgs.supports(Thyra::MEB::OUT_ARG_DgDp, responseIndex, parameterIndex).none());
  TEST_ASSERT(outArgs.supports(Thyra::MEB::OUT_ARG_DgDp, solutionResponseIndex, parameterIndex).none());*/

  TEST_NOTHROW(solver->evalModel(inArgs, outArgs));
}


TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_NoDgDp_NoResponseSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model(
      new WeakenedModelEvaluator_NoDgDp(defaultModelNew()));

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "None");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int responseIndex = 0;
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;

  TEST_ASSERT(outArgs.supports(Thyra::MEB::OUT_ARG_DgDp, responseIndex, parameterIndex).none());
  TEST_ASSERT(!outArgs.supports(Thyra::MEB::OUT_ARG_DgDp, solutionResponseIndex, parameterIndex).none());

  TEST_NOTHROW(solver->evalModel(inArgs, outArgs));
}


TEUCHOS_UNIT_TEST(Piro_TempusSolver, ObserveInitialCondition)
{
  const auto model = defaultModelNew();
  const RCP<MockObserver<double> > observer(new MockObserver<double>);
  const double timeStamp = 2.0;

  const RCP<TempusSolver<double> > solver = solverNew(model, timeStamp, timeStamp, observer, "None");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  const Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  solver->evalModel(inArgs, outArgs);

  {
    const RCP<const Thyra::VectorBase<double> > solution =
      observer->lastSolution();

    const RCP<const Thyra::VectorBase<double> > initialCondition =
      model->getNominalValues().get_x();

    TEST_COMPARE_FLOATING_ARRAYS(
        arrayFromVector(*solution),
        arrayFromVector(*initialCondition),
        tol);
  }

  TEST_FLOATING_EQUALITY(observer->lastStamp(), timeStamp, tol);
}


TEUCHOS_UNIT_TEST(Piro_TempusSolver, ObserveFinalSolution)
{
  const auto model = defaultModelNew();
  const RCP<MockObserver<double> > observer(new MockObserver<double>);
  const double initialTime = 0.0;
  const double finalTime = 0.1;
  const double timeStepSize = 0.05;

  const RCP<ParameterList> appParams(new ParameterList("params"));
  Teuchos::updateParametersFromYamlFile("input_tempus_be_nox_solver.yaml",appParams.ptr());

  auto& tsc = appParams->sublist("Tempus").sublist("Demo Integrator").sublist("Time Step Control");
  tsc.set("Initial Time", initialTime);
  tsc.set("Final Time", finalTime);
  tsc.set("Minimum Time Step", timeStepSize/10);
  tsc.set("Initial Time Step", timeStepSize);
  tsc.set("Maximum Time Step", timeStepSize);

  RCP<Thyra::ModelEvaluator<double>> adjointModel = Teuchos::null;
  RCP<ObserverBase<double>> piro_obs = observer;
  auto solver = tempusSolver(appParams,model,adjointModel,piro_obs);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const RCP<Thyra::VectorBase<double> > solution =
    Thyra::createMember(solver->get_g_space(solutionResponseIndex));
  outArgs.set_g(solutionResponseIndex, solution);

  solver->evalModel(inArgs, outArgs);

  TEST_COMPARE_FLOATING_ARRAYS(
      arrayFromVector(*observer->lastSolution()),
      arrayFromVector(*solution),
      tol);

  TEST_FLOATING_EQUALITY(observer->lastStamp(), finalTime, tol);
}

TEUCHOS_UNIT_TEST(Piro_TempusSolver, ObserveFinalSolutionNoWrapper)
{
  const auto model = defaultModelNew();
  auto observer = Teuchos::rcp(new MockTempusObserver<double>());

  const double initialTime = 0.0;
  const double finalTime = 0.1;
  const double timeStepSize = 0.05;

  const RCP<ParameterList> appParams(new ParameterList("params"));
  Teuchos::updateParametersFromYamlFile("input_tempus_be_nox_solver.yaml",appParams.ptr());

  auto& tsc = appParams->sublist("Tempus").sublist("Demo Integrator").sublist("Time Step Control");
  tsc.set("Initial Time", initialTime);
  tsc.set("Final Time", finalTime);
  tsc.set("Minimum Time Step", timeStepSize/10);
  tsc.set("Initial Time Step", timeStepSize);
  tsc.set("Maximum Time Step", timeStepSize);

  RCP<Thyra::ModelEvaluator<double>> adjointModel = Teuchos::null;
  RCP<ObserverBase<double>> piro_obs = observer;
  auto solver = tempusSolver(appParams,model,adjointModel,piro_obs);

  // Check that the piro integrator did NOT wrap our observer in a tempus adapter
  auto piro_integrator = solver->getPiroTempusIntegrator();
  auto tempus_integrator = piro_integrator->getIntegrator();
  auto basic_integrator = Teuchos::rcp_dynamic_cast<const Tempus::IntegratorBasic<double>>(tempus_integrator);
  auto nonconst_basic_integrator = Teuchos::rcp_const_cast<Tempus::IntegratorBasic<double>>(basic_integrator);
  auto tempus_observer = nonconst_basic_integrator->getObserver();
  TEST_ASSERT (tempus_observer==observer);

  // Proceed like the test above
  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const RCP<Thyra::VectorBase<double> > solution =
    Thyra::createMember(solver->get_g_space(solutionResponseIndex));
  outArgs.set_g(solutionResponseIndex, solution);

  solver->evalModel(inArgs, outArgs);

  TEST_COMPARE_FLOATING_ARRAYS(
      arrayFromVector(*observer->lastSolution()),
      arrayFromVector(*solution),
      tol);

  TEST_FLOATING_EQUALITY(observer->lastStamp(), finalTime, tol);
}

// The initial state is made consistent (e.g., the algebraic components of a DAE) by the
// function set with setConsistentInitialStateFunction, for the parameters of evalModel
TEUCHOS_UNIT_TEST(Piro_TempusSolver, ConsistentInitialState)
{
  const auto model = defaultModelNew();
  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "None");

  // Parameters different from the nominal ones
  Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  const RCP<Thyra::VectorBase<double> > p = inArgs.get_p(0)->clone_v();
  Thyra::assign(p.ptr(), 1.5);
  inArgs.set_p(0, p);

  const RCP<const Thyra::VectorBase<double> > nominalState = model->getNominalValues().get_x();

  // "Consistent" initial state: x0 <- 2 x0 + 1
  int numCalls = 0;
  solver->setConsistentInitialStateFunction(
    [&](const RCP<Thyra::ModelEvaluator<double> >& icModel,
        const Thyra::MEB::InArgs<double>& icInArgs,
        Thyra::VectorBase<double>& x0) {
      ++numCalls;
      TEST_ASSERT(icModel.get() == model.get());
      TEST_EQUALITY_CONST(icInArgs.get_t(), 0.0);
      TEST_COMPARE_FLOATING_ARRAYS(arrayFromVector(*icInArgs.get_p(0)), arrayFromVector(*p), tol);
      // Each evaluation starts from the nominal initial state
      TEST_COMPARE_FLOATING_ARRAYS(arrayFromVector(x0), arrayFromVector(*nominalState), tol);
      Thyra::Vt_S(Teuchos::ptrFromRef(x0), 2.0);
      Thyra::Vp_S(Teuchos::ptrFromRef(x0), 1.0);
    });

  const RCP<Thyra::VectorBase<double> > expected = nominalState->clone_v();
  Thyra::Vt_S(expected.ptr(), 2.0);
  Thyra::Vp_S(expected.ptr(), 1.0);

  const int solutionResponseIndex = solver->Ng() - 1;
  for (int eval = 1; eval <= 2; ++eval) {
    Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
    const RCP<Thyra::VectorBase<double> > solution =
      Thyra::createMember(solver->get_g_space(solutionResponseIndex));
    outArgs.set_g(solutionResponseIndex, solution);

    solver->evalModel(inArgs, outArgs);

    TEST_EQUALITY(numCalls, eval);
    TEST_COMPARE_FLOATING_ARRAYS(arrayFromVector(*solution), arrayFromVector(*expected), tol);
  }
}


// An initial state set with setInitialState is used as it is
TEUCHOS_UNIT_TEST(Piro_TempusSolver, ConsistentInitialState_SetInitialState)
{
  const auto model = defaultModelNew();
  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "None");

  int numCalls = 0;
  solver->setConsistentInitialStateFunction(
    [&numCalls](const RCP<Thyra::ModelEvaluator<double> >&,
                const Thyra::MEB::InArgs<double>&,
                Thyra::VectorBase<double>& x0) {
      ++numCalls;
      Thyra::assign(Teuchos::ptrFromRef(x0), 5.0);
    });

  const RCP<Thyra::VectorBase<double> > initialState = model->getNominalValues().get_x()->clone_v();
  Thyra::assign(initialState.ptr(), 7.0);
  solver->setInitialState(0.0, initialState);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const RCP<Thyra::VectorBase<double> > solution =
    Thyra::createMember(solver->get_g_space(solutionResponseIndex));
  outArgs.set_g(solutionResponseIndex, solution);

  solver->evalModel(inArgs, outArgs);

  TEST_EQUALITY_CONST(numCalls, 0);
  TEST_COMPARE_FLOATING_ARRAYS(arrayFromVector(*solution), arrayFromVector(*initialState), tol);
}


// The Tempus-driven objective with a consistent initial state is tested in
// Piro_ThyraProductME_Tempus_FinalObjective_UnitTests.cpp (one process, explicit adjoint model).

#endif /* HAVE_PIRO_TEMPUS */
