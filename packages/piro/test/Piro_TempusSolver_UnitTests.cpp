// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_TEMPUS
#include "Piro_TempusSolver.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Piro_ObserverToTempusIntegrationObserverAdapter.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_NOXSolver.hpp"
#endif /* HAVE_PIRO_NOX */

#include "Piro_Test_ThyraSupport.hpp"
#include "Piro_Test_WeakenedModelEvaluator.hpp"
#include "Piro_Test_MockObserver.hpp"
#include "Piro_TempusIntegrator.hpp"

#include "MockModelEval_A.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"

#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"

#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

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

const RCP<EpetraExt::ModelEvaluator> epetraModelNew()
{
#ifdef HAVE_MPI
  const MPI_Comm comm = MPI_COMM_WORLD;
#else /*HAVE_MPI*/
  const int comm = 0;
#endif /*HAVE_MPI*/
  return rcp(new MockModelEval_A(comm));
}

const RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyraModelNew(const RCP<EpetraExt::ModelEvaluator> &epetraModel)
{
  const RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory(new Thyra::AmesosLinearOpWithSolveFactory);
  return epetraModelEvaluator(epetraModel, lowsFactory);
}

RCP<Thyra::ModelEvaluatorDefaultBase<double> > defaultModelNew()
{
  return thyraModelNew(epetraModelNew());
}

const RCP<TempusSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
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
  Teuchos::RCP<Piro::TempusIntegrator<double> > integrator 
      = Teuchos::rcp(new Piro::TempusIntegrator<double>(tempusPL, thyraModel, sens_method));
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver = Teuchos::null;

  RCP<ParameterList> stepperPL = Teuchos::rcp(&(tempusPL->sublist("Demo Stepper")), false);

  RCP<Tempus::StepperFactory<double> > sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  const RCP<Tempus::Stepper<double> > stepper = sf->createStepper(stepperPL, thyraModel);
  return rcp(new TempusSolver<double>(integrator, stepper, stepSolver, thyraModel, finalTime, sens_method_string));
}

const RCP<TempusSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
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
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
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
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
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
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();

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
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
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
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const RCP<MockObserver<double> > observer(new MockObserver<double>);
  const double initialTime = 0.0;
  const double finalTime = 0.1;
  const double timeStepSize = 0.05;

  const RCP<TempusSolver<double> > solver =
    solverNew(model, initialTime, finalTime, timeStepSize, observer, "None");

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

#endif /* HAVE_PIRO_TEMPUS */
