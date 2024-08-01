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
#include "Piro_Helpers.hpp" 

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
#include<mpi.h>

using namespace Teuchos;
using namespace Piro;
using namespace Piro::Test;

//#define DEBUG_OUTPUT

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
  auto x0 =  thyraModel->getNominalValues().get_x(); 
  const int num_param = thyraModel->get_p_space(0)->dim();
  RCP<Thyra::MultiVectorBase<double> > DxDp0 =
      Thyra::createMembers(thyraModel->get_x_space(), num_param);
  DxDp0->assign(0.0); 
  integrator->initializeSolutionHistory(0.0, x0, Teuchos::null, Teuchos::null,
                                          DxDp0, Teuchos::null, Teuchos::null);
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
  auto x0 =  thyraModel->getNominalValues().get_x(); 
  const int num_param = thyraModel->get_p_space(0)->dim();
  RCP<Thyra::MultiVectorBase<double> > DxDp0 =
      Thyra::createMembers(thyraModel->get_x_space(), num_param);
  DxDp0->assign(0.0); 
  integrator->initializeSolutionHistory(initialTime, x0, Teuchos::null, Teuchos::null,
                                          DxDp0, Teuchos::null, Teuchos::null);
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
  SENS_METHOD sens_method; 
  if (sens_method_string == "None") sens_method = Piro::NONE; 
  else if (sens_method_string == "Forward") sens_method = Piro::FORWARD; 
  else if (sens_method_string == "Adjoint") sens_method = Piro::ADJOINT; 
  Teuchos::RCP<Piro::TempusIntegrator<double> > integrator 
      = Teuchos::rcp(new Piro::TempusIntegrator<double>(tempusPL, thyraModel, sens_method));
  const RCP<const Tempus::SolutionHistory<double> > solutionHistory = integrator->getSolutionHistory();
  const RCP<const Tempus::TimeStepControl<double> > timeStepControl = integrator->getTimeStepControl();

  const Teuchos::RCP<Tempus::IntegratorObserver<double> > tempusObserver 
	  = Teuchos::rcp(new ObserverToTempusIntegrationObserverAdapter<double>(solutionHistory, timeStepControl, observer));
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

TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_DefaultSolutionForwardSensitivity)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::VerboseObjectBase::getDefaultOStream();

  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const double finalTime = 0.0;

  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "Forward");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  solver->evalModel(inArgs, outArgs);
  
  // Test if at 'Final Time'
  double time = solver->getPiroTempusIntegrator()->getTime();
  TEST_FLOATING_EQUALITY(time, 0.0, 1.0e-14);

  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.0, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 0.0, 0.0, 0.0)));
  RCP<const Thyra::MultiVectorBase<double> > DxDp = solver->getPiroTempusIntegrator()->getDxDp();
  TEST_EQUALITY(DxDp->domain()->dim(), expected.size());

  for (int i = 0; i < expected.size(); ++i) {
    TEST_EQUALITY(DxDp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromVector(*DxDp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_DefaultSolutionForwardSensitivityOp)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::VerboseObjectBase::getDefaultOStream();
  
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const double finalTime = 0.0;

  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "Forward");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    solver->create_DgDp_op(solutionResponseIndex, parameterIndex);
  const RCP<Thyra::LinearOpBase<double> > dxdp = dxdp_deriv.getLinearOp();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  solver->evalModel(inArgs, outArgs);
  
  // Test if at 'Final Time'
  double time = solver->getPiroTempusIntegrator()->getTime();
  TEST_FLOATING_EQUALITY(time, 0.0, 1.0e-14);
  
  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.0, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 0.0, 0.0, 0.0)));
  RCP<const Thyra::MultiVectorBase<double> > DxDp = solver->getPiroTempusIntegrator()->getDxDp();
  TEST_EQUALITY(DxDp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
  TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromLinOp(*DxDp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_DefaultResponseForwardSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();

  const int responseIndex = 0;
  const int parameterIndex = 0;

  const Thyra::MEB::Derivative<double> dgdp_deriv_expected =
    Thyra::create_DgDp_mv(*model, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp_expected = dgdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> modelInArgs = getStaticNominalValues(*model);
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv_expected);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "Forward");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  TEST_EQUALITY(dgdp->domain()->dim(), dgdp_expected->domain()->dim());
  TEST_EQUALITY(dgdp->range()->dim(), dgdp_expected->range()->dim());
  for (int i = 0; i < dgdp_expected->domain()->dim(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    const Array<double> expected = arrayFromVector(*dgdp_expected->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_DefaultResponseForwardSensitivity_NoDgDxMv)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::VerboseObjectBase::getDefaultOStream();
  
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model(
    new WeakenedModelEvaluator_NoDgDxMv(defaultModelNew()));

  const int responseIndex = 0;
  const int parameterIndex = 0;

  const Thyra::MEB::Derivative<double> dgdp_deriv_expected =
    Thyra::create_DgDp_mv(*model, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp_expected = dgdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> modelInArgs = getStaticNominalValues(*model);
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv_expected);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "Forward");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);
  
  TEST_EQUALITY(dgdp->domain()->dim(), dgdp_expected->domain()->dim());
  TEST_EQUALITY(dgdp->range()->dim(), dgdp_expected->range()->dim());
  for (int i = 0; i < dgdp_expected->domain()->dim(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    const Array<double> expected = arrayFromVector(*dgdp_expected->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}


TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_ResponseAndDefaultForwardSensitivities)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::VerboseObjectBase::getDefaultOStream();
  
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();

  const int responseIndex = 0;
  const int parameterIndex = 0;

  const RCP<Thyra::VectorBase<double> > expectedResponse =
    Thyra::createMember(model->get_g_space(responseIndex));
  {
    const Thyra::MEB::InArgs<double> modelInArgs = getStaticNominalValues(*model);
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_g(responseIndex, expectedResponse);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const Thyra::MEB::Derivative<double> dgdp_deriv_expected =
    Thyra::create_DgDp_mv(*model, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp_expected = dgdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> modelInArgs = getStaticNominalValues(*model);
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv_expected);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "Forward");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  // Requesting response
  const RCP<Thyra::VectorBase<double> > response =
    Thyra::createMember(solver->get_g_space(responseIndex));
  outArgs.set_g(responseIndex, response);

  // Requesting response sensitivity
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  // Run solver
  solver->evalModel(inArgs, outArgs);

  // Checking response
  {
    const Array<double> expected = arrayFromVector(*expectedResponse);
    const Array<double> actual = arrayFromVector(*response);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }

  // Checking sensitivity
  {
    TEST_EQUALITY(dgdp->domain()->dim(), dgdp_expected->domain()->dim());
    TEST_EQUALITY(dgdp->range()->dim(), dgdp_expected->range()->dim());
    for (int i = 0; i < dgdp_expected->domain()->dim(); ++i) {
      const Array<double> actual = arrayFromVector(*dgdp->col(i));
      const Array<double> expected = arrayFromVector(*dgdp_expected->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
    }
  }
}

TEUCHOS_UNIT_TEST(Piro_TempusSolver, ObserveInitialConditionWhenForwardSensitivitiesRequested)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::VerboseObjectBase::getDefaultOStream();
  
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const RCP<MockObserver<double> > observer(new MockObserver<double>);
  const double timeStamp = 2.0;

  const RCP<TempusSolver<double> > solver = solverNew(model, timeStamp, timeStamp, observer, "Forward");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
      Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
    outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);
  solver->evalModel(inArgs, outArgs);

  const RCP<const Thyra::VectorBase<double> > solution =
    observer->lastSolution();
    
  const RCP<const Thyra::MultiVectorBase<double> > solution_dxdp =
    observer->lastSolution_dxdp();

  const RCP<const Thyra::VectorBase<double> > initialCondition =
    model->getNominalValues().get_x();

  TEST_COMPARE_FLOATING_ARRAYS(
      arrayFromVector(*solution),
      arrayFromVector(*initialCondition),
      tol);

  //Test observer output of lastSolution_dxdp 
  for (int np = 0; np < dxdp->domain()->dim(); np++) { 
    Teuchos::RCP<const Thyra::VectorBase<double>> solution_dxdp_vec = solution_dxdp->col(np);
    Teuchos::RCP<Thyra::VectorBase<double>> zero_vec = 
        Thyra::createMember(solution_dxdp_vec->space());
    Thyra::put_scalar(0.0, zero_vec.ptr());  
#ifdef DEBUG_OUTPUT 
    auto out_ =Teuchos::VerboseObjectBase::getDefaultOStream(); 
    *out_ << "\n*** Piro::TransientSolver dxdp for p = " << np << " ***\n";
    Teuchos::Range1D range;
    RTOpPack::ConstSubVectorView<double> dxdpv;
    solution_dxdp_vec->acquireDetachedView(range, &dxdpv);
    auto dxdpa = dxdpv.values();
    for (auto j = 0; j < dxdpa.size(); ++j) *out_ << dxdpa[j] << " ";
    *out_ << "\n*** Piro::TransientSolver dxdp for p = " << np << " ***\n";
#endif
    TEST_COMPARE_FLOATING_ARRAYS(
        arrayFromVector(*solution_dxdp->col(np)),
        arrayFromVector(*zero_vec),
        tol);
  }

  TEST_FLOATING_EQUALITY(observer->lastStamp(), timeStamp, tol);
}

TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_DefaultResponseAdjointSensitivity)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::VerboseObjectBase::getDefaultOStream();

  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const double finalTime = 0.0;

  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "Adjoint");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  solver->evalModel(inArgs, outArgs);
  
  // Test if at 'Final Time'
  double time = solver->getPiroTempusIntegrator()->getTime();
  TEST_FLOATING_EQUALITY(time, 0.0, 1.0e-14);

  std::vector<double> expected(2);
  expected[0] = 2.0; expected[1] = 2.0; //Exact solution for DgDp = [2,2]
  //Specifically, we have dg/dp0 = -sum(x_i) + 2*p0 + p1 + 11 and dg/dp1 = -sum(x_i) + p0 + p1 + 12
  //evaluated at the initial condition x = [3,3,3,3] and p = [1,1] 
  RCP<const Thyra::MultiVectorBase<double> > DgDp = solver->getPiroTempusIntegrator()->getDgDp();
  int expected_size = expected.size();  
  TEST_EQUALITY(DgDp->range()->dim(), expected_size);

  const Array<double> actual = arrayFromVector(*DgDp->col(0));
  const Array<double> expected_array = Array<double>(expected); 
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected_array, tol);
}



TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_ResponseAndDefaultAdjointSensitivities)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::VerboseObjectBase::getDefaultOStream();
  
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();

  const int responseIndex = 0;
  const int parameterIndex = 0;

  const RCP<Thyra::VectorBase<double> > expectedResponse =
    Thyra::createMember(model->get_g_space(responseIndex));
  {
    const Thyra::MEB::InArgs<double> modelInArgs = getStaticNominalValues(*model);
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_g(responseIndex, expectedResponse);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const Thyra::MEB::Derivative<double> dgdp_deriv_expected =
    Thyra::create_DgDp_mv(*model, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp_expected = dgdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> modelInArgs = getStaticNominalValues(*model);
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv_expected);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime, "Adjoint");

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  // Requesting response
  const RCP<Thyra::VectorBase<double> > response =
    Thyra::createMember(solver->get_g_space(responseIndex));
  outArgs.set_g(responseIndex, response);

  // Requesting response sensitivity
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  // Run solver
  solver->evalModel(inArgs, outArgs);

  // Checking response
  {
    const Array<double> expected = arrayFromVector(*expectedResponse);
    const Array<double> actual = arrayFromVector(*response);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }

  // Checking sensitivity
  {
    TEST_EQUALITY(dgdp->domain()->dim(), dgdp_expected->domain()->dim());
    TEST_EQUALITY(dgdp->range()->dim(), dgdp_expected->range()->dim());
    for (int i = 0; i < dgdp_expected->domain()->dim(); ++i) {
      const Array<double> actual = arrayFromVector(*dgdp->col(i));
      const Array<double> expected = arrayFromVector(*dgdp_expected->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
    }
  }
}

#endif /* HAVE_PIRO_TEMPUS */
