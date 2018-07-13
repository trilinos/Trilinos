/*
// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER
*/

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_TEMPUS
#include "Piro_TempusSolver.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
#include "Piro_ObserverToTempusIntegrationObserverAdapter.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_NOXSolver.hpp"
#endif /* HAVE_PIRO_NOX */

#include "Piro_Test_ThyraSupport.hpp"
#include "Piro_Test_WeakenedModelEvaluator.hpp"
#include "Piro_Test_MockObserver.hpp"

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
    double finalTime)
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
  Teuchos::RCP<Tempus::IntegratorBasic<double> > integrator = Tempus::integratorBasic<double>(tempusPL, thyraModel);
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver = Teuchos::null;

  RCP<ParameterList> stepperPL = Teuchos::rcp(&(tempusPL->sublist("Demo Stepper")), false);
  const RCP<Tempus::Stepper<double> > stepper = rcp(new Tempus::StepperBackwardEuler<double>(thyraModel, stepperPL));
  return rcp(new TempusSolver<double>(integrator, stepper, stepSolver, thyraModel, finalTime));
}

const RCP<TempusSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
    double initialTime,
    double finalTime,
    const RCP<Piro::ObserverBase<double> > &observer)
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
  Teuchos::RCP<Tempus::IntegratorBasic<double> > integrator = Tempus::integratorBasic<double>(tempusPL, thyraModel);
  const RCP<const Tempus::SolutionHistory<double> > solutionHistory = integrator->getSolutionHistory();
  const RCP<const Tempus::TimeStepControl<double> > timeStepControl = integrator->getTimeStepControl();

  const Teuchos::RCP<Tempus::IntegratorObserver<double> > tempusObserver = Teuchos::rcp(new ObserverToTempusIntegrationObserverAdapter<double>(solutionHistory, timeStepControl, observer));
  integrator->setObserver(tempusObserver);
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver = Teuchos::null;
  RCP<ParameterList> stepperPL = Teuchos::rcp(&(tempusPL->sublist("Demo Stepper")), false);
  const RCP<Tempus::Stepper<double> > stepper = rcp(new Tempus::StepperBackwardEuler<double>(thyraModel, stepperPL));

  return rcp(new TempusSolver<double>(integrator, stepper, stepSolver, thyraModel, initialTime, finalTime));
}

const RCP<TempusSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
    double initialTime,
    double finalTime,
    double fixedTimeStep,
    const RCP<Piro::ObserverBase<double> > &observer)
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
  Teuchos::RCP<Tempus::IntegratorBasic<double> > integrator = Tempus::integratorBasic<double>(tempusPL, thyraModel);
  const RCP<const Tempus::SolutionHistory<double> > solutionHistory = integrator->getSolutionHistory();
  const RCP<const Tempus::TimeStepControl<double> > timeStepControl = integrator->getTimeStepControl();

  const Teuchos::RCP<Tempus::IntegratorObserver<double> > tempusObserver = Teuchos::rcp(new ObserverToTempusIntegrationObserverAdapter<double>(solutionHistory, timeStepControl, observer));
  integrator->setObserver(tempusObserver);
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver = Teuchos::null;
  RCP<ParameterList> stepperPL = Teuchos::rcp(&(tempusPL->sublist("Demo Stepper")), false);
  const RCP<Tempus::Stepper<double> > stepper = rcp(new Tempus::StepperBackwardEuler<double>(thyraModel, stepperPL));

  return rcp(new TempusSolver<double>(integrator, stepper, stepSolver, thyraModel, initialTime, finalTime));
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

  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime);

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
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime);

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
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int responseIndex = 0;
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;

  TEST_ASSERT(outArgs.supports(Thyra::MEB::OUT_ARG_DgDp, responseIndex, parameterIndex).none());
  TEST_ASSERT(outArgs.supports(Thyra::MEB::OUT_ARG_DgDp, solutionResponseIndex, parameterIndex).none());

  TEST_NOTHROW(solver->evalModel(inArgs, outArgs));
}


TEUCHOS_UNIT_TEST(Piro_TempusSolver, TimeZero_NoDgDp_NoResponseSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model(
      new WeakenedModelEvaluator_NoDgDp(defaultModelNew()));

  const double finalTime = 0.0;
  const RCP<TempusSolver<double> > solver = solverNew(model, finalTime);

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
  //IKT, FIXME: set to 2.0 instead of 0.0 -- does not work correctly in this case
  const double timeStamp = 0.0;

  const RCP<TempusSolver<double> > solver = solverNew(model, timeStamp, timeStamp, observer);

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
    solverNew(model, initialTime, finalTime, timeStepSize, observer);

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
