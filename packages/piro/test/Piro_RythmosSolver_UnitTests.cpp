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

#ifdef Piro_ENABLE_Rythmos
#include "Piro_RythmosSolver.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_NOXSolver.hpp"
#endif /* Piro_ENABLE_NOX */

#include "Piro_Test_ThyraSupport.hpp"

#include "MockModelEval_A.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"

#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"

#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

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

RCP<Thyra::ModelEvaluatorDefaultBase<double> > defaultModelNew() {
  return thyraModelNew(epetraModelNew());
}

const RCP<RythmosSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
    double finalTime)
{
  const RCP<Rythmos::DefaultIntegrator<double> > integrator =
    Rythmos::defaultIntegrator<double>();
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver =
    Rythmos::timeStepNonlinearSolver<double>();
  const RCP<Rythmos::SolverAcceptingStepperBase<double> > stepper =
    Rythmos::backwardEulerStepper<double>(thyraModel, stepSolver);

  return rcp(new RythmosSolver<double>(integrator, stepper, stepSolver, thyraModel, finalTime));
}

const RCP<RythmosSolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
    double finalTime,
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &steadyStateModel)
{
  const RCP<Rythmos::DefaultIntegrator<double> > integrator =
    Rythmos::defaultIntegrator<double>();
  const RCP<Thyra::NonlinearSolverBase<double> > stepSolver =
    Rythmos::timeStepNonlinearSolver<double>();
  const RCP<Rythmos::SolverAcceptingStepperBase<double> > stepper =
    Rythmos::backwardEulerStepper<double>(thyraModel, stepSolver);

  return rcp(new RythmosSolver<double>(integrator, stepper, stepSolver, thyraModel, finalTime, steadyStateModel));
}

const RCP<RythmosSolver<double> > solverNew(
    const RCP<EpetraExt::ModelEvaluator> &epetraModel,
    double finalTime)
{
  return solverNew(thyraModelNew(epetraModel), finalTime);
}

void vectorAssign(Thyra::VectorBase<double> &target, const ArrayView<const double> &values) {
  for (int i = 0; i < values.size(); ++i) {
    Thyra::set_ele(i, values[i], ptrFromRef(target));
  }
}

// Floating point tolerance
const double tol = 1.0e-8;

TEUCHOS_UNIT_TEST(Piro_RythmosSolver, TimeZero_Solution)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const double finalTime = 0.0;

  const RCP<RythmosSolver<double> > solver = solverNew(model, finalTime);

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

TEUCHOS_UNIT_TEST(Piro_RythmosSolver, TimeZero_Response)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();

  const int responseIndex = 0;

  const RCP<Thyra::VectorBase<double> > expectedResponse =
    Thyra::createMember(model->get_g_space(responseIndex));
  {
    const Thyra::MEB::InArgs<double> modelInArgs = model->getNominalValues();
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_g(responseIndex, expectedResponse);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<RythmosSolver<double> > solver = solverNew(model, finalTime);

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

TEUCHOS_UNIT_TEST(Piro_RythmosSolver, TimeZero_DefaultSolutionSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const double finalTime = 0.0;

  const RCP<RythmosSolver<double> > solver = solverNew(model, finalTime);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.0, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 0.0, 0.0, 0.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
  TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromVector(*dxdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_RythmosSolver, TimeZero_DefaultResponseSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();

  const int responseIndex = 0;
  const int parameterIndex = 0;

  const Thyra::MEB::Derivative<double> dgdp_deriv_expected =
    Thyra::create_DgDp_mv(*model, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp_expected = dgdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> modelInArgs = model->getNominalValues();
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv_expected);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<RythmosSolver<double> > solver = solverNew(model, finalTime);

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

TEUCHOS_UNIT_TEST(Piro_RythmosSolver, TimeZero_ResponseAndDefaultSensitivities)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();

  const int responseIndex = 0;
  const int parameterIndex = 0;

  const RCP<Thyra::VectorBase<double> > expectedResponse =
    Thyra::createMember(model->get_g_space(responseIndex));
  {
    const Thyra::MEB::InArgs<double> modelInArgs = model->getNominalValues();
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_g(responseIndex, expectedResponse);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const Thyra::MEB::Derivative<double> dgdp_deriv_expected =
    Thyra::create_DgDp_mv(*model, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp_expected = dgdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> modelInArgs = model->getNominalValues();
    Thyra::MEB::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv_expected);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<RythmosSolver<double> > solver = solverNew(model, finalTime);

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

#ifdef Piro_ENABLE_NOX
TEUCHOS_UNIT_TEST(Piro_RythmosSolver, SteadyState_SolutionSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > steadyStateSolver(
      new NOXSolver<double>(rcp(new ParameterList), model));

  const int solutionResponseIndex = steadyStateSolver->Ng() - 1;
  const int parameterIndex = 0;

  const Thyra::MEB::Derivative<double> dxdp_deriv_expected =
    Thyra::create_DgDp_mv(*steadyStateSolver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dxdp_expected = dxdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> steadyInArgs = steadyStateSolver->getNominalValues();
    Thyra::MEB::OutArgs<double> steadyOutArgs = steadyStateSolver->createOutArgs();
    steadyOutArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv_expected);
    steadyStateSolver->evalModel(steadyInArgs, steadyOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<RythmosSolver<double> > solver = solverNew(model, finalTime, steadyStateSolver);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  solver->evalModel(inArgs, outArgs);

  TEST_EQUALITY(dxdp->domain()->dim(), dxdp_expected->domain()->dim());
  TEST_EQUALITY(dxdp->range()->dim(), dxdp_expected->range()->dim());
  for (int i = 0; i < dxdp_expected->domain()->dim(); ++i) {
    const Array<double> actual = arrayFromVector(*dxdp->col(i));
    const Array<double> expected = arrayFromVector(*dxdp_expected->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_RythmosSolver, SteadyState_ResponseSensitivity)
{
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > model = defaultModelNew();
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > steadyStateSolver(
      new NOXSolver<double>(rcp(new ParameterList), model));

  const int responseIndex = 0;
  const int parameterIndex = 0;

  const Thyra::MEB::Derivative<double> dgdp_deriv_expected =
    Thyra::create_DgDp_mv(*steadyStateSolver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<const Thyra::MultiVectorBase<double> > dgdp_expected = dgdp_deriv_expected.getMultiVector();
  {
    const Thyra::MEB::InArgs<double> steadyInArgs = steadyStateSolver->getNominalValues();
    Thyra::MEB::OutArgs<double> steadyOutArgs = steadyStateSolver->createOutArgs();
    steadyOutArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv_expected);
    steadyStateSolver->evalModel(steadyInArgs, steadyOutArgs);
  }

  const double finalTime = 0.0;
  const RCP<RythmosSolver<double> > solver = solverNew(model, finalTime, steadyStateSolver);

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
#endif /* Piro_ENABLE_NOX */

#endif /* Piro_ENABLE_Rythmos */
