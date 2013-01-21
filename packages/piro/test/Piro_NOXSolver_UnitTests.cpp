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

#ifdef Piro_ENABLE_NOX
#include "Piro_NOXSolver.hpp"

#include "MockModelEval_A.hpp"

#include "Piro_Test_WeakenedModelEvaluator.hpp"
#include "Piro_Test_ThyraSupport.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

#include "Teuchos_UnitTestHarness.hpp"

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

const RCP<NOXSolver<double> > solverNew(const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel)
{
  const RCP<ParameterList> piroParams(new ParameterList("Piro Parameters"));
  return rcp(new NOXSolver<double>(piroParams, thyraModel));
}

const RCP<NOXSolver<double> > solverNew(const RCP<EpetraExt::ModelEvaluator> &epetraModel)
{
  return solverNew(thyraModelNew(epetraModel));
}

// Floating point tolerance
const double tol = 1.0e-8;

// Tests

TEUCHOS_UNIT_TEST(Piro_NOXSolver, Spaces)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  TEST_ASSERT(solver->Np() == 1);
  TEST_ASSERT(solver->Ng() == 2);

  const int parameterIndex = 0;
  const int responseIndex = 0;
  const int solutionResponseIndex = solver->Ng() - 1;
  TEST_ASSERT(nonnull(solver->get_p_space(parameterIndex)));
  TEST_ASSERT(nonnull(solver->get_g_space(responseIndex)));
  TEST_ASSERT(nonnull(solver->get_g_space(solutionResponseIndex)));

  // TODO
  //TEST_THROW(solver->get_x_space(), std::exception);
  //TEST_THROW(solver->get_f_space(), std::exception);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, Solution)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  outArgs.set_g(solutionResponseIndex, Thyra::createMember(*solver->get_g_space(solutionResponseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(solutionResponseIndex));
  const Array<double> expected = tuple(1.0, 2.0, 3.0, 4.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SolutionForMissingParameterValues)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->createInArgs();
  const int parameterIndex = 0;
  inArgs.set_p(parameterIndex, null);

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  outArgs.set_g(solutionResponseIndex, Thyra::createMember(*solver->get_g_space(solutionResponseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(solutionResponseIndex));
  const Array<double> expected = tuple(1.0, 2.0, 3.0, 4.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SolutionForAlternateParameterValues)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  {
    const int parameterIndex = 0;
    const RCP<Thyra::VectorBase<double> > p_in = Thyra::createMember(*solver->get_p_space(0));
    TEST_EQUALITY(p_in->space()->dim(), 2);
    Thyra::set_ele(0, 1.0, p_in.ptr());
    Thyra::set_ele(1, 0.0, p_in.ptr());
    inArgs.set_p(parameterIndex, p_in);
  }

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  outArgs.set_g(solutionResponseIndex, Thyra::createMember(*solver->get_g_space(solutionResponseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(solutionResponseIndex));
  const Array<double> expected = tuple(1.0, 1.0, 2.0, 3.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, Response)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  outArgs.set_g(responseIndex, Thyra::createMember(*solver->get_g_space(responseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(responseIndex));
  const Array<double> expected = tuple(8.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, ResponseForMissingParameterValues)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->createInArgs();
  const int parameterIndex = 0;
  inArgs.set_p(parameterIndex, null);

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  outArgs.set_g(responseIndex, Thyra::createMember(*solver->get_g_space(responseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(responseIndex));
  const Array<double> expected = tuple(8.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, ResponseForAlternateParameterValues)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  {
    const int parameterIndex = 0;
    const RCP<Thyra::VectorBase<double> > p_in = Thyra::createMember(*solver->get_p_space(0));
    TEST_EQUALITY(p_in->space()->dim(), 2);
    Thyra::set_ele(0, 1.0, p_in.ptr());
    Thyra::set_ele(1, 0.0, p_in.ptr());
    inArgs.set_p(parameterIndex, p_in);
  }

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  outArgs.set_g(responseIndex, Thyra::createMember(*solver->get_g_space(responseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(responseIndex));
  const Array<double> expected = tuple(18.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SolutionSensitivityMvJac)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

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
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
  TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromVector(*dxdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SolutionSensitivityOp)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dxdp =
    solver->create_DgDp_op(solutionResponseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dxdp));

  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp);

  solver->evalModel(inArgs, outArgs);

  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromLinOp(*dxdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SolutionSensitivityOp_NoDfDpMv)
{
  // Disable support for MultiVector-based DfDp derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDfDpMv(thyraModelNew(epetraModelNew())));
  const RCP<NOXSolver<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dxdp =
    solver->create_DgDp_op(solutionResponseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dxdp));

  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp);

  solver->evalModel(inArgs, outArgs);

  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromLinOp(*dxdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvJac)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityAndResponseSensitivityMvJac)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int parameterIndex = 0;

  const int solutionResponseIndex = solver->Ng() - 1;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  const int responseIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Solution sensitivity
  {
    const Array<Array<double> > expected = tuple(
        Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
        Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
    TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
      const Array<double> actual = arrayFromVector(*dxdp->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
    }
  }

  // Response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      const Array<double> actual = arrayFromVector(*dgdp->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
    }
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvGrad)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  const Array<double> actual = arrayFromVector(*dgdp->col(parameterIndex));
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityOp)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dgdp =
    solver->create_DgDp_op(responseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dgdp));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromLinOp(*dgdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvJac_NoDgDxMv)
{
  // Disable support for MultiVector-based DgDx derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDxMv(thyraModelNew(epetraModelNew())));
  const RCP<NOXSolver<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvGrad_NoDgDpMvJac)
{
  // Disable support for Jacobian-oriented MultiVector DgDx derivative
  // (Only gradient layout is available)
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDpMvJac(thyraModelNew(epetraModelNew())));
  const RCP<NOXSolver<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  const Array<double> actual = arrayFromVector(*dgdp);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityOp_NoDgDpMv)
{
  // Disable support for MultiVector-based DgDp derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDpMv(thyraModelNew(epetraModelNew())));
  const RCP<NOXSolver<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dgdp =
    solver->create_DgDp_op(responseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dgdp));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromLinOp(*dgdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvGradWithSolutionSensitivityOp)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int parameterIndex = 0;

  // Request solution sensitivity (LINOP layout)
  const int solutionResponseIndex = solver->Ng() - 1;
  const RCP<Thyra::LinearOpBase<double> > dxdp =
    solver->create_DgDp_op(solutionResponseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dxdp));
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp);

  // Request response sensitivity (MV_GRAD layout)
  const int responseIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Verify solution sensitivity
  {
    const Array<Array<double> > expected = tuple(
        Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
        Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
    TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
      const Array<double> actual = arrayFromLinOp(*dxdp, i);
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
    }
  }

  // Verifiy response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    const Array<double> actual = arrayFromVector(*dgdp);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityOpWithSolutionSensitivityMvJac)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int parameterIndex = 0;

  // Request solution sensitivity (MV_JAC layout)
  const int solutionResponseIndex = solver->Ng() - 1;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  // Request response sensitivity (LINOP layout)
  const int responseIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dgdp =
    solver->create_DgDp_op(responseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dgdp));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Verify solution sensitivity
  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
  TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromVector(*dxdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }

  // Verify response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    for (int i = 0; i < expected.size(); ++i) {
      const Array<double> actual = arrayFromLinOp(*dgdp, i);
      TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
    }
  }

}

#endif /*Piro_ENABLE_NOX*/
