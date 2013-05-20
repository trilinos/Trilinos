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

#include "Piro_MatrixFreeDecorator.hpp"

#include "Piro_MatrixFreeLinearOp.hpp"

#include "Piro_Test_ThyraSupport.hpp"

#include "Piro_Test_WeakenedModelEvaluator.hpp"

#include "MockModelEval_A.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

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

// Evaluation support

void eval_f_default(const Thyra::ModelEvaluator<double> &model, Thyra::VectorBase<double> &residual)
{
  const Thyra::MEB::InArgs<double> inArgs = model.createInArgs();
  Thyra::MEB::OutArgs<double> outArgs = model.createOutArgs();
  outArgs.set_f(rcpFromRef(residual));
  model.evalModel(inArgs, outArgs);
}

void eval_W_op_default(const Thyra::ModelEvaluator<double> &model, Thyra::LinearOpBase<double> &jacobian)
{
  const Thyra::MEB::InArgs<double> inArgs = model.createInArgs();
  Thyra::MEB::OutArgs<double> outArgs = model.createOutArgs();
  outArgs.set_W_op(rcpFromRef(jacobian));
  model.evalModel(inArgs, outArgs);
}

void eval_f_W_op_default(
    const Thyra::ModelEvaluator<double> &model,
    const RCP<Thyra::VectorBase<double> > &residual,
    const RCP<Thyra::LinearOpBase<double> > &jacobian)
{
  const Thyra::MEB::InArgs<double> inArgs = model.createInArgs();
  Thyra::MEB::OutArgs<double> outArgs = model.createOutArgs();
  outArgs.set_f(residual);
  outArgs.set_W_op(jacobian);
  model.evalModel(inArgs, outArgs);
}

void eval_W_op_dynamic(
    const Thyra::ModelEvaluator<double> &model,
    double alpha,
    const Thyra::VectorBase<double> &x_dot,
    double beta,
    const Thyra::VectorBase<double> &x,
    Thyra::LinearOpBase<double> &jacobian)
{
  Thyra::MEB::InArgs<double> inArgs = model.createInArgs();
  inArgs.set_x_dot(rcpFromRef(x_dot));
  inArgs.set_x(rcpFromRef(x));
  inArgs.set_alpha(alpha);
  inArgs.set_beta(beta);

  Thyra::MEB::OutArgs<double> outArgs = model.createOutArgs();
  outArgs.set_W_op(rcpFromRef(jacobian));

  model.evalModel(inArgs, outArgs);
}

// Checking support

#include "Thyra_MultiVectorStdOps.hpp"

template <typename Scalar>
RCP<Thyra::MultiVectorBase<Scalar> >
identityMultiVectorNew(const RCP<const Thyra::VectorSpaceBase<Scalar> > &space)
{
  const RCP<Thyra::MultiVectorBase<Scalar> > result = Thyra::createMembers(space, space);
  Thyra::assign(result.ptr(), ScalarTraits<Scalar>::zero());

  const Thyra::Ordinal colCount = result->range()->dim();
  for (Thyra::Ordinal j = Thyra::Ordinal(); j < colCount; ++j) {
    const RCP<Thyra::VectorBase<Scalar> > v = result->col(j);
    Thyra::set_ele(j, ScalarTraits<Scalar>::one(), v.ptr());
  }

  return result;
}

// Floating point tolerances
const double tightTol = 1.0e-8;
const double relaxedTol = 1.0e-5;

// Tests: Operator

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, JacobianOperatorInitialization)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel(new WeakenedModelEvaluator_NoW(model));
  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(weakenedModel));
  const RCP<Thyra::LinearOpBase<double> > jacobian = decorator->create_W_op();

  const RCP<Thyra::VectorBase<double> > rhs = Thyra::createMember(jacobian->domain());
  Thyra::assign(rhs.ptr(), 1.0);
  const RCP<Thyra::VectorBase<double> > result = Thyra::createMember(jacobian->range());
  Thyra::assign(result.ptr(), 0.0);

  TEST_ASSERT(isPartiallyInitialized(*jacobian));
  TEST_ASSERT(!opSupported(*jacobian, Thyra::NOTRANS));
  TEST_THROW(Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr()), Thyra::Exceptions::OpNotSupported);

  eval_W_op_default(*decorator, *jacobian);

  TEST_ASSERT(isFullyInitialized(*jacobian));
  TEST_ASSERT(opSupported(*jacobian, Thyra::NOTRANS));
  TEST_NOTHROW(Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr()));
}

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, JacobianOperatorAppliesWithCoefs)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel(new WeakenedModelEvaluator_NoW(model));
  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(weakenedModel));
  const RCP<Thyra::LinearOpBase<double> > jacobian = decorator->create_W_op();
  eval_W_op_default(*decorator, *jacobian);

  {
    const RCP<Thyra::VectorBase<double> > rhs = Thyra::createMember(jacobian->domain());
    const RCP<Thyra::VectorBase<double> > result = Thyra::createMember(jacobian->range());

    const double rhs_coef = 0.0;
    const double result_coef = 0.0;
    Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr(), rhs_coef, result_coef);

    const Array<double> actual = arrayFromVector(*result);
    const Array<double> expected = tuple(0.0, 0.0, 0.0, 0.0);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tightTol);
  }

  {
    const RCP<Thyra::VectorBase<double> > rhs = Thyra::createMember(jacobian->domain());
    const RCP<Thyra::VectorBase<double> > result = Thyra::createMember(jacobian->range());

    const double rhs_coef = 0.0;
    const double result_coef = -2.0;
    Thyra::assign(result.ptr(), 1.0);

    const Array<double> expected = tuple(-2.0, -2.0, -2.0, -2.0);

    {
      Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr(), rhs_coef, result_coef);
      const Array<double> actual = arrayFromVector(*result);
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tightTol);
    }

    Thyra::assign(rhs.ptr(), 1.0);
    Thyra::assign(result.ptr(), 1.0);

    {
      Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr(), rhs_coef, result_coef);
      const Array<double> actual = arrayFromVector(*result);
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tightTol);
    }
  }

  {
    const RCP<Thyra::VectorBase<double> > rhs = Thyra::createMember(jacobian->domain());
    const RCP<Thyra::VectorBase<double> > result = Thyra::createMember(jacobian->range());

    const double rhs_coef = 2.0;
    const double result_coef = 0.0;
    Thyra::assign(rhs.ptr(), 1.0);

    const Array<double> expected = tuple(12.0, 12.0, 12.0, 12.0);

    {
      Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr(), rhs_coef, result_coef);
      const Array<double> actual = arrayFromVector(*result);
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected, relaxedTol);
    }

    Thyra::assign(rhs.ptr(), 1.0);
    Thyra::assign(result.ptr(), -1.0);

    {
      Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr(), rhs_coef, result_coef);
      const Array<double> actual = arrayFromVector(*result);
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected, relaxedTol);
    }
  }

  {
    const RCP<Thyra::VectorBase<double> > rhs = Thyra::createMember(jacobian->domain());
    const RCP<Thyra::VectorBase<double> > result = Thyra::createMember(jacobian->range());

    const double rhs_coef = 2.0;
    const double result_coef = -1.0;
    Thyra::assign(rhs.ptr(), 1.0);
    Thyra::assign(result.ptr(), -1.0);

    const Array<double> expected = tuple(13.0, 13.0, 13.0, 13.0);

    {
      Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, result.ptr(), rhs_coef, result_coef);

      const Array<double> actual = arrayFromVector(*result);
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected, relaxedTol);
    }
  }
}

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, JacobianOperatorAppliesToMv)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel(new WeakenedModelEvaluator_NoW(model));
  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(weakenedModel));
  const RCP<Thyra::LinearOpBase<double> > jacobian = decorator->create_W_op();
  eval_W_op_default(*decorator, *jacobian);

  const RCP<Thyra::MultiVectorBase<double> > mv =
    Thyra::createMembers(jacobian->range(), jacobian->domain());
  {
    const RCP<const Thyra::MultiVectorBase<double> > rhs =
      identityMultiVectorNew(jacobian->domain());
    Thyra::apply(*jacobian, Thyra::NOTRANS, *rhs, mv.ptr());
  }

  for (int i = 0; i < jacobian->domain()->dim(); ++i) {
    const Array<double> actual = arrayFromVector(*mv->col(i));
    const Array<double> expected = arrayFromLinOp(*jacobian, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tightTol);
  }
}

// Tests: Decorator

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, DefaultResidual)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::VectorBase<double> > expectedResidual =
    Thyra::createMember(model->get_f_space());
  eval_f_default(*model, *expectedResidual);

  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(model));
  const RCP<Thyra::VectorBase<double> > residual =
    Thyra::createMember(decorator->get_f_space());
  eval_f_default(*decorator, *residual);

  const Array<double> expected = arrayFromVector(*expectedResidual);
  const Array<double> actual = arrayFromVector(*residual);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tightTol);
}

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, DefaultJacobian)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::LinearOpBase<double> > expectedJacobian = model->create_W_op();
  eval_W_op_default(*model, *expectedJacobian);

  const RCP<Thyra::ModelEvaluator<double> > weakenedModel(new WeakenedModelEvaluator_NoW(model));
  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(weakenedModel));
  const RCP<Thyra::LinearOpBase<double> > jacobian = decorator->create_W_op();
  TEST_ASSERT(nonnull(jacobian));
  eval_W_op_default(*decorator, *jacobian);

  TEST_ASSERT(jacobian->domain()->isCompatible(*expectedJacobian->domain()));
  TEST_ASSERT(jacobian->range()->isCompatible(*expectedJacobian->range()));
  for (int i = 0; i < jacobian->domain()->dim(); ++i) {
    const Array<double> actual = arrayFromLinOp(*jacobian, i);
    const Array<double> expected = arrayFromLinOp(*expectedJacobian, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, relaxedTol);
  }
}

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, DefaultJacobian_ModelSupportsW)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::LinearOpBase<double> > expectedJacobian = model->create_W_op();
  eval_W_op_default(*model, *expectedJacobian);

  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(model));
  const RCP<Thyra::LinearOpBase<double> > jacobian = decorator->create_W_op();
  TEST_ASSERT(nonnull(jacobian));
  eval_W_op_default(*decorator, *jacobian);

  TEST_ASSERT(jacobian->domain()->isCompatible(*expectedJacobian->domain()));
  TEST_ASSERT(jacobian->range()->isCompatible(*expectedJacobian->range()));
  for (int i = 0; i < jacobian->domain()->dim(); ++i) {
    const Array<double> actual = arrayFromLinOp(*jacobian, i);
    const Array<double> expected = arrayFromLinOp(*expectedJacobian, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, relaxedTol);
  }
}

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, DynamicJacobian)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::LinearOpBase<double> > expectedJacobian = model->create_W_op();

  const RCP<Thyra::ModelEvaluator<double> > weakenedModel(new WeakenedModelEvaluator_NoW(model));
  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(weakenedModel));
  const RCP<Thyra::LinearOpBase<double> > jacobian = decorator->create_W_op();

  const Array<double> alphas = tuple(0.0, -1.0, 0.1, 10.0, 1.0e-6);
  const Array<double> betas = tuple(1.0, -1.0, 2.0, 10.0, 1.0e-6);

  for (Array<double>::const_iterator alphaIt = alphas.begin(); alphaIt != alphas.end(); ++alphaIt) {
    for (Array<double>::const_iterator betaIt = betas.begin(); betaIt != betas.end(); ++betaIt) {
      const double alpha = *alphaIt;
      const double beta = *betaIt;

      {
        const RCP<const Thyra::VectorBase<double> > model_x_dot = model->getNominalValues().get_x_dot();
        const RCP<const Thyra::VectorBase<double> > model_x = model->getNominalValues().get_x();
        eval_W_op_dynamic(*model, alpha, *model_x_dot, beta, *model_x, *expectedJacobian);
      }

      {
        const RCP<const Thyra::VectorBase<double> > x_dot = decorator->getNominalValues().get_x_dot();
        const RCP<const Thyra::VectorBase<double> > x = decorator->getNominalValues().get_x();
        eval_W_op_dynamic(*decorator, alpha, *x_dot, beta, *x, *jacobian);
      }

      for (int i = 0; i < jacobian->domain()->dim(); ++i) {
        const Array<double> actual = arrayFromLinOp(*jacobian, i);
        const Array<double> expected = arrayFromLinOp(*expectedJacobian, i);
        TEST_COMPARE_FLOATING_ARRAYS(actual, expected, 6.0 * relaxedTol);
      }
    }
  }
}

TEUCHOS_UNIT_TEST(Piro_MatrixFreeDecorator, ResidualAliasing)
{
  const RCP<Thyra::ModelEvaluator<double> > model = defaultModelNew();
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel(new WeakenedModelEvaluator_NoW(model));
  const RCP<Thyra::ModelEvaluator<double> > decorator(new MatrixFreeDecorator<double>(weakenedModel));

  const RCP<Thyra::LinearOpBase<double> > jacobian = decorator->create_W_op();

  {
    const RCP<Thyra::VectorBase<double> > residual =
      Thyra::createMember(decorator->get_f_space());

    eval_f_W_op_default(*decorator, residual, jacobian);

    TEST_EQUALITY(residual.total_count(), 1);
  }

  {
    const RCP<MatrixFreeLinearOp<double> > jacobian_downcasted =
      rcp_dynamic_cast<MatrixFreeLinearOp<double> >(jacobian);
    TEST_ASSERT(nonnull(jacobian_downcasted));

    const RCP<const Thyra::VectorBase<double> > f_base = jacobian_downcasted->f_base();
    TEST_ASSERT(f_base.is_valid_ptr());
    TEST_ASSERT(f_base.has_ownership());
    TEST_EQUALITY(f_base.strength(), RCP_STRONG);
  }
}
