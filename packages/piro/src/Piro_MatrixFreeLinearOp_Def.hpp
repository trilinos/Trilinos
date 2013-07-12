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

#include "Piro_MatrixFreeLinearOp.hpp"

#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

#include "Teuchos_TestForException.hpp"
#include "Teuchos_ScalarTraits.hpp"

template <typename Scalar>
Piro::MatrixFreeLinearOp<Scalar>::MatrixFreeLinearOp(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model) :
  model_(model),
  basePoint_(),
  f_base_(Teuchos::null)
{
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::MatrixFreeLinearOp<Scalar>::range() const
{
  return model_->get_f_space();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::MatrixFreeLinearOp<Scalar>::domain() const
{
  return model_->get_x_space();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
Piro::MatrixFreeLinearOp<Scalar>::model() const
{
  return model_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
Piro::MatrixFreeLinearOp<Scalar>::f_base() const
{
  return f_base_;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::MatrixFreeLinearOp<Scalar>::basePoint() const
{
  return basePoint_;
}

template <typename Scalar>
void
Piro::MatrixFreeLinearOp<Scalar>::setBase(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &f_base)
{
  // Shallow copies of thr input parameters for the base point
  basePoint_ = basePoint;
  // Shallow copy of the precomputed residual at the base point
  f_base_ = f_base;
}

template <typename Scalar>
bool
Piro::MatrixFreeLinearOp<Scalar>::opSupportedImpl(Thyra::EOpTransp M_trans) const
{
  return (M_trans == Thyra::NOTRANS) && (Teuchos::nonnull(f_base_));
}

template <typename Scalar>
void
Piro::MatrixFreeLinearOp<Scalar>::applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar> &X,
    const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta) const
{
  using Teuchos::RCP;
  using Teuchos::Ptr;

  TEUCHOS_TEST_FOR_EXCEPTION(
      !this->opSupported(M_trans),
      Thyra::Exceptions::OpNotSupported,
      this->description() << " does not support operation " << Thyra::toString(M_trans));

  TEUCHOS_TEST_FOR_EXCEPTION(
      !X.range()->isCompatible(*this->domain()),
      Thyra::Exceptions::IncompatibleVectorSpaces,
      "Domain of " << this->description() << ": " << this->domain()->description() <<
      " is not compatible with column space of " << X.description() << ": " << X.range()->description());

  TEUCHOS_TEST_FOR_EXCEPTION(
      !Y->range()->isCompatible(*this->range()),
      Thyra::Exceptions::IncompatibleVectorSpaces,
      "Range of " << this->description() << ": " << this->range()->description() <<
      " is not compatible with column space of " << Y->description() << ": " << Y->range()->description());

  TEUCHOS_TEST_FOR_EXCEPTION(
      !Y->domain()->isCompatible(*X.domain()),
      Thyra::Exceptions::IncompatibleVectorSpaces,
      "Row space of " << Y->description() << ": " << Y->domain()->description() <<
      " is not compatible with row space of " << X.description() << ": " << X.domain()->description());

  TEUCHOS_TEST_FOR_EXCEPTION(
      &X == Y.get(),
      std::logic_error,
      "X and Y arguments are both aliases of " << X.description());

  if (alpha == Teuchos::ScalarTraits<Scalar>::zero()) {
    // Y <- beta * Y
    Thyra::Vt_S(Y, beta);
    return;
  }

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMagnitude;

  const RCP<const Thyra::VectorBase<Scalar> > x_dot_base = basePoint_.get_x_dot();

  RCP<const Thyra::VectorBase<Scalar> > x_base = basePoint_.get_x();
  if (Teuchos::is_null(x_base)) {
    x_base = model_->getNominalValues().get_x();
  }
  x_base.assert_not_null();

  const ScalarMagnitude norm_x_base = Thyra::norm_2(*x_base);

  // Number of columns common to both vectors X and Y
  // (X and Y have compatible row spaces)
  const Thyra::Ordinal colCount = X.domain()->dim();
  for (Teuchos::Ordinal j = Teuchos::Ordinal(); j < colCount; ++j) {
    const RCP<const Thyra::VectorBase<Scalar> > X_vec = X.col(j);
    const RCP<Thyra::VectorBase<Scalar> > Y_vec = Y->col(j);

    const ScalarMagnitude norm_dx = Thyra::norm_2(*X_vec);

    if (norm_dx == Teuchos::ScalarTraits<ScalarMagnitude>::zero()) {
      if (beta == Teuchos::ScalarTraits<Scalar>::zero()) {
        // Y_vec <- 0
        Thyra::put_scalar(Teuchos::ScalarTraits<ScalarMagnitude>::zero(), Y_vec.ptr());
      } else {
        // Y_vec <- beta * Y_vec
        Thyra::scale(beta, Y_vec.ptr());
      }
    } else {
      // Scalar perturbation
      const ScalarMagnitude relative_pert_ratio = static_cast<ScalarMagnitude>(1.0e-6);
      const ScalarMagnitude eta = (relative_pert_ratio * ((norm_x_base / norm_dx) + relative_pert_ratio));

      // Compute perturbed residual
      // Dynamic: f_pert <- f(x_dot_base + eta * (W_alpha * X), x_base + eta * (W_beta * X))
      // Static: f_pert <- f(x_base + eta * X)
      const RCP<Thyra::VectorBase<Scalar> > f_pert = Thyra::createMember(this->range());
      {
        Thyra::ModelEvaluatorBase::InArgs<Scalar> pertInArgs = model_->createInArgs();
        {
          pertInArgs.setArgs(basePoint_);

          const bool isDynamic = Teuchos::nonnull(x_dot_base);

          if (isDynamic) {
            const RCP<Thyra::VectorBase<Scalar> > x_dot_pert = Thyra::createMember(this->domain());
            const Scalar W_alpha = pertInArgs.get_alpha();
            Thyra::V_VpStV<Scalar>(x_dot_pert.ptr(), *x_dot_base, W_alpha * eta, *X_vec);
            pertInArgs.set_x_dot(x_dot_pert);
          }

          const RCP<Thyra::VectorBase<Scalar> > x_pert = Thyra::createMember(this->domain());
          const Scalar W_beta = isDynamic ? pertInArgs.get_beta() : Teuchos::ScalarTraits<Scalar>::one();
          Thyra::V_VpStV<Scalar>(x_pert.ptr(), *x_base, W_beta * eta, *X_vec);
          pertInArgs.set_x(x_pert);
        }

        Thyra::ModelEvaluatorBase::OutArgs<Scalar> pertOutArgs = model_->createOutArgs();
        {
          pertOutArgs.set_f(f_pert);
        }

        model_->evalModel(pertInArgs, pertOutArgs);
      }

      // Y <- alpha * (1/eta) * (f_pert - f_base) + beta * Y
      const Scalar alpha_over_eta = alpha / eta;

      if (beta == Teuchos::ScalarTraits<Scalar>::zero()) {
        // Y <- alpha * (1/eta) * (f_pert - f_base)
        Thyra::V_StVpStV<Scalar>(Y_vec.ptr(), alpha_over_eta, *f_pert, -alpha_over_eta, *f_base_);
      } else {
        // Aliasing f_pert and alpha_op_X (f_pert == alpha_op_X)
        const RCP<Thyra::VectorBase<Scalar> > alpha_op_X = f_pert;

        // alpha_op_X <- alpha * (1/eta) * (f_pert - f_base)
        Thyra::Vp_StV(alpha_op_X.ptr(), -Teuchos::ScalarTraits<Scalar>::one(), *f_base_);
        const Scalar alpha_over_eta = alpha / eta;
        Thyra::Vt_S(alpha_op_X.ptr(), alpha_over_eta);

        // Y <- alpha_op_X + beta * Y
        Thyra::Vp_V<Scalar>(Y_vec.ptr(), *alpha_op_X, beta);
      }
    }
  }
}
