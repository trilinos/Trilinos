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

#include "Piro_TransientSolver.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"

#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"

#define DEBUG_OUTPUT

#include <string>
#include <stdexcept>
#include <iostream>

template <typename Scalar>
Piro::TransientSolver<Scalar>::TransientSolver(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model, 
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel):
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  model_(model), 
  initialConditionModel_(icModel),
  num_p_(model->Np()), 
  num_g_(model->Ng())
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
}

template <typename Scalar>
Piro::TransientSolver<Scalar>::TransientSolver(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model, int numParameters, 
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel) : 
    out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    model_(model),
    initialConditionModel_(icModel),
    num_p_(numParameters),
    num_g_(model->Ng())
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TransientSolver<Scalar>::get_p_space(int l) const
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_p_ || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TransientSolver::get_p_map():  " <<
      "Invalid parameter index l = " <<
      l << "\n");

  return model_->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TransientSolver<Scalar>::get_g_space(int j) const
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > num_g_ || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TransientSolver::get_g_map():  " <<
      "Invalid response index j = " <<
      j << "\n");

  if (j < num_g_) {
    return model_->get_g_space(j);
  } else {
    // j == num_g_
    return model_->get_x_space();
  }
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TransientSolver<Scalar>::getNominalValues() const
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model_->getNominalValues();
  for (int l = 0; l < num_p_; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TransientSolver<Scalar>::createInArgs() const
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p_);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::TransientSolver<Scalar>::createOutArgsImpl() const
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(num_p_, num_g_ + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model_->createOutArgs();

  if (num_p_ > 0) {
    // Only one parameter supported
    const int l = 0;

    if (Teuchos::nonnull(initialConditionModel_)) {
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs =
        initialConditionModel_->createOutArgs();
      const Thyra::ModelEvaluatorBase::DerivativeSupport init_dxdp_support =
        initCondOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, initCondOutArgs.Ng() - 1, l);
      if (!init_dxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
        // Ok to return early since only one parameter supported
        return outArgs;
      }
    }

    // Computing the DxDp sensitivity for a transient problem currently requires the evaluation of
    // the mutilivector-based, Jacobian-oriented DfDp derivatives of the underlying transient model.
    const Thyra::ModelEvaluatorBase::DerivativeSupport model_dfdp_support =
      modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
    if (!model_dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
      // Ok to return early since only one parameter supported
      return outArgs;
    }

    // Solution sensitivity
    outArgs.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
        num_g_,
        l,
        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

    if (num_g_ > 0) {
      // Only one response supported
      const int j = 0;

      const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdx_support =
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j);
      if (!model_dgdx_support.none()) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdp_support =
          modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        // Response sensitivity
        Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support;
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        }
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        }
        outArgs.setSupports(
            Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
            j,
            l,
            dgdp_support);
      }
    }
  }

  return outArgs;
}


template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::TransientSolver<Scalar>::create_DgDp_op_impl(int j, int l) const
{
  TEUCHOS_ASSERT(j != num_g_);
  const Teuchos::Array<Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > > dummy =
    Teuchos::tuple(Thyra::zero<Scalar>(this->get_g_space(j), this->get_p_space(l)));
  return Teuchos::rcp(new Thyra::DefaultAddedLinearOp<Scalar>(dummy));
}

template <typename Scalar>
void 
Piro::TransientSolver<Scalar>::evalConvergedModel(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
#ifdef DEBUG_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //IKT FIXME: FILL IN! 
}
