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

#ifndef PIRO_STEADYSTATESOLVER_DEF_HPP
#define PIRO_STEADYSTATESOLVER_DEF_HPP

#include "Piro_SteadyStateSolver.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"

#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_FancyOStream.hpp"
#include <stdexcept>
#include <cstddef>
#include <ostream>

template <typename Scalar>
Piro::SteadyStateSolver<Scalar>::
SteadyStateSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model) :
  model_(model),
  num_p_(model->Np()),
  num_g_(model->Ng())
{}

template <typename Scalar>
Piro::SteadyStateSolver<Scalar>::
SteadyStateSolver(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    int numParameters) :
  model_(model),
  num_p_(numParameters),
  num_g_(model->Ng())
{}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::SteadyStateSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p_ || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model_->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::SteadyStateSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j > num_g_ || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if (j < num_g_) {
    return model_->get_g_space(j);
  } else {
    // j == num_g_, corresponding to the state by convention
    return model_->get_x_space();
  }
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::SteadyStateSolver<Scalar>::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> result;
  result.setModelEvalDescription(this->description());
  result.set_Np(num_p_);
  return result;
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::SteadyStateSolver<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgsImpl();
  result.setArgs(
      model_->getNominalValues(),
      /* ignoreUnsupported = */ true,
      /* cloneObjects = */ false);
  return result;
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::SteadyStateSolver<Scalar>::getLowerBounds() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgsImpl();
  result.setArgs(
      model_->getLowerBounds(),
      /* ignoreUnsupported = */ true,
      /* cloneObjects = */ false);
  return result;
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::SteadyStateSolver<Scalar>::getUpperBounds() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgsImpl();
  result.setArgs(
      model_->getUpperBounds(),
      /* ignoreUnsupported = */ true,
      /* cloneObjects = */ false);
  return result;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> Piro::SteadyStateSolver<Scalar>::createInArgs() const
{
  return this->createInArgsImpl();
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> Piro::SteadyStateSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> result;
  result.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  result.set_Np_Ng(num_p_, num_g_ + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model_->createOutArgs();

  // Sensitivity support (Forward approach only)
  // Jacobian solver required for all sensitivities
  if (modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W)) {
    for (int l = 0; l < num_p_; ++l) {
      // Solution sensitivities: DxDp(l)
      // DfDp(l) required
      const Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_support =
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
      const bool dxdp_linOpSupport =
        dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
      const bool dxdp_mvJacSupport =
        dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
      {
        Thyra::ModelEvaluatorBase::DerivativeSupport dxdp_support;
        if (dxdp_linOpSupport) {
          dxdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        }
        if (dxdp_mvJacSupport) {
          dxdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        }
        result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, num_g_, l, dxdp_support);
      }

      // Response sensitivities: DgDp(j, l)
      // DxDp(l) required
      if (dxdp_linOpSupport || dxdp_mvJacSupport) {
        for (int j = 0; j < num_g_; ++j) {
          // DgDx(j) required
          const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
            modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j);
          const bool dgdx_linOpSupport =
            dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
          const bool dgdx_mvGradSupport =
            dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
          if (dgdx_linOpSupport || dgdx_mvGradSupport) {
            // Dgdp(j, l) required
            const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
              modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
            Thyra::ModelEvaluatorBase::DerivativeSupport total_dgdp_support;
            if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP) &&
                dgdx_linOpSupport && dxdp_linOpSupport) {
              total_dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
            }
            if (dxdp_mvJacSupport) {
              if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
                total_dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
              }
              if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) &&
                  dgdx_mvGradSupport) {
                total_dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
              }
            }
            result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l, total_dgdp_support);
          }
        }
      }
    }
  }

  for (int i=0; i<num_g_; i++) {
    for (int j=0; j<num_p_; j++) {
      Thyra::ModelEvaluatorBase::DerivativeSupport ds = modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j);
      if (!ds.none()) {
        ds.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j, ds);
      }
    }
  }

  return result;
}

template <typename Scalar>
const Thyra::ModelEvaluator<Scalar> &
Piro::SteadyStateSolver<Scalar>::getModel() const
{
  return *model_;
}

template <typename Scalar>
int
Piro::SteadyStateSolver<Scalar>::num_p() const
{
  return num_p_;
}

template <typename Scalar>
int
Piro::SteadyStateSolver<Scalar>::num_g() const
{
  return num_g_;
}

template <typename Scalar>
void Piro::SteadyStateSolver<Scalar>::evalConvergedModel(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  //int g_size = 0;  // Commenting out since g_size is not used.
  // Solution at convergence is the response at index num_g_
  {
    const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g_);
    if (Teuchos::nonnull(gx_out)) {
      //g_size = gx_out->space()->dim();
      Thyra::copy(*modelInArgs.get_x(), gx_out.ptr());
    }
  }

  // Setup output for final evalution of underlying model
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model_->createOutArgs();
  {
    // Responses
    for (int j = 0; j < num_g_; ++j) {
      const RCP<Thyra::VectorBase<Scalar> > g_out = outArgs.get_g(j);
      // Forward to underlying model
      modelOutArgs.set_g(j, g_out);
    }

    // Jacobian
    {
      bool jacobianRequired = false;
      for (int j = 0; j <= num_g_; ++j) { // resize
        for (int l = 0; l < num_p_; ++l) {
          const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
            outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
          if (!dgdp_support.none()) {
            const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
              outArgs.get_DgDp(j, l);
            if (!dgdp_deriv.isEmpty()) {
              jacobianRequired = true;
            }
          }
        }
      }
      if (jacobianRequired) {
        const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
          model_->create_W();
        modelOutArgs.set_W(jacobian);
      }
    }

    // DfDp derivatives
    for (int l = 0; l < num_p_; ++l) {
      Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_request;
      for (int j = 0; j <= num_g_; ++j) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
          outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        if (!dgdp_support.none()) {
          const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
            outArgs.get_DgDp(j, l);
          if (Teuchos::nonnull(dgdp_deriv.getLinearOp())) {
            dfdp_request.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
          } else if (Teuchos::nonnull(dgdp_deriv.getMultiVector())) {
            dfdp_request.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
          }
        }
      }

      if (!dfdp_request.none()) {
        Thyra::ModelEvaluatorBase::Derivative<Scalar> dfdp_deriv;
        if (dfdp_request.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
          dfdp_deriv = Thyra::create_DfDp_mv(*model_, l, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        } else if (dfdp_request.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          dfdp_deriv = model_->create_DfDp_op(l);
        }
        modelOutArgs.set_DfDp(l, dfdp_deriv);
      }
    }

    // DgDx derivatives
    for (int j = 0; j < num_g_; ++j) {
      Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_request;
      for (int l = 0; l < num_p_; ++l) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
          outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        if (!dgdp_support.none()) {
          const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
            outArgs.get_DgDp(j, l);
          if (!dgdp_deriv.isEmpty()) {
            const bool dgdp_mvGrad_required =
              Teuchos::nonnull(dgdp_deriv.getMultiVector()) &&
              dgdp_deriv.getMultiVectorOrientation() == Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
            if (dgdp_mvGrad_required) {
              dgdx_request.plus(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
            } else {
              dgdx_request.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
            }
          }
        }
      }

      if (!dgdx_request.none()) {
        Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_deriv;
        if (dgdx_request.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)) {
          dgdx_deriv = Thyra::create_DgDx_mv(*model_, j, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        } else if (dgdx_request.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          dgdx_deriv = model_->create_DgDx_op(j);
        }
        modelOutArgs.set_DgDx(j, dgdx_deriv);
      }
    }

    // DgDp derivatives
    for (int l = 0; l < num_p_; ++l) {
      for (int j = 0; j < num_g_; ++j) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
          outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        if (!dgdp_support.none()) {
          const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
            outArgs.get_DgDp(j, l);
          Thyra::ModelEvaluatorBase::Derivative<Scalar> model_dgdp_deriv;
          const RCP<Thyra::LinearOpBase<Scalar> > dgdp_op = dgdp_deriv.getLinearOp();
          if (Teuchos::nonnull(dgdp_op)) {
            model_dgdp_deriv = model_->create_DgDp_op(j, l);
          } else {
            model_dgdp_deriv = dgdp_deriv;
          }
          if (!model_dgdp_deriv.isEmpty()) {
            modelOutArgs.set_DgDp(j, l, model_dgdp_deriv);
          }
        }
      }
    }
  }

  // Evaluate underlying model
  model_->evalModel(modelInArgs, modelOutArgs);

  // Assemble user-requested sensitivities
  if (modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W)) {
    const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
      modelOutArgs.get_W();
    if (Teuchos::nonnull(jacobian)) {
      for (int l = 0; l < num_p_; ++l) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_support =
          modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
        if (!dfdp_support.none()) {
          const Thyra::ModelEvaluatorBase::Derivative<Scalar> dfdp_deriv =
            modelOutArgs.get_DfDp(l);
          const RCP<Thyra::MultiVectorBase<Scalar> > dfdp_mv =
            dfdp_deriv.getMultiVector();
          RCP<Thyra::LinearOpBase<Scalar> > dfdp_op =
            dfdp_deriv.getLinearOp();
          if (Teuchos::is_null(dfdp_op)) {
            dfdp_op = dfdp_mv;
          }

          const Thyra::ModelEvaluatorBase::Derivative<Scalar> dxdp_deriv =
            outArgs.get_DgDp(num_g_, l);
          const RCP<Thyra::LinearOpBase<Scalar> > dxdp_op =
            dxdp_deriv.getLinearOp();
          const RCP<Thyra::MultiVectorBase<Scalar> > dxdp_mv =
            dxdp_deriv.getMultiVector();

          RCP<const Thyra::LinearOpBase<Scalar> > minus_dxdp_op;
          RCP<Thyra::MultiVectorBase<Scalar> > minus_dxdp_mv;
          if (Teuchos::nonnull(dfdp_mv)) {
            if (Teuchos::nonnull(dxdp_mv)) {
              minus_dxdp_mv = dxdp_mv; // Use user-provided object as temporary
            } else {
              minus_dxdp_mv =
                Thyra::createMembers(model_->get_x_space(), model_->get_p_space(l));
              minus_dxdp_op = minus_dxdp_mv;
            }
          }

          if (Teuchos::is_null(minus_dxdp_op) && Teuchos::nonnull(dfdp_op)) {
            const RCP<const Thyra::LinearOpBase<Scalar> > dfdx_inv_op =
              Thyra::inverse<Scalar>(jacobian);
            minus_dxdp_op = Thyra::multiply<Scalar>(dfdx_inv_op, dfdp_op);
          }

          if (Teuchos::nonnull(minus_dxdp_mv) && Teuchos::nonnull(dfdp_mv)) {
            Thyra::assign(minus_dxdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

            const Thyra::SolveCriteria<Scalar> defaultSolveCriteria;
            const Thyra::SolveStatus<Scalar> solveStatus =
              Thyra::solve(
                  *jacobian,
                  Thyra::NOTRANS,
                  *dfdp_mv,
                  minus_dxdp_mv.ptr(),
                  Teuchos::ptr(&defaultSolveCriteria));

            //  AGS: Made this a warning instead of exception since it is 'just' post-processing
            if (solveStatus.solveStatus == Thyra::SOLVE_STATUS_UNCONVERGED)
              *(Teuchos::VerboseObjectBase::getDefaultOStream() ) <<
              "\nWARNING: Linear Solver in sensitivity computation failed to fully converge\n"
              << "         Accuracy of sensitivity calculations is less then requested." << std::endl;
          }

          // Solution sensitivities
          if (Teuchos::nonnull(dxdp_mv)) {
            minus_dxdp_mv = Teuchos::null; // Invalidates temporary
            Thyra::scale(-Teuchos::ScalarTraits<Scalar>::one(), dxdp_mv.ptr());
          } else if (Teuchos::nonnull(dxdp_op)) {
            const RCP<Thyra::DefaultMultipliedLinearOp<Scalar> > dxdp_op_downcasted =
              Teuchos::rcp_dynamic_cast<Thyra::DefaultMultipliedLinearOp<Scalar> >(dxdp_op);
            TEUCHOS_TEST_FOR_EXCEPTION(
                Teuchos::is_null(dxdp_op_downcasted),
                std::invalid_argument,
                "Illegal operator for DgDp(" <<
                "j = " << num_g_ << ", " <<
                "index l = " << l << ")\n");

            const RCP<const Thyra::LinearOpBase<Scalar> > minus_id_op =
              Thyra::scale<Scalar>(-Teuchos::ScalarTraits<Scalar>::one(), Thyra::identity(dfdp_op->domain()));

            dxdp_op_downcasted->initialize(Teuchos::tuple(minus_dxdp_op, minus_id_op));
          }

          // Response sensitivities
          for (int j = 0; j < num_g_; ++j) {
            const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
              outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
            if (!dgdp_support.none()) {
              const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
                outArgs.get_DgDp(j, l);
              if (!dgdp_deriv.isEmpty()) {
                const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_deriv =
                  modelOutArgs.get_DgDx(j);
                const RCP<const Thyra::MultiVectorBase<Scalar> > dgdx_mv =
                  dgdx_deriv.getMultiVector();
                RCP<const Thyra::LinearOpBase<Scalar> > dgdx_op =
                  dgdx_deriv.getLinearOp();
                if (Teuchos::is_null(dgdx_op)) {
                  dgdx_op = Thyra::adjoint<Scalar>(dgdx_mv);
                }

                const RCP<Thyra::LinearOpBase<Scalar> > dgdp_op =
                  dgdp_deriv.getLinearOp();
                if (Teuchos::nonnull(dgdp_op)) {
                  const RCP<Thyra::DefaultAddedLinearOp<Scalar> > dgdp_op_downcasted =
                    Teuchos::rcp_dynamic_cast<Thyra::DefaultAddedLinearOp<Scalar> >(dgdp_op);
                  TEUCHOS_TEST_FOR_EXCEPTION(
                      Teuchos::is_null(dgdp_op_downcasted),
                      std::invalid_argument,
                      "Illegal operator for DgDp(" <<
                      "j = " << j << ", " <<
                      "index l = " << l << ")\n");

                  dgdp_op_downcasted->uninitialize();

                  const RCP<const Thyra::LinearOpBase<Scalar> > implicit_dgdp_op =
                    Thyra::multiply<Scalar>(
                      Thyra::scale<Scalar>(-Teuchos::ScalarTraits<Scalar>::one(), dgdx_op),
                      minus_dxdp_op);

                  const RCP<const Thyra::LinearOpBase<Scalar> > model_dgdp_op =
                    modelOutArgs.get_DgDp(j, l).getLinearOp();

                  Teuchos::Array<RCP<const Thyra::LinearOpBase<Scalar> > > op_args(2);
                  op_args[0] = model_dgdp_op;
                  op_args[1] = implicit_dgdp_op;
                  dgdp_op_downcasted->initialize(op_args);
                }

                const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv =
                  dgdp_deriv.getMultiVector();
                if (Teuchos::nonnull(dgdp_mv)) {
                  if (dgdp_deriv.getMultiVectorOrientation() == Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) {
                    if (Teuchos::nonnull(dxdp_mv)) {
                      Thyra::apply(
                          *dxdp_mv,
                          Thyra::TRANS,
                          *dgdx_mv,
                          dgdp_mv.ptr(),
                          Teuchos::ScalarTraits<Scalar>::one(),
                          Teuchos::ScalarTraits<Scalar>::one());
                    } else {
                      Thyra::apply(
                          *minus_dxdp_mv,
                          Thyra::TRANS,
                          *dgdx_mv,
                          dgdp_mv.ptr(),
                          -Teuchos::ScalarTraits<Scalar>::one(),
                          Teuchos::ScalarTraits<Scalar>::one());
                    }
                  } else {
                    if (Teuchos::nonnull(dxdp_mv)) {
                      Thyra::apply(
                          *dgdx_op,
                          Thyra::NOTRANS,
                          *dxdp_mv,
                          dgdp_mv.ptr(),
                          Teuchos::ScalarTraits<Scalar>::one(),
                          Teuchos::ScalarTraits<Scalar>::one());
                    } else {
                      Thyra::apply(
                          *dgdx_op,
                          Thyra::NOTRANS,
                          *minus_dxdp_mv,
                          dgdp_mv.ptr(),
                          -Teuchos::ScalarTraits<Scalar>::one(),
                          Teuchos::ScalarTraits<Scalar>::one());
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::SteadyStateSolver<Scalar>::create_DgDp_op_impl(int j, int l) const
{
  const Teuchos::Array<Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > > dummy =
    Teuchos::tuple(Thyra::zero<Scalar>(this->get_g_space(j), this->get_p_space(l)));
  if (j == num_g_)  {
    return Thyra::defaultMultipliedLinearOp<Scalar>(dummy);
  } else {
    return Teuchos::rcp(new Thyra::DefaultAddedLinearOp<Scalar>(dummy));
  }
}
#endif /*PIRO_STEADYSTATESOLVER_DEF_HPP*/
