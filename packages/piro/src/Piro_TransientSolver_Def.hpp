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

#include <string>
#include <stdexcept>
#include <iostream>

//#define DEBUG_OUTPUT

template <typename Scalar>
Piro::TransientSolver<Scalar>::TransientSolver(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model) :  
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  model_(model), 
  num_p_(model->Np()), 
  num_g_(model->Ng()),
  sensitivityMethod_(NONE)
{
  //Nothing to do
}

template <typename Scalar>
Piro::TransientSolver<Scalar>::TransientSolver(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model, int numParameters) :  
    out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
    model_(model),
    num_p_(numParameters),
    num_g_(model->Ng()),
    sensitivityMethod_(NONE) 
{
  //Nothing to do
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TransientSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_p_ || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TransientSolver::get_p_space():  " <<
      "Invalid parameter index l = " <<
      l << "\n");

  return model_->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TransientSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > num_g_ || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TransientSolver::get_g_space():  " <<
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
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p_);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::TransientSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(num_p_, num_g_ + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model_->createOutArgs();

  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_f));

  // Sensitivity support (Forward approach only, for now)
  // Jacobian solver required for all sensitivities
  if (modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W)) {
    for (int l = 0; l < num_p_; ++l) {
      // Solution sensitivities: DxDp(l)
      // DfDp(l) required
      const Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_support =
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
      if (!dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
        return outArgs;
      }
      const bool dxdp_linOpSupport =
        dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
      const bool dxdp_mvJacSupport =
        dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
      Thyra::ModelEvaluatorBase::DerivativeSupport dxdp_support;
      if (dxdp_linOpSupport) {
        dxdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
      }
      if (dxdp_mvJacSupport) {
        dxdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
      }
      outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, num_g_, l, dxdp_support);

      // Response sensitivities: DgDp(j, l)
      // DxDp(l) required
      if (dxdp_linOpSupport || dxdp_mvJacSupport) {
        for (int j = 0; j < num_g_; ++j) {
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
            outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l, dgdp_support);
          }
        }
      }
    }
  }

  return outArgs;
}


template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::TransientSolver<Scalar>::create_DgDp_op_impl(int j, int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > num_g_ || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TransientSolver::create_DgDp_op_impl():  " <<
      "Invalid response index j = " <<
      j << "\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_p_ || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TransientSolver::create_DgDp_op_impl():  " <<
      "Invalid parameter index l = " <<
      l << "\n");
  const Teuchos::Array<Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > > dummy =
    Teuchos::tuple(Thyra::zero<Scalar>(this->get_g_space(j), this->get_p_space(l)));
  if (j == num_g_)  {
    return Thyra::defaultMultipliedLinearOp<Scalar>(dummy);
  } else {
    return Teuchos::rcp(new Thyra::DefaultAddedLinearOp<Scalar>(dummy));
  }
}


template <typename Scalar>
const Thyra::ModelEvaluator<Scalar> &
Piro::TransientSolver<Scalar>::getModel() const 
{
  return model_; 
}

template <typename Scalar>
int 
Piro::TransientSolver<Scalar>::num_p() const 
{
  return num_p_; 
}

template <typename Scalar> 
int 
Piro::TransientSolver<Scalar>::num_g() const 
{
  return num_g_; 
}

template <typename Scalar>
void 
Piro::TransientSolver<Scalar>::setSensitivityMethod(const std::string& sensitivity_method_string)
{
  if (sensitivity_method_string == "None") sensitivityMethod_ = NONE; 
  else if (sensitivity_method_string == "Forward") sensitivityMethod_ = FORWARD;
  else if (sensitivity_method_string == "Adjoint") sensitivityMethod_ = ADJOINT;
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
        "\n Error! Piro::TransientSolver: invalid Sensitivity Method = " << sensitivity_method_string << "! \n" 
        << " Valid options for Sensitivity Method are 'None', 'Forward' and 'Adjoint'.\n");
  }
}

template <typename Scalar>
Piro::SENS_METHOD 
Piro::TransientSolver<Scalar>::getSensitivityMethod()
{
  return sensitivityMethod_; 
}

template <typename Scalar>
void 
Piro::TransientSolver<Scalar>::setPiroTempusIntegrator(Teuchos::RCP<const Piro::TempusIntegrator<Scalar>> piroTempusIntegrator)
{
  piroTempusIntegrator_ = piroTempusIntegrator;
}


template <typename Scalar>
void 
Piro::TransientSolver<Scalar>::evalConvergedModelResponsesAndSensitivities(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  *out_ << "\nF) Calculate responses ...\n";

  // Solution at convergence is the response at index num_g_
  RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g_);
  if (Teuchos::nonnull(gx_out)) {
    Thyra::copy(*modelInArgs.get_x(), gx_out.ptr());
  }

  // Setup output for final evalution of underlying model
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model_->createOutArgs();
   
  //Responses 
  for (int j=0; j<num_g_; j++) { 
    auto g_out = outArgs.get_g(j); 
    if (g_out != Teuchos::null) {
      Thyra::put_scalar(Teuchos::ScalarTraits<Scalar>::zero(), g_out.ptr());
      modelOutArgs.set_g(j, g_out);
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
          } 
          else {
            dgdx_request.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
          }
        }
      }
    }

    if (!dgdx_request.none()) {
      Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_deriv;
      if (dgdx_request.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)) {
        dgdx_deriv = Thyra::create_DgDx_mv(*model_, j, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
      } 
      else if (dgdx_request.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
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
        } 
        else {
          model_dgdp_deriv = dgdp_deriv;
        }
        if (!model_dgdp_deriv.isEmpty()) {
          modelOutArgs.set_DgDp(j, l, model_dgdp_deriv);
        }
      }
    }
  }
    
  model_->evalModel(modelInArgs, modelOutArgs);

  // Check if sensitivities are requested 
  bool requestedSensitivities = false;
  for (int i=0; i<num_p_; i++) {
    for (int j=0; j<=num_g_; j++) {
      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none() && !outArgs.get_DgDp(j,i).isEmpty()) {
        requestedSensitivities = true;
        break;
      }
    }
  }

  // Calculate response sensitivities
  if (requestedSensitivities == true) {

    if (sensitivityMethod_ == NONE) {

      //If sensitivities are requested but 'Sensitivity Method' is set to 'None', throw an error.
      TEUCHOS_TEST_FOR_EXCEPTION(
          true,
          Teuchos::Exceptions::InvalidParameter,
          "\n Error! Piro::TransientSolver: you have specified 'Sensitivity Method = None' yet the model supports suggest " 
          << "sensitivities are requested.  Please change 'Sensitivity Method' to 'Forward' or 'Adjoint'\n");
    }
    //
    *out_ << "\nF) Calculate response sensitivities...\n";
 
    switch(sensitivityMethod_) {
      case NONE: //no sensitivities
        break; 

      case FORWARD : //forward sensitivities
      {
        //Get dxdp_mv from Tempus::ForwardIntegratorSensitivity class  
        const RCP<const Thyra::MultiVectorBase<Scalar> > dxdp_mv = piroTempusIntegrator_->getDxDp();
#ifdef DEBUG_OUTPUT
        *out_ << "\n*** Piro::TransientSolver: num_p, num vecs in dxdp = " << num_p_ << ", " << dxdp_mv->domain()->dim() << " ***\n";
#endif
        for (int i=0; i < dxdp_mv->domain()->dim(); ++i) { 
          Teuchos::RCP<const Thyra::VectorBase<Scalar>> dxdp = dxdp_mv->col(i);
#ifdef DEBUG_OUTPUT
          *out_ << "\n*** Piro::TransientSolver dxdp for p = " << i << " ***\n";
          Teuchos::Range1D range;
          RTOpPack::ConstSubVectorView<Scalar> dxdpv;
          dxdp->acquireDetachedView(range, &dxdpv);
          auto dxdpa = dxdpv.values();
          for (auto j = 0; j < dxdpa.size(); ++j) *out_ << dxdpa[j] << " ";
          *out_ << "\n*** Piro::TransientSolver dxdp for p = " << i << " ***\n";
#endif
        }
        //IMPORTANT REMARK: we are currently NOT using DxdotDp and DxdotdotDp in transient sensitivities!  
        //The capability to use them can be added at a later point in time, if desired. 
        //IKT, 5/10/20: throw error if dxdp_mv returned by Tempus is null.  Not sure if this can happen in practice or not...
        if (dxdp_mv == Teuchos::null) {
           TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
              "\n Error! Piro::TransientSolver: DxDp returned by Tempus::IntegratorForwardSensitivity::getDxDp() routine is null!\n"); 
        } 
        //IKT FIXME: probably a lot of this can be reused for adjoint sensitivities.  Will clean up later, 
        //when adjoint sensitivities are implemented. 
        for (int l = 0; l < num_p_; ++l) {
          for (int j = 0; j < num_g_; ++j) {
            //Get DgDp and DgDx 
            const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
               outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
            if (!dgdp_support.none()) {
              const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv = outArgs.get_DgDp(j, l);
              if (!dgdp_deriv.isEmpty()) {
                const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_deriv = modelOutArgs.get_DgDx(j);
                const RCP<const Thyra::MultiVectorBase<Scalar> > dgdx_mv = dgdx_deriv.getMultiVector();
                RCP<const Thyra::LinearOpBase<Scalar> > dgdx_op = dgdx_deriv.getLinearOp();
                if (Teuchos::is_null(dgdx_op)) {
                  //NOTE: dgdx_mv is the transpose, so by calling Thyra::adjoint on dgdx_mv, 
                  //we get the untransposed operator back as dgdx_op
                  dgdx_op = Thyra::adjoint<Scalar>(dgdx_mv);
                }
                const RCP<Thyra::LinearOpBase<Scalar> > dgdp_op = dgdp_deriv.getLinearOp();
                if (Teuchos::nonnull(dgdp_op)) {
                  //Case 1: DgDp, DgDx and DxDp are linear ops.  This corresponds to a non-scalar
                  //response and distributed parameters.  Tempus::ForwardSensitivityIntegrator 
                  //cannot return a LinearOp for DxDp.  Therefore this case is not relevant here.
                  TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                    "\n Error! Piro::TransientSolver: DgDp = DERIV_LINEAR_OP is not supported with forward sensitivities!");
                }
                //Cases 2 and 3 below correspond to scalar responses and scalar parameters.  These can happen
                //with dgdp = DERIV_MV_GRADIENT_FORM and dgpg = DERIV_MV_JACOBIAN_FORM.  For 
                //DERIV_MV_GRADIENT_FORM, the map is the responses, and the columns are the parameters; for 
                //DERIV_MV_JACOBIAN_FORM, the map is the parameters, and the columns are the responses.
                const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv = dgdp_deriv.getMultiVector();
                if (Teuchos::nonnull(dgdp_mv)) {
                  if (dgdp_deriv.getMultiVectorOrientation() == Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) {
                    //Case 2: DgDp = DERIV_MV_GRADIENT_FORM, DgDx is MV, DxDp is MV.
                    //This corresponds to a scalar response and distributed parameters 
                    //[dgdp_mv]^T = [dx/dp_mv]^T*dg/dx_mv + [dg/dp_mv]^T
                    //Note: Gradient form stores transpose of derivative in dgdx_mv (DERIV_MV_GRADIENT_FORM is transposed!) 
                    Thyra::apply(*dxdp_mv, Thyra::TRANS, *dgdx_mv, dgdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::one(),
                                 Teuchos::ScalarTraits<Scalar>::one());
                  } 
                  else {
                    //Case 3: DgDp = DERIV_MV_JACOBIAN_FORM (the alternate to DERIV_MV_GRADIENT_FORM for getMultiVectorOrientation),
                    //DgDx = DERIV_LINEAR_OP and DxDp is MV.  Note that DgDx implementes a DERIV_LINEAR_OP for MVs, 
                    //so there is no contradiction here in the type.   
                    //dgdp_mv = dg/dx_op*dx/dp_mv + dg/dp_mv
                    Thyra::apply(*dgdx_op, Thyra::NOTRANS, *dxdp_mv, dgdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::one(),
                                 Teuchos::ScalarTraits<Scalar>::one());
                  }
                }
              }
            }
          }
        }
        break; 
      }
    case ADJOINT: //adjoint sensitivities - not yet implemented 
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
        "\n Error! Piro::TransientSolver: adjoint sentivities (Sensitivity Method = "
        << "Adjoint) are not yet supported!  Please set 'Sensitivity Method' to 'None' or 'Forward'.\n");
      break; 
    }
  }
}
