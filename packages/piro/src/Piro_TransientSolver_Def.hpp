// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model, 
  const Teuchos::RCP<Teuchos::ParameterList> &appParams,
  const Teuchos::RCP<Piro::ObserverBase<Scalar> > &piroObserver) :
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  model_(model),
  piroObserver_(piroObserver),
  num_p_(model->Np()), 
  num_g_(model->Ng()),
  sensitivityMethod_(NONE)
{
  Teuchos::RCP<Teuchos::ParameterList> tempusPL = sublist(appParams, "Tempus", true);
  if (tempusPL->isSublist("Sensitivities")){
    Teuchos::ParameterList& tempusSensPL = tempusPL->sublist("Sensitivities", true);
    if (sensitivityMethod_ != ADJOINT) {
      response_fn_index_ = tempusSensPL.get<int>("Response Function Index", 0);
    }
    if (sensitivityMethod_ != NONE) {
      sens_param_index_ = tempusSensPL.get<int>("Sensitivity Parameter Index", 0);
    }
  }
}

template <typename Scalar>
Piro::TransientSolver<Scalar>::TransientSolver(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
  const int sens_param_index, 
  const int response_fn_index) :
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  model_(model), 
  num_p_(model->Np()), 
  num_g_(model->Ng()),
  sensitivityMethod_(NONE)
{
  if (num_g_ == 1) {
    response_fn_index_ = 0; 
  }
  else {
    response_fn_index_ = response_fn_index; 
    if ((sensitivityMethod_ == ADJOINT) && (response_fn_index_ == -1)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true,
          Teuchos::Exceptions::InvalidParameter,
          "\n Error in Piro::TransientSolver constructor: 'Response Function Index' must be specified for ADJOINT sensitivity method "
          << "with >1 response!\n"); 
    }
  }
  if (num_p_ == 1) {
    sens_param_index_ = 0; 
  }
  else {
    sens_param_index_ = sens_param_index; 
    if ((sensitivityMethod_ != NONE) && (sens_param_index_ == -1)) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true,
          Teuchos::Exceptions::InvalidParameter,
          "\n Error in Piro::TransientSolver constructor: 'Parameter Sensitivity Index' must be specified for ADJOINT and FORWARD sensitivity method "
          << "with >1 parameter!\n"); 
    }
  }
}

template<typename Scalar>
void Piro::TransientSolver<Scalar>::resetSensitivityParamIndex(const int sens_param_index)
{
  sens_param_index_ = sens_param_index; 
}

template<typename Scalar>
void Piro::TransientSolver<Scalar>::resetResponseFnIndex(const int response_fn_index)
{
  response_fn_index_ = response_fn_index;  
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
      //Special case that comes up when not doing forward sensitivity calculations
      if (sensitivityMethod_ == FORWARD) {
        if (!dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) { 
          return outArgs;
        }
      }
      // DxDp has same support as DfDp 
      const bool dxdp_linOpSupport =
          dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
      const bool dxdp_mvJacobSupport =
          dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
      const bool dxdp_mvGradSupport =
          dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
      {
        Thyra::ModelEvaluatorBase::DerivativeSupport dxdp_support;
        if (dxdp_linOpSupport) {
          dxdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        }
        if (dxdp_mvJacobSupport) {
          dxdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        }
        if (dxdp_mvGradSupport) {
          dxdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        }
        outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, num_g_, l, dxdp_support);
      }

      // Response sensitivities: DgDp(j, l)
      // DxDp(l) required
      if (dxdp_linOpSupport || dxdp_mvJacobSupport || dxdp_mvGradSupport) {
        for (int j = 0; j < num_g_; ++j) {
          // DgDx(j) required
          const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
              modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j);
	  // IKT, 6/30/2021: note that DERIV_LINEAR_OP is not relevant for DgDx for the 
	  // current use cases in Albany but keeping it nonetheless for completeness and
	  // consistency with Piro::SteadySolver
          const bool dgdx_linOpSupport =
              dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
          const bool dgdx_mvGradSupport =
              dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
          const bool dgdx_mvJacobSupport =
              dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
          if (dgdx_linOpSupport || dgdx_mvGradSupport || dgdx_mvJacobSupport) {
            // Dgdp(j, l) required
            const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
                modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
            Thyra::ModelEvaluatorBase::DerivativeSupport total_dgdp_support;
	    //IKT. 6/30/2021: DERIV_LINEAR_OP is not relevant (and actually not supported)
	    //for Tempus sensitivities but keeping it nonethless for completeness and 
	    //consistency with Piro::SteadySolver. There is a throw for this case in 
	    //evalConvergedModelResponsesAndSensitivities(). 
            if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP) &&
                dgdx_linOpSupport && dxdp_linOpSupport) {
              total_dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
            }
            //if (dxdp_mvJacobSupport && dgdx_mvJacobSupport) {
            if (sensitivityMethod_ == FORWARD) {
              if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
                total_dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
              }
            }
            //if (dxdp_mvGradSupport && dgdx_mvGradSupport) {
            if (sensitivityMethod_ == ADJOINT) {
              if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)) {
                total_dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
              }
            }
            outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l, total_dgdp_support);
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
  return *model_;
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
    *out_ << "\nG) Calculate response sensitivities...\n";
 
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
        const int l = sens_param_index_; 
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
                  "\n Error! Piro::TransientSolver: DgDp = DERIV_LINEAR_OP (relevant for distributed responses) is not supported with forward sensitivities!");
              }
              //Cases 2 and 3.  These can happen with dgdp = DERIV_MV_GRADIENT_FORM and dgpg = DERIV_MV_JACOBIAN_FORM.  For 
              //DERIV_MV_GRADIENT_FORM, the map is the responses, and the columns are the parameters; for 
              //DERIV_MV_JACOBIAN_FORM, the map is the parameters, and the columns are the responses.
	      //Note that Case 2, which assumes distributed parameters, would not get called for forward sensitivities
	      //as it would be slow, but it could be called in theory. 
              const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv = dgdp_deriv.getMultiVector();
	      //IKT, question: is it worth throwing if dgdp_mv == null, or this cannot happen?
              if (Teuchos::nonnull(dgdp_mv)) {
                if (dgdp_deriv.getMultiVectorOrientation() == Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) { //case 2
                  //Case 2: DgDp = DERIV_MV_GRADIENT_FORM, DgDx is MV, DxDp is MV.
                  //This corresponds to a scalar response and distributed parameters
                  //[dgdp_mv]^T = [dx/dp_mv]^T*dg/dx_mv + [dg/dp_mv]^T
                  //Note: Gradient form stores transpose of derivative in dgdx_mv (DERIV_MV_GRADIENT_FORM is transposed!) 
	          //Note that forward sensitivities will not be called for distributed parameters, therefore this case
	          //is not really relevant here. 
                  Thyra::apply(*dxdp_mv, Thyra::TRANS, *dgdx_mv, dgdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::one(),
                               Teuchos::ScalarTraits<Scalar>::one());
                } 
                else { //case 3
                  //Case 3: DgDp = DERIV_MV_JACOBIAN_FORM (the alternate to DERIV_MV_GRADIENT_FORM for getMultiVectorOrientation),
                  //DgDx = DERIV_LINEAR_OP (for distributed responses) or DERIV_MV_JACOBIAN_FORM (for scalar responses), 
	          //and DxDp is MV.  Note that DgDx implementes a DERIV_LINEAR_OP for MVs, so there is no contradiction here in the type,
		  //and this case encompasses both distributed and scalar responses.   
                  //dgdp_mv = dg/dx_op*dx/dp_mv + dg/dp_mv
                  Thyra::apply(*dgdx_op, Thyra::NOTRANS, *dxdp_mv, dgdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::one(),
                               Teuchos::ScalarTraits<Scalar>::one());
                }
              }
            }
          }
        }
        break; 
      }
    case ADJOINT: //adjoint sensitivities
      const int l = sens_param_index_; 
      //Get DgDp from outArgs and set it based on adjoint integrator from Tempus 
      //Note that one could return DgDp for a single parameter and response by setting
      //const int j = response_fn_index_, but this is not done now.
      for (int j = 0; j < num_g_; ++j) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
           outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        if (!dgdp_support.none()) {
          const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv = outArgs.get_DgDp(j, l);
          if (!dgdp_deriv.isEmpty()) {
            const RCP<Thyra::LinearOpBase<Scalar> > dgdp_op = dgdp_deriv.getLinearOp();
            if (Teuchos::nonnull(dgdp_op)) {
              //Case 1: DgDp is a linear ops.  This corresponds to a non-scalar
              //response and distributed parameters.  Tempus::AdjointSensitivityIntegrator 
              //cannot return a LinearOp for DgDp.  Therefore this case is not relevant here.
              TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
              "\n Error! Piro::TransientSolver: DgDp = DERIV_LINEAR_OP (relevant for distributed responses) is not supported with adjoint sensitivities!");
            }
            //Cases 2 and 3.  These can happen with dgdp = DERIV_MV_GRADIENT_FORM and dgpg = DERIV_MV_JACOBIAN_FORM.  For 
            //DERIV_MV_GRADIENT_FORM, the map is the responses, and the columns are the parameters; for 
            //DERIV_MV_JACOBIAN_FORM, the map is the parameters, and the columns are the responses.
            //Both cases are relevant for adjoint sensitivities: Case 2 corresponds to distributed parameters, whereas
            //case 3 correspondes to scalar parameters.
            const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv = dgdp_deriv.getMultiVector();
	    //IKT, question: is it worth throwing if dgdp_mv == null, or this cannot happen?
            if (Teuchos::nonnull(dgdp_mv)) {
	      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> dgdp_mv_from_tempus = piroTempusIntegrator_->getDgDp();
              if (dgdp_mv_from_tempus == Teuchos::null) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                  "\n Error! Piro::TransientSolver: DgDp returned by Tempus::IntegratorAdjointSensitivity::getDgDp() routine is null!\n"); 
              } 
	      //Copy dgdp_mv_from_tempus into dgdp_mv - IKT, there may be better way to do this
	      dgdp_mv->assign(*dgdp_mv_from_tempus);
	      //Uncomment to observe DgDp from within Piro
	      /*if (piroObserver_ != Teuchos::null) {
	        std::cout << "IKT start observing dgdp\n";
	        //Observe also the solution, since observeSolution requires passing this field
	        //This would be relevant if observing DgDp to a separate file, in which it may be useful to 
	        //also have the solution.
                Teuchos::RCP<const Thyra::VectorBase<Scalar>> solution = piroTempusIntegrator_->getX(); 
                piroObserver_->observeSolution(*solution, *dgdp_mv_from_tempus, piroTempusIntegrator_->getTime());
	        std::cout << "IKT end observing dgdp\n";
	      }*/
	    }
	  }
        }
      }
      break; 
    }
  }
}
