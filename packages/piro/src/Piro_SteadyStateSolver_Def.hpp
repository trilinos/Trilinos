// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_STEADYSTATESOLVER_DEF_HPP
#define PIRO_STEADYSTATESOLVER_DEF_HPP

#include "Piro_SteadyStateSolver.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Teuchos_TimeMonitor.hpp"

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

#include <Piro_config.hpp>

#ifdef HAVE_PIRO_ROL
#include "ROL_ThyraVector.hpp"
#include "ROL_ScaledThyraVector.hpp"
#include "ROL_Thyra_BoundConstraint.hpp"
#include "ROL_ThyraME_Objective.hpp"
#include "ROL_ThyraProductME_Objective.hpp"
#include "Piro_ThyraProductME_Objective_SimOpt.hpp"
#include "Piro_ThyraProductME_Constraint_SimOpt.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#endif

template <typename Scalar>
Piro::SteadyStateSolver<Scalar>::
SteadyStateSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &adjointModel) :
  model_(model),
  adjointModel_(adjointModel),
  num_p_(model->Np()),
  num_g_(model->Ng()),
  sensitivityMethod_(NONE)
{}

template <typename Scalar>
Piro::SteadyStateSolver<Scalar>::
SteadyStateSolver(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &adjointModel,
    int numParameters) :
  model_(model),
  adjointModel_(adjointModel),
  num_p_(numParameters),
  num_g_(model->Ng()),
  sensitivityMethod_(NONE)
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

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::SteadyStateSolver<Scalar>::get_x_space() const
{
  return model_->get_x_space();
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::SteadyStateSolver<Scalar>::get_f_space() const
{
  return model_->get_f_space();
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::SteadyStateSolver<Scalar>::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> result;
  result.setModelEvalDescription(this->description());
  result.set_Np_Ng(num_p_, num_g_);
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = model_->createInArgs();
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, modelInArgs.supports(Thyra::ModelEvaluatorBase::IN_ARG_x));
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

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_f));

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

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
  for (int i=0; i<num_g_; i++) {
    for (int j1=0; j1<num_p_; j1++) {
      for (int j2=0; j2<num_p_; j2++) {
        result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp,
        i,
        j1,
        j2,
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, i, j1, j2));
      }
    }
  }
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

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
void
Piro::SteadyStateSolver<Scalar>::setSensitivityMethod(const std::string& sensitivity_method_string)
{
  if (sensitivity_method_string == "None") sensitivityMethod_ = NONE;
  else if (sensitivity_method_string == "Forward") sensitivityMethod_ = FORWARD;
  else if (sensitivity_method_string == "Adjoint") sensitivityMethod_ = ADJOINT;
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
        "\n Error! Piro::SteadyStateSolver: invalid Sensitivity Method = " << sensitivity_method_string << "! \n"
        << " Valid options for Sensitivity Method are 'None', 'Forward' and 'Adjoint'.\n");
  }
}

template <typename Scalar>
Piro::SENS_METHOD
Piro::SteadyStateSolver<Scalar>::getSensitivityMethod()
{
  return sensitivityMethod_;
}

template <typename Scalar>
void Piro::SteadyStateSolver<Scalar>::evalConvergedModelResponsesAndSensitivities(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs,
    Teuchos::ParameterList& appParams) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Solution at convergence is the response at index num_g_
  {
    const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g_);
    if (Teuchos::nonnull(gx_out)) {
      Thyra::copy(*modelInArgs.get_x(), gx_out.ptr());
    }
  }

  // Setup output for final evalution of underlying model
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model_->createOutArgs();

  // Responses
  for (int j = 0; j < num_g_; ++j) {
    const RCP<Thyra::VectorBase<Scalar> > g_out = outArgs.get_g(j);
    if (Teuchos::nonnull(g_out)) {
      Thyra::put_scalar(0.0, g_out.ptr());
      modelOutArgs.set_g(j, g_out);
    }
  }

  bool compute_sensitivities = false;
  for (int j = 0; j <= num_g_; ++j) { // resize
    for (int l = 0; l < num_p_; ++l) {
      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l).none() && !outArgs.get_DgDp(j, l).isEmpty()) {
        compute_sensitivities = true;
      }
    }
  }

  bool computeForwardSensitivities = compute_sensitivities && (sensitivityMethod_==FORWARD);
  bool computeAdjointSensitivities = compute_sensitivities && (sensitivityMethod_==ADJOINT);

#ifdef HAVE_PIRO_ROL
  if(computeAdjointSensitivities) {
    double tol = 1e-8;

    RCP<Thyra::VectorSpaceBase<Scalar> const> p_space = this->getModel().get_p_space(0);
    RCP<Thyra::VectorBase<Scalar>> thyra_p = Thyra::createMember(p_space);

    RCP<const Thyra::VectorBase<Scalar> > p_init = modelInArgs.get_p(0) != Teuchos::null ? modelInArgs.get_p(0) : this->getModel().getNominalValues().get_p(0);
    Thyra::copy(*p_init, thyra_p.ptr());

    ROL::ThyraVector<Scalar> rol_p(thyra_p);

    Teuchos::RCP<Thyra::VectorSpaceBase<Scalar> const> x_space = this->getModel().get_x_space();
    Teuchos::RCP<Thyra::VectorBase<Scalar>> x = Thyra::createMember(x_space);
    Thyra::copy(*modelInArgs.get_x(), x.ptr());

    ROL::ThyraVector<Scalar> rol_x(x);
    Teuchos::RCP<Thyra::VectorBase<Scalar>> lambda_vec = Thyra::createMember(x_space);
    ROL::ThyraVector<Scalar> rol_lambda(lambda_vec);

    RCP<Thyra::VectorBase<Scalar> > current_g = Thyra::createMember(p_space);
    ROL::ThyraVector<Scalar> rol_current_g(current_g);

    ROL::Ptr<ROL::Vector<Scalar> > rol_p_ptr = ROL::makePtrFromRef(rol_p);
    ROL::Ptr<ROL::Vector<Scalar> > rol_x_ptr = ROL::makePtrFromRef(rol_x);
    ROL::Ptr<ROL::Vector<Scalar> > rol_lambda_ptr = ROL::makePtrFromRef(rol_lambda);


    Piro::ThyraProductME_Constraint_SimOpt<Scalar> constr(model_, adjointModel_, appParams, Teuchos::VERB_NONE);
    auto  stateStore = ROL::makePtr<ROL::VectorController<Scalar>>();
      
    for (int i=0; i<num_g_; ++i) {      

      Piro::ThyraProductME_Objective_SimOpt<Scalar> obj(model_, i, appParams, Teuchos::VERB_NONE);

      ROL::Ptr<ROL::Objective_SimOpt<Scalar> > obj_ptr = ROL::makePtrFromRef(obj);
      ROL::Ptr<ROL::Constraint_SimOpt<Scalar> > constr_ptr = ROL::makePtrFromRef(constr);

      //create the ROL reduce objective initializing it with the current state and parameter
      ROL::Reduced_Objective_SimOpt<Scalar> reduced_obj(obj_ptr,constr_ptr,stateStore,rol_x_ptr,rol_p_ptr,rol_lambda_ptr);
      reduced_obj.update(rol_p,ROL::UpdateType::Temp);
      stateStore->set(*rol_x_ptr, std::vector<Scalar>());  //second argument not meaningful for deterministic problems 

      //reduced_obj.set_precomputed_state(rol_x,rol_p);
      if(i>0) { //a bit hacky, but the jacobian and, if needed, its adjoint have been computed at iteration 0
        constr.computeJacobian1_ = false;
        constr.computeAdjointJacobian1_ = false;
      }

      Scalar tmp = reduced_obj.value(rol_p,tol);
      reduced_obj.gradient(rol_current_g, rol_p, tol);

      RCP<Thyra::VectorBase<Scalar> > g_out = outArgs.get_g(i);
      Thyra::set_ele(0,tmp,g_out.ptr());


      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, 0).none() &&
          !outArgs.get_DgDp(i,0).isEmpty()) {

        RCP<Thyra::MultiVectorBase<Scalar> > dgdp_out = outArgs.get_DgDp(i,0).getMultiVector();
        Thyra::assign(dgdp_out->col(0).ptr(), *current_g);
      }
    }
    return;
  }
#endif

  if(computeForwardSensitivities)
  {
    // Jacobian
    if (compute_sensitivities) {
      const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
          model_->create_W();
      modelOutArgs.set_W(jacobian);
    }

    // DfDp derivatives
    for (int l = 0; l < num_p_; ++l) {
      Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_request;
      for (int j = 0; j <= num_g_; ++j) {
        if(!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l).none()) {
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
        } else if (dfdp_request.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)) {
          dfdp_deriv = Thyra::create_DfDp_mv(*model_, l, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        } else if (dfdp_request.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          dfdp_deriv = model_->create_DfDp_op(l);
          RCP<Thyra::LinearOpBase<Scalar> > dfdp_op =
              dfdp_deriv.getLinearOp();
          if (Teuchos::is_null(dfdp_op)) {
            TEUCHOS_TEST_FOR_EXCEPTION(
                true, std::logic_error,
                std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
                "For df/dp, getLinearOp is null. " << std::endl);
          }
        }
        modelOutArgs.set_DfDp(l, dfdp_deriv);
      }
    }

    // DgDx derivatives
    for (int j = 0; j < num_g_; ++j) {
      Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_request;
      for (int l = 0; l < num_p_; ++l) {
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l).none()) {
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
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l).none()) {
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
  } else if(computeAdjointSensitivities) {

    // Compute adjoint layouts of df/dp, dg/dx depending on
    for (int i=0; i<num_p_; i++) {
      // p
      //modelInArgs.set_p(i, inArgs.get_p(i));

      // df/dp
      bool compute_sensitivities_wrt_pi = false;
      for (int j=0; j<=num_g_; j++) {
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none() &&
            !outArgs.get_DgDp(j,i).isEmpty()) {
          compute_sensitivities_wrt_pi = true;
        }
      }

      if (compute_sensitivities_wrt_pi) {
        auto p_space = this->getModel().get_p_space(i);
        auto f_space = this->getModel().get_f_space();
        auto p_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(p_space);
        auto f_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(f_space);
        Thyra::ModelEvaluatorBase::DerivativeSupport ds =  modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,i);

        // Determine which layout to use for df/dp.  Ideally one would look
        // at num_params, num_resids, what is supported by the underlying
        // model evaluator, and the sensitivity method, and make the best
        // choice to minimze the number of solves.  However this choice depends
        // also on what layout of dg/dx is supported (e.g., if only the operator
        // form is supported for forward sensitivities, then df/dp must be
        // DERIV_MV_JACOBIAN_FORM).  For simplicity, we order the conditional tests
        // to get the right layout in most situations.

        if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          Thyra::ModelEvaluatorBase::Derivative<Scalar> dfdp_deriv = this->getModel().create_DfDp_op(i);
          RCP<Thyra::LinearOpBase<Scalar> > dfdp_op = dfdp_deriv.getLinearOp();
          TEUCHOS_TEST_FOR_EXCEPTION(
              dfdp_op == Teuchos::null, std::logic_error,
              std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
              "Needed df/dp operator (" << i << ") is null!" << std::endl);
          modelOutArgs.set_DfDp(i,dfdp_deriv);
        } else {
          /*
          TEUCHOS_TEST_FOR_EXCEPTION(
            !ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP),
            std::logic_error,
            std::endl <<
            "ROL::Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
            "The code related to df/dp multivector has been commented out because never tested.  " <<
            std::endl);

          if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) &&
              f_space_plus->isLocallyReplicated()) {
            auto dfdp = Thyra::createMembers(p_space, f_space->dim());

            Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
            dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
            modelOutArgs.set_DfDp(i,dmv_dfdp);
          } else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) &&
                     p_space_plus->isLocallyReplicated()) {
            auto dfdp = Thyra::createMembers(f_space, p_space->dim());
            Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
            dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
            modelOutArgs.set_DfDp(i,dmv_dfdp);
          } else
            TEUCHOS_TEST_FOR_EXCEPTION(
                true, std::logic_error,
                std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
                "For df/dp(" << i <<") with adjoint sensitivities, " <<
                "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
                "DERIV_MV_JACOBIAN_FORM with p not distributed, or "
                "DERIV_MV_GRADIENT_FORM with f not distributed." <<
                std::endl);
          */
        }
      }
    }

    for (int j=0; j<num_g_; j++) {
      // dg/dx
      bool compute_gj_sensitivities = false;
      for (int i=0; i<num_p_; i++) {
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none() &&
            !outArgs.get_DgDp(j,i).isEmpty()) {
          compute_gj_sensitivities = true;
        }
      }
      if (compute_gj_sensitivities) {
        auto g_space = this->getModel().get_g_space(j);
        auto x_space = this->getModel().get_x_space();
        int num_responses = g_space->dim();
        int num_solution = x_space->dim();
        auto g_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(g_space);
        auto x_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(x_space);
        bool g_dist = !g_space_plus->isLocallyReplicated();//g_space->DistributedGlobal();
        bool x_dist = !x_space_plus->isLocallyReplicated();//x_space->DistributedGlobal();


        Thyra::ModelEvaluatorBase::DerivativeSupport ds =  modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx,j);
        if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && !g_dist) {
          auto dgdx = Thyra::createMembers(x_space, num_responses);
          Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
          dmv_dgdx(dgdx, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
          modelOutArgs.set_DgDx(j,dmv_dgdx);
        } else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) && !x_dist) {
          auto dgdx = Thyra::createMembers(g_space, num_solution);
          Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
          dmv_dgdx(dgdx, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
          modelOutArgs.set_DgDx(j,dmv_dgdx);
        } else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          auto dgdx_op = this->getModel().create_DgDx_op(j);
          TEUCHOS_TEST_FOR_EXCEPTION(
              dgdx_op == Teuchos::null, std::logic_error,
              std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
              "Needed dg/dx operator (" << j << ") is null!" << std::endl);
          modelOutArgs.set_DgDx(j,dgdx_op);
        } else
          TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
              "For dg/dx(" << j <<") with adjoint sensitivities, " <<
              "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
              "DERIV_MV_JACOBIAN_FORM with x not distributed, or "
              "DERIV_MV_GRADIENT_FORM with g not distributed." <<
              std::endl);

        // dg/dp
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,j,0).none()) {
          Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp = outArgs.get_DgDp(j,0);
          if (dgdp.getLinearOp() != Teuchos::null) {
            auto p_space = this->getModel().get_p_space(0);
            int num_params = p_space->dim();
            auto p_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(p_space);
            bool p_dist = !p_space_plus->isLocallyReplicated();//p_space->DistributedGlobal();
            Thyra::ModelEvaluatorBase::DerivativeSupport ds_dgdp = modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,j,0);
            if (ds_dgdp.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
              auto dgdp_op =
                  this->getModel().create_DgDp_op(j,0);
              TEUCHOS_TEST_FOR_EXCEPTION(
                  dgdp_op == Teuchos::null, std::logic_error,
                  std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
                  "Needed dg/dp operator (" << j << "," << 0 << ") is null!" <<
                  std::endl);
              modelOutArgs.set_DgDp(j,0,dgdp_op);
            }
            else if (ds_dgdp.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) && !p_dist) {
              auto tmp_dgdp = createMembers(g_space, num_params);
              Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
              dmv_dgdp(tmp_dgdp, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
              modelOutArgs.set_DgDp(j,0,dmv_dgdp);
            }
            else if (ds_dgdp.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && !g_dist) {
              auto tmp_dgdp = createMembers(p_space, num_responses);
              Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
              dmv_dgdp(tmp_dgdp, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
              modelOutArgs.set_DgDp(j,0,dmv_dgdp);
            }
            else
              TEUCHOS_TEST_FOR_EXCEPTION(
                  true, std::logic_error,
                  std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
                  "For dg/dp(" << j << "," << 0 <<
                  ") with operator sensitivities, "<<
                  "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
                  "DERIV_MV_JACOBIAN_FORM with p not distributed, or "
                  "DERIV_MV_GRADIENT_FORM with g not distributed." <<
                  std::endl);
          }
          else
            modelOutArgs.set_DgDp(j,0,outArgs.get_DgDp(j,0));
        }
      }
    }
  }

  // Calculate g, df/dp, dg/dp, dg/dx by evaluating the underlying model
  {
    const auto timer = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities::calc g, df/dp, dg/dp, dg/dx")));
    model_->evalModel(modelInArgs, modelOutArgs);
  }

  if(computeForwardSensitivities) {

    // Assemble user-requested sensitivities
    if (modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W)) {
      const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
          modelOutArgs.get_W();
      if (Teuchos::nonnull(jacobian)) {
        for (int l = 0; l < num_p_; ++l) {
          const Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_support =
              modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
          if (!dfdp_support.none()) {
            const bool dfdp_linOpSupport =
                dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
            const bool dfdp_mvJacSupport =
                dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
            const bool dfdp_mvGradSupport =
                dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
            const Thyra::ModelEvaluatorBase::Derivative<Scalar> dfdp_deriv =
                modelOutArgs.get_DfDp(l);
            const RCP<Thyra::MultiVectorBase<Scalar> > dfdp_mv =
                dfdp_deriv.getMultiVector();
            RCP<Thyra::LinearOpBase<Scalar> > dfdp_op =
                dfdp_deriv.getLinearOp();
            if (Teuchos::is_null(dfdp_op)) {
              dfdp_op = dfdp_mv;
            }

            if (Teuchos::is_null(dfdp_mv) && Teuchos::is_null(dfdp_op))
              TEUCHOS_TEST_FOR_EXCEPTION(
                  true, std::logic_error,
                  std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
                  "For df/dp, both getMultiVector and getLinearOp are null. " << dfdp_linOpSupport << " " << dfdp_mvJacSupport << " " << dfdp_mvGradSupport <<
                  std::endl);

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
              if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l).none()) {
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
                      } else if (Teuchos::nonnull(minus_dxdp_mv)) {
                        Thyra::apply(
                            *minus_dxdp_mv,
                            Thyra::TRANS,
                            *dgdx_mv,
                            dgdp_mv.ptr(),
                            -Teuchos::ScalarTraits<Scalar>::one(),
                            Teuchos::ScalarTraits<Scalar>::one());
                      } else {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                            true,
                            std::invalid_argument,
                            "Both dxdp_mv and minus_dxdp_mv are null.\n");
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
                      } else if (Teuchos::nonnull(minus_dxdp_mv)) {
                        Thyra::apply(
                            *dgdx_op,
                            Thyra::NOTRANS,
                            *minus_dxdp_mv,
                            dgdp_mv.ptr(),
                            -Teuchos::ScalarTraits<Scalar>::one(),
                            Teuchos::ScalarTraits<Scalar>::one());
                      } else {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                            true,
                            std::invalid_argument,
                            "Both dxdp_mv and minus_dxdp_mv are null.\n");
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
  } else if (computeAdjointSensitivities) {

    bool availableAdjointModel = Teuchos::nonnull(adjointModel_);

    // Create implicitly transpose Jacobian and preconditioner
    Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory = availableAdjointModel ?  adjointModel_->get_W_factory() : model_->get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
    Teuchos::RCP< Thyra::LinearOpBase<double> > lop = availableAdjointModel ?  adjointModel_->create_W_op() : model_->create_W_op();
    Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<double> > losb = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop));
    Teuchos::RCP< ::Thyra::PreconditionerBase<double> > prec;

    auto in_args = availableAdjointModel ?  adjointModel_->createInArgs() : model_->createInArgs();
    auto out_args = availableAdjointModel ?  adjointModel_->createOutArgs() : model_->createOutArgs();

    Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> > prec_factory =  lows_factory->getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory)) {
      prec = prec_factory->createPrec();
    } else if (out_args.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
      prec = availableAdjointModel ?  adjointModel_->create_W_prec() : model_->create_W_prec();
    }

    const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
        lows_factory->createOp();

    {
      const auto timer = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities::evalModel")));
      in_args.set_x(modelInArgs.get_x());
      out_args.set_W_op(lop);
      if(availableAdjointModel)
        adjointModel_->evalModel(in_args, out_args);
      else 
        model_->evalModel(in_args, out_args);
      in_args.set_x(Teuchos::null);
      out_args.set_W_op(Teuchos::null);
    }

    if (Teuchos::nonnull(prec_factory))
      prec_factory->initializePrec(losb, prec.get());
    else if ( Teuchos::nonnull(prec) && (out_args.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) ) {
      in_args.set_x(modelInArgs.get_x());
      out_args.set_W_prec(prec);
      if(availableAdjointModel)
        adjointModel_->evalModel(in_args, out_args);
      else
        model_->evalModel(in_args, out_args);
    }

    if(Teuchos::nonnull(prec)) {
      if(availableAdjointModel) {
        Thyra::initializePreconditionedOp<double>(*lows_factory,
          lop,
          prec,
          jacobian.ptr());
      } else {
        Thyra::initializePreconditionedOp<double>(*lows_factory,
            Thyra::transpose<double>(lop),
            Thyra::unspecifiedPrec<double>(::Thyra::transpose<double>(prec->getUnspecifiedPrecOp())),
            jacobian.ptr());
      }
    } else {
      if(availableAdjointModel)
        Thyra::initializeOp<double>(*lows_factory, lop, jacobian.ptr());
      else
        Thyra::initializeOp<double>(*lows_factory, Thyra::transpose<double>(lop), jacobian.ptr());
    }


    for (int j=0; j<num_g_; j++) {

      // See if there are any forward sensitivities we need to do
      // that aren't handled by the operator
      bool compute_gj_sensitivities = false;
      for (int i=0; i<num_p_; i++) {
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none()) {
          if (outArgs.get_DgDp(j,i).getMultiVector() != Teuchos::null)
            compute_gj_sensitivities = true;
        }
      }
      if (!compute_gj_sensitivities)
        continue;

      if (!modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j).none()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
            modelOutArgs.get_DgDx(j).getLinearOp()!=Teuchos::null,
            std::logic_error,
            std::endl << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
            "Can\'t use dg/dx operator " << j << " with non-operator " <<
            "adjoint sensitivities." << std::endl);
        auto dgdx  = modelOutArgs.get_DgDx(j).getMultiVector();
        if (dgdx != Teuchos::null) {
          int num_cols = dgdx->domain()->dim();

          // (2) Calculate xbar multivector from -(J^{-T}*dg/dx)

          auto xbar = createMembers(dgdx->range(), num_cols);
          {
            Teuchos::TimeMonitor timer(*Teuchos::TimeMonitor::getNewTimer("Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities::calc xbar"));
            xbar->assign(0);
            Thyra::solve(
                *jacobian,
                Thyra::NOTRANS,
                *dgdx,
                xbar.ptr());
                //,ptr(&solve_criteria));

            Thyra::scale(-1.0, xbar.ptr());
          }
          // (3) Calculate dg/dp^T = df/dp^T*xbar + dg/dp^T
          for (int i=0; i<num_p_; i++) {
            std::string paramName = "Parameter " + std::to_string(i);
            std::ostringstream ss; ss << "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities::calc dg/dp^T(" << paramName << ")";
            const auto timer = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(ss.str())));
            if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none()) {
              auto dgdp_out = outArgs.get_DgDp(j,i).getMultiVector();
              if (dgdp_out != Teuchos::null) {
                Thyra::ModelEvaluatorBase::Derivative<Scalar> dfdp_dv = modelOutArgs.get_DfDp(i);
                auto dfdp_op = dfdp_dv.getLinearOp();
                auto dfdp = dfdp_dv.getMultiVector();
                Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_dv = modelOutArgs.get_DgDx(j);
                if (dfdp_op != Teuchos::null) {
                  Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient =
                      outArgs.get_DgDp(j,i).getMultiVectorOrientation();
                  bool transpose = false;
                  if (dgdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)
                    transpose = true;
                  auto tmp = createMembers(dfdp_op->domain(),
                      xbar->domain()->dim());

                  dfdp_op->apply(Thyra::TRANS,*xbar, tmp.ptr(),1.0, 0.0);
                  auto dgdp_range = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(dgdp_out->range());

                  if (!transpose) {
                    Thyra::update(1.0,  *tmp, dgdp_out.ptr());
                  } else {
                    TEUCHOS_TEST_FOR_EXCEPTION(
                      transpose,
                      std::logic_error,
                      std::endl <<
                      "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
                      "The code related to df/dp operator and dg/dp with DERIV_MV_JACOBIAN_FORM layout has been commented out because never tested.  " <<
                      std::endl);
                  }
                }
                else {
                  TEUCHOS_TEST_FOR_EXCEPTION(
                    dfdp_op == Teuchos::null,
                    std::logic_error,
                    std::endl <<
                    "Piro::SteadyStateSolver::evalConvergedModelResponsesAndSensitivities():  " <<
                    "The code related to df/dp multivector has been commented out because never tested.  " <<
                    std::endl);
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
void Piro::SteadyStateSolver<Scalar>::evalReducedHessian(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs,
    Teuchos::ParameterList& appParams) const
{
#ifdef HAVE_PIRO_ROL
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Compute the number of directions
  int n_directions = 0;

  for (auto p_index = 0; p_index < num_p_; ++p_index) {
    if (Teuchos::nonnull(modelInArgs.get_p_direction(p_index))) {
      const int current_n_directions = modelInArgs.get_p_direction(p_index)->domain()->dim();
      if (current_n_directions > n_directions)
        n_directions = current_n_directions;
    }
  }

  if (n_directions==0)
    return;

  double tol = 1e-8;

  RCP<Thyra::VectorSpaceBase<Scalar> const> p_space = this->getModel().get_p_space(0);
  RCP<Thyra::VectorBase<Scalar>> thyra_p = Thyra::createMember(p_space);
  RCP<Thyra::VectorBase<Scalar>> thyra_direction_p = Thyra::createMember(p_space);

  RCP<const Thyra::VectorBase<Scalar> > p_init = modelInArgs.get_p(0) != Teuchos::null ? modelInArgs.get_p(0) : this->getModel().getNominalValues().get_p(0);
  Thyra::copy(*p_init, thyra_p.ptr());

  ROL::ThyraVector<Scalar> rol_p(thyra_p);
  ROL::ThyraVector<Scalar> rol_direction_p(thyra_direction_p);

  Teuchos::RCP<Thyra::VectorSpaceBase<Scalar> const> x_space = this->getModel().get_x_space();
  Teuchos::RCP<Thyra::VectorBase<Scalar>> x = Thyra::createMember(x_space);
  Thyra::copy(*modelInArgs.get_x(), x.ptr());

  ROL::ThyraVector<Scalar> rol_x(x);
  Teuchos::RCP<Thyra::VectorBase<Scalar>> lambda_vec = Thyra::createMember(x_space);
  ROL::ThyraVector<Scalar> rol_lambda(lambda_vec);

  RCP<Thyra::VectorBase<Scalar>> current_hv = Thyra::createMember(p_space);
  ROL::ThyraVector<Scalar> rol_current_hv(current_hv);

  ROL::Ptr<ROL::Vector<Scalar> > rol_p_ptr = ROL::makePtrFromRef(rol_p);
  ROL::Ptr<ROL::Vector<Scalar> > rol_x_ptr = ROL::makePtrFromRef(rol_x);
  ROL::Ptr<ROL::Vector<Scalar> > rol_lambda_ptr = ROL::makePtrFromRef(rol_lambda);

  if (Teuchos::nonnull(modelInArgs.get_f_multiplier())) {
    RCP<const Thyra::VectorBase<Scalar> > lambda_init = modelInArgs.get_f_multiplier();
    Thyra::copy(*lambda_init, lambda_vec.ptr());
  }


  Piro::ThyraProductME_Constraint_SimOpt<Scalar> constr(model_, adjointModel_, appParams, Teuchos::VERB_NONE);
  auto stateStore = ROL::makePtr<ROL::VectorController<Scalar>>(); 
  
  for (int g_index=0; g_index<num_g_; ++g_index) {
    Piro::ThyraProductME_Objective_SimOpt<Scalar> obj(model_, g_index, appParams, Teuchos::VERB_NONE);
    
    ROL::Ptr<ROL::Constraint_SimOpt<Scalar> > constr_ptr = ROL::makePtrFromRef(constr);
    ROL::Ptr<ROL::Objective_SimOpt<Scalar> > obj_ptr = ROL::makePtrFromRef(obj);    

    //create the ROL reduce objective initializing it with the current state and parameter
    ROL::Reduced_Objective_SimOpt<Scalar> reduced_obj(obj_ptr,constr_ptr,stateStore,rol_x_ptr,rol_p_ptr,rol_lambda_ptr);
    reduced_obj.update(rol_p,ROL::UpdateType::Temp);
    stateStore->set(*rol_x_ptr, std::vector<Scalar>());  //second argument not meaningful for deterministic problems 

    if(g_index>0) { //a bit hacky, but the jacobian and, if needed, its adjoint have been computed at iteration 0
      constr.computeJacobian1_ = false;
      constr.computeAdjointJacobian1_ = false;
    }

    for (auto j = 0; j < n_directions; ++j) {
      if (Teuchos::nonnull(modelInArgs.get_p_direction(0))) {
        const int current_n_directions = modelInArgs.get_p_direction(0)->domain()->dim();
        if (j < current_n_directions) {
          RCP<const Thyra::VectorBase<Scalar> > current_p_direction = modelInArgs.get_p_direction(0)->col(j);
          Thyra::copy(*current_p_direction, thyra_direction_p.ptr());
        }
      }

      reduced_obj.hessVec(rol_current_hv, rol_direction_p, rol_p, tol);

      
      if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, 0, 0)) {
        RCP<Thyra::MultiVectorBase<Scalar> > reduced_hv_out = outArgs.get_hess_vec_prod_g_pp(g_index, 0, 0);
        Thyra::assign(reduced_hv_out->col(j).ptr(), *current_hv);
      }
    }
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    std::endl <<
    "Piro::SteadyStateSolver::evalReducedHessian():  " <<
    "This fuction requires ROL." <<
    std::endl);
#endif
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
