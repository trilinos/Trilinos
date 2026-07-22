// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TEMPUS_REDUCED_OBJECTIVE_HPP
#define ROL_TEMPUS_REDUCED_OBJECTIVE_HPP

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorForwardSensitivity.hpp"
#include "Tempus_IntegratorAdjointSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientForwardSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientAdjointSensitivity.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_ThyraVector.hpp"

namespace ROL {

template <typename Real>
class TempusReducedObjective : public virtual ROL::Objective<Real> {
public:

  TempusReducedObjective(
    const Teuchos::RCP<Thyra::ModelEvaluator<Real> >& model,
    const Teuchos::RCP<Teuchos::ParameterList>& tempus_params,
    const Teuchos::RCP<Teuchos::ParameterList>& objective_params);

  virtual ~TempusReducedObjective() {}

  //! Compute value of objective
  Real value( const ROL::Vector<Real> &x, Real &tol );

  //! Compute gradient of objective
  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol );

  //! Set response target for computing objective
  void set_target(const Teuchos::RCP<Thyra::VectorBase<Real> >& target);
  void set_target(const Teuchos::RCP<ROL::Vector<Real> >& target);

  //! Helper function to create optimization vector
  Teuchos::RCP<ROL::Vector<Real> > create_design_vector() const;

  //! Helper function to create a response vector
  Teuchos::RCP<ROL::Vector<Real> > create_response_vector() const;

  //! Helper function to run tempus, computing responses and derivatives
  void run_tempus(ROL::Vector<Real>& r, const ROL::Vector<Real>& p);
  void run_tempus(const Thyra::ModelEvaluatorBase::InArgs<Real>&  inArgs,
                  const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs);

private:

  Teuchos::RCP<Thyra::ModelEvaluator<Real> > model_;
  Teuchos::RCP<Teuchos::ParameterList> tempus_params_;
  Teuchos::RCP<Thyra::VectorBase<Real> > target_;
  std::string objective_type_;
  std::string sensitivity_method_;
  int param_index_;
  int response_index_;
  bool use_fd_gradient_;

}; // class TempusReducedObjective

template <typename Real>
TempusReducedObjective<Real>::
TempusReducedObjective(
  const Teuchos::RCP<Thyra::ModelEvaluator<Real> >& model,
  const Teuchos::RCP<Teuchos::ParameterList>& tempus_params,
  const Teuchos::RCP<Teuchos::ParameterList>& objective_params) :
  model_(model),
  tempus_params_(tempus_params),
  objective_type_("Sum of Squares"),
  sensitivity_method_("Forward"),
  param_index_(0),
  response_index_(0),
  use_fd_gradient_(false)
{
  objective_type_ =
    objective_params->get("Objective Type", objective_type_);
  sensitivity_method_ =
    objective_params->get("Sensitivity Method", sensitivity_method_);
  param_index_ =
    objective_params->get("Sensitivity Parameter Index", param_index_);
  response_index_ =
    objective_params->get("Response Function Index", response_index_);
  use_fd_gradient_ =
    objective_params->get("Use FD Gradient", use_fd_gradient_);
}

template <typename Real>
Real
TempusReducedObjective<Real>::
value( const ROL::Vector<Real> &p, Real &tol )
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  // Run tempus and compute response for specified parameter values
  MEB::InArgs<Real> inArgs = model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  inArgs.set_p(param_index_, thyra_p.getVector());
  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(model_->get_g_space(response_index_));
  outArgs.set_g(response_index_, g);
  run_tempus(inArgs, outArgs);

  // Compute objective
  if (target_ != Teuchos::null)
    Thyra::Vp_StV(g.ptr(), -1.0, *target_); // g <- g - target
  Real val = 0.0;
  if (objective_type_ == "Sum of Squares") {
    Thyra::ele_wise_scale(*g, g.ptr()); // g <- g.*g
    val = Thyra::sum(*g); // val = \sum_i g_i
  }
  else if (objective_type_ == "Sum")
    val = Thyra::sum(*g); // val = \sum_i g_i
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Invalid objective type " << objective_type_ << "!" << std::endl
      << "Valid choices are:  \"Sum\" and \"Sum of Squares\".");

  return val;
}

template <typename Real>
void
TempusReducedObjective<Real>::
gradient(ROL::Vector<Real> &grad, const ROL::Vector<Real> &p, Real &tol)
{
  if (use_fd_gradient_) {
    ROL::Objective<Real>::gradient(grad, p, tol);
    return;
  }

  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  // Run tempus and compute response gradient for specified parameter values
  const int num_p = model_->get_p_space(param_index_)->dim();
  const int num_g = model_->get_g_space(response_index_)->dim();
  MEB::InArgs<Real> inArgs = model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  inArgs.set_p(param_index_, thyra_p.getVector());
  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(model_->get_g_space(response_index_));
  RCP<Thyra::MultiVectorBase<Real> > dgdp;
  MEB::EDerivativeMultiVectorOrientation dgdp_orientation =
    MEB::DERIV_MV_JACOBIAN_FORM;
  if (sensitivity_method_ == "Adjoint"||
      sensitivity_method_ == "Pseudotransient Adjoint") {
    // Adjoint, PseudoTransientAdjoint integrators require gradient form
    // (but does not necessarily require the model to support gradient form)
    dgdp =
      Thyra::createMembers<Real>(model_->get_p_space(param_index_), num_g);
    dgdp_orientation = MEB::DERIV_MV_GRADIENT_FORM;
  }
  else {
    // Form determined by what model supports
    MEB::DerivativeSupport support =
      outArgs.supports(MEB::OUT_ARG_DgDp, response_index_, param_index_);
    if (support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
      dgdp =
        Thyra::createMembers<Real>(model_->get_g_space(response_index_), num_p);
      dgdp_orientation = MEB::DERIV_MV_JACOBIAN_FORM;
    }
    else if (support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
      dgdp =
        Thyra::createMembers<Real>(model_->get_p_space(param_index_), num_g);
      dgdp_orientation = MEB::DERIV_MV_GRADIENT_FORM;
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Model must support Jacobian or gradient forms for dg/dp!");
  }
  outArgs.set_DgDp(response_index_, param_index_,
                   MEB::DerivativeMultiVector<Real>(dgdp, dgdp_orientation));
  outArgs.set_g(response_index_, g);
  run_tempus(inArgs, outArgs);

  // Compute objective gradient
  RCP<Thyra::VectorBase<Real> > final_grad =
    Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector()->col(0);
  if (target_ != Teuchos::null)
    Thyra::Vp_StV(g.ptr(), -1.0, *target_); // g <- g - target
  if (dgdp_orientation == MEB::DERIV_MV_JACOBIAN_FORM) {
    for (int j=0; j<num_p; ++j) {
      RCP<Thyra::VectorBase<Real> > dgdp_j = dgdp->col(j);
      Real grad_val = 0.0;
      if (objective_type_ == "Sum of Squares") {
        Thyra::ele_wise_prod_update(2.0, *g, dgdp_j.ptr());
        grad_val = Thyra::sum(*dgdp_j);
      }
      else if (objective_type_ == "Sum")
        grad_val = Thyra::sum(*dgdp_j);
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Invalid objective type " << objective_type_ << "!" << std::endl
          << "Valid choices are:  \"Sum\" and \"Sum of Squares\".");
      Thyra::set_ele(j, grad_val, final_grad.ptr());
    }
  }
  else { // gradient form
    Thyra::assign(final_grad.ptr(), 0.0);
    for (int i=0; i<num_g; ++i) {
      if (objective_type_ == "Sum of Squares") {
        Real g_val = Thyra::get_ele(*g, i);
        Thyra::Vp_StV(final_grad.ptr(), 2.0*g_val, *(dgdp->col(i)));
      }
      else if (objective_type_ == "Sum") {
        Thyra::Vp_StV(final_grad.ptr(), 1.0, *(dgdp->col(i)));
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Invalid objective type " << objective_type_ << "!" << std::endl
          << "Valid choices are:  \"Sum\" and \"Sum of Squares\".");
    }
  }
}

template <typename Real>
void
TempusReducedObjective<Real>::
set_target(const Teuchos::RCP<Thyra::VectorBase<Real> >& target) {
  target_ = target;
}

template <typename Real>
void
TempusReducedObjective<Real>::
set_target(const Teuchos::RCP<ROL::Vector<Real> >& target) {
  using Teuchos::rcp_dynamic_cast;
  target_ =
    rcp_dynamic_cast<ROL::ThyraVector<double> >(target, true)->getVector();
}

template <typename Real>
Teuchos::RCP<ROL::Vector<Real> >
TempusReducedObjective<Real>::
create_design_vector() const {
  typedef Thyra::ModelEvaluatorBase MEB;

  Teuchos::RCP<Thyra::VectorBase<Real> > p =
    Thyra::createMember<Real>(model_->get_p_space(param_index_));
  MEB::InArgs<Real> nominalValues = model_->getNominalValues();
  if (nominalValues.get_p(param_index_) != Teuchos::null)
    Thyra::assign(p.ptr(), *(nominalValues.get_p(param_index_)));
  else
    Thyra::assign(p.ptr(), Teuchos::ScalarTraits<Real>::zero());
  return Teuchos::rcp(new ROL::ThyraVector<Real>(p));
}

template <typename Real>
Teuchos::RCP<ROL::Vector<Real> >
TempusReducedObjective<Real>::
create_response_vector() const {
  Teuchos::RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(model_->get_g_space(response_index_));
  Thyra::assign(g.ptr(), Teuchos::ScalarTraits<Real>::zero());
  return Teuchos::rcp(new ROL::ThyraVector<Real>(g));
}

template <typename Real>
void
TempusReducedObjective<Real>::
run_tempus(ROL::Vector<Real>& r, const ROL::Vector<Real>& p)
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Real> inArgs = model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  ROL::ThyraVector<Real>& thyra_r =
    Teuchos::dyn_cast<ROL::ThyraVector<Real> >(r);
  inArgs.set_p(param_index_, thyra_p.getVector());
  outArgs.set_g(response_index_, thyra_r.getVector());
  run_tempus(inArgs, outArgs);
}

template <typename Real>
void
TempusReducedObjective<Real>::
run_tempus(const Thyra::ModelEvaluatorBase::InArgs<Real>&  inArgs,
           const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs)
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultNominalBoundsOverrideModelEvaluator<Real> DNBOME;

  // Override nominal values in model to supplied inArgs
  RCP<DNBOME> wrapped_model = rcp(new DNBOME(model_, rcpFromRef(inArgs)));

  Real t;
  RCP<const Thyra::VectorBase<Real> > x, x_dot;
  RCP<const Thyra::MultiVectorBase<double> > dxdp, dxdotdp;
  RCP<Thyra::VectorBase<Real> > g = outArgs.get_g(response_index_);
  RCP<Thyra::MultiVectorBase<Real> > dgdp =
    outArgs.get_DgDp(response_index_, param_index_).getMultiVector();
  MEB::EDerivativeMultiVectorOrientation dgdp_orientation =
    outArgs.get_DgDp(response_index_, param_index_).getMultiVectorOrientation();

  // Create and run integrator
  if (dgdp != Teuchos::null && sensitivity_method_ == "Forward") {
    RCP<Tempus::IntegratorForwardSensitivity<Real> > integrator =
      Tempus::createIntegratorForwardSensitivity<Real>(tempus_params_, wrapped_model);
    const bool integratorStatus = integrator->advanceTime();
    TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorStatus, std::logic_error, "Integrator failed!");

    // Get final state
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
    dxdp = integrator->getDxDp();
    dxdotdp = integrator->getDXDotDp();
  }
  else if (dgdp != Teuchos::null && sensitivity_method_ == "Adjoint") {
    RCP<Tempus::IntegratorAdjointSensitivity<Real> > integrator =
      Tempus::createIntegratorAdjointSensitivity<Real>(tempus_params_, wrapped_model);
    const bool integratorStatus = integrator->advanceTime();
    TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorStatus, std::logic_error, "Integrator failed!");

    // Get final state
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
    Thyra::assign(dgdp.ptr(), *(integrator->getDgDp()));
  }
  else if (dgdp != Teuchos::null &&
           sensitivity_method_ == "Pseudotransient Forward") {
    RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Real> > integrator =
      Tempus::createIntegratorPseudoTransientForwardSensitivity<Real>(tempus_params_,
                                                                wrapped_model);
    const bool integratorStatus = integrator->advanceTime();
    TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorStatus, std::logic_error, "Integrator failed!");

    // Get final state
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
    dxdp = integrator->getDxDp();
    dxdotdp = integrator->getDXDotDp();
  }
  else if (dgdp != Teuchos::null &&
           sensitivity_method_ == "Pseudotransient Adjoint") {
    RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Real> > integrator =
      Tempus::integratorPseudoTransientAdjointSensitivity<Real>(tempus_params_,
                                                                wrapped_model);
    const bool integratorStatus = integrator->advanceTime();
    TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorStatus, std::logic_error, "Integrator failed!");

    // Get final state
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
    Thyra::assign(dgdp.ptr(), *(integrator->getDgDp()));
  }
  else if (dgdp != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Invalid sensitivity method " << sensitivity_method_ << "!" << std::endl
      << "Valid choices are:  Forward, Adjoint, Pseudotransient Forward, "
      << " and Pseudotransient Adjoint.");
  }
  else {
    RCP<Tempus::IntegratorBasic<Real> > integrator =
      Tempus::createIntegratorBasic<Real>(tempus_params_, wrapped_model);
    const bool integratorStatus = integrator->advanceTime();
    TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorStatus, std::logic_error, "Integrator failed!");

    // Get final state
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
  }

  // Evaluate response at final state
  const int num_g = model_->get_g_space(response_index_)->dim();
  MEB::InArgs<Real> modelInArgs   = inArgs;
  MEB::OutArgs<Real> modelOutArgs = outArgs;
  modelInArgs.set_x(x);
  if (modelInArgs.supports(MEB::IN_ARG_x_dot)) modelInArgs.set_x_dot(x_dot);
  if (modelInArgs.supports(MEB::IN_ARG_t)) modelInArgs.set_t(t);
  RCP<Thyra::MultiVectorBase<Real> > dgdx, dgdxdot;
  MEB::EDerivativeMultiVectorOrientation dgdx_orientation =
    MEB::DERIV_MV_JACOBIAN_FORM;
  MEB::EDerivativeMultiVectorOrientation dgdxdot_orientation =
    MEB::DERIV_MV_JACOBIAN_FORM;
  if (dgdp != Teuchos::null &&
      (sensitivity_method_ == "Forward" ||
       sensitivity_method_ == "Pseudotransient Forward")) {
    MEB::DerivativeSupport dgdx_support =
      outArgs.supports(MEB::OUT_ARG_DgDx, response_index_);
    if (!dgdx_support.none()) {
      if (dgdx_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
        dgdx = Thyra::createMembers<Real>(model_->get_x_space(), num_g);
        dgdx_orientation = MEB::DERIV_MV_GRADIENT_FORM;
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Model must support gradient forms for dg/dx!");
      modelOutArgs.set_DgDx(
        response_index_,
        MEB::DerivativeMultiVector<Real>(dgdx, dgdx_orientation));
    }

    MEB::DerivativeSupport dgdxdot_support =
      modelOutArgs.supports(MEB::OUT_ARG_DgDx_dot, response_index_);
    if (!dgdxdot_support.none()) {
      if (dgdxdot_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
        dgdxdot = Thyra::createMembers<Real>(model_->get_x_space(), num_g);
        dgdxdot_orientation = MEB::DERIV_MV_GRADIENT_FORM;
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Model must support Jacobian or gradient forms for dg/dx!");
      modelOutArgs.set_DgDx_dot(
        response_index_,
        MEB::DerivativeMultiVector<Real>(dgdxdot, dgdxdot_orientation));
    }
  }
  else if (dgdp != Teuchos::null &&
           (sensitivity_method_ == "Adjoint" ||
            sensitivity_method_ == "Pseudotransient Adjoint")) {
    // Clear dg/dp as an out arg since it was already computed by the adjoint
    // integrator
    modelOutArgs.set_DgDp(response_index_, param_index_,
                          MEB::Derivative<Real>());
  }

  model_->evalModel(modelInArgs, modelOutArgs);

  // dg/dp = dg/dp + dg/dx*dx/dp + dg/dx_dot*dx_dot/dp
  // We assume dg/dx, dg/dxdot are in gradient form while dxdp, dxdotdp are in
  // Jacobian form
  if (dgdp != Teuchos::null && dgdx != Teuchos::null) {
    if (dgdp_orientation == MEB::DERIV_MV_JACOBIAN_FORM)
      dgdx->apply(Thyra::TRANS, *dxdp, dgdp.ptr(), Real(1.0), Real(1.0));
    else
      dxdp->apply(Thyra::TRANS, *dgdx, dgdp.ptr(), Real(1.0), Real(1.0));
  }
  if (dgdp != Teuchos::null && dgdxdot != Teuchos::null) {
    if (dgdp_orientation == MEB::DERIV_MV_JACOBIAN_FORM)
      dgdxdot->apply(Thyra::TRANS, *dxdotdp, dgdp.ptr(), Real(1.0), Real(1.0));
    else
      dxdotdp->apply(Thyra::TRANS, *dgdxdot, dgdp.ptr(), Real(1.0), Real(1.0));
  }
}

} // namespace ROL

#endif // ROL_TEMPUS_REDUCED_OBJECTIVE_HPP
