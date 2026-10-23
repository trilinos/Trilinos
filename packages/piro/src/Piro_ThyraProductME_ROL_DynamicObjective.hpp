// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_THYRAPRODUCTME_ROL_DYNAMICOBJECTIVE_HPP
#define PIRO_THYRAPRODUCTME_ROL_DYNAMICOBJECTIVE_HPP

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
#include "ROL_DynamicObjective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_ThyraVector.hpp"

namespace Piro {

//! Backward-difference approximation of the time derivative over a time step:
//! u_dot = (u_new - u_old) / dt.
template <typename Real>
  Teuchos::RCP<Thyra::VectorBase<Real> > compute_u_dot(Teuchos::RCP<const Thyra::VectorBase<Real> > u_new, Teuchos::RCP<const Thyra::VectorBase<Real> > u_old, Real dt) {
    Teuchos::RCP<Thyra::VectorBase<Real> > u_dot = u_new->clone_v();
    Thyra::V_StVpStV(u_dot.ptr(), Real(1.0)/dt, *u_new, Real(-1.0)/dt, *u_old);
    return u_dot;
  }

/** \brief ThyraProductME_ROL_DynamicObjective
 *
 * Objective contribution of one ROL time stamp, used by ROL::ReducedDynamicObjective.
 *
 * ROL needs derivatives with respect to u_old, u_new, and z, whereas the
 * model evaluator provides derivatives with respect to x, x_dot, and p.
 * On a time stamp [t_old, t_new], with dt = t_new - t_old, the time derivative
 * is approximated by the backward difference
 * \f[
 *      \dot{u} = \frac{u_{new}-u_{old}}{dt}.
 * \f]
 *
 * If the response depends only on the final time, the last time stamp contributes
 * \f[
 *      G = g(u_{new}, \dot{u}, z, t_{new})
 * \f]
 * and all other time stamps contribute zero.
 *
 * Otherwise the response is integrated in time with the trapezoidal rule
 * \f[
 *      G = \frac{dt}{2} \left[ g(u_{new}, \dot{u}, z, t_{new}) + g(u_{old}, \dot{u}, z, t_{old}) \right].
 * \f]
 *
 * The gradients are the exact derivatives of G, including the dependence of
 * \f$\dot{u}\f$ on both u_old (factor -1/dt) and u_new (factor 1/dt). Models that
 * do not support DgDx_dot are treated as responses that do not depend on x_dot.
*/

template <typename Real>
class ThyraProductME_ROL_DynamicObjective : public virtual ROL::DynamicObjective<Real> {
public:

  ThyraProductME_ROL_DynamicObjective(
    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & model,
    const Teuchos::RCP<Tempus::Integrator<Real> >& integrator,
    const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & modelAdjoin,
    int g_index,
    Teuchos::ParameterList& piroParams,
    const int Nt,
    const bool onlyFinalTime = true,
    const bool useTrapezoidalTimeIntegration = true,
    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_HIGH,
    Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null);

  virtual ~ThyraProductME_ROL_DynamicObjective() {}

  //! Compute value of objective
  Real value( const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
              const ROL::Vector<Real> &z, const ROL::TimeStamp<Real> &timeStamp ) const;

  void gradient_uo( ROL::Vector<Real> &g, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                    const ROL::Vector<Real> &z, const ROL::TimeStamp<Real> &timeStamp ) const;

  void gradient_un( ROL::Vector<Real> &g, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                    const ROL::Vector<Real> &z, const ROL::TimeStamp<Real> &timeStamp ) const;

  void gradient_z( ROL::Vector<Real> &g, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                    const ROL::Vector<Real> &z, const ROL::TimeStamp<Real> &timeStamp ) const;

  //! Response g(x, x_dot, p, t)
  Real value( const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot,
              const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t) const;

  //! Partial derivative dg/dp(x, x_dot, p, t)
  void gradient_z( Teuchos::RCP<Thyra::VectorBase<Real>> &g, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot,
                    const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t) const;

  using ROL::DynamicObjective<Real>::update;

  void update( const ROL::Vector<Real> &x, ROL::UpdateType type, int iter = -1 ) {
    (void) x;
    (void) type;
    (void) iter;
  }

  // Update old state
  void update_uo( const ROL::Vector<Real> &x, const ROL::TimeStamp<Real> &ts ) {
    (void) x;
    (void) ts;
  }

  // Update new state
  void update_un( const ROL::Vector<Real> &x, const ROL::TimeStamp<Real> &ts ) {
    (void) x;
    (void) ts;
  }

  // Update control
  void update_z( const ROL::Vector<Real> &x, const ROL::TimeStamp<Real> &ts ) {
    (void) x;
    (void) ts;
  }

private:
  //! Whether the time stamp contributes to the objective.
  bool isActive(const ROL::TimeStamp<Real> &timeStamp) const {
    return !onlyFinalTime_ || (int) timeStamp.k == Nt_-1;
  }

  //! Partial derivatives dg/dx and dg/dx_dot (gradient form) at (x, x_dot, p, t).
  //! dgdx_dot is set to zero when the model does not support DgDx_dot.
  void evalStateDerivatives( const Teuchos::RCP<Thyra::VectorBase<Real>> &dgdx,
                             const Teuchos::RCP<Thyra::VectorBase<Real>> &dgdx_dot,
                             const Teuchos::RCP<const Thyra::VectorBase<Real>> &x,
                             const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot,
                             const Teuchos::RCP<const Thyra::VectorBase<Real>> &p,
                             const Real &t) const;

  const Teuchos::RCP<Tempus::Integrator<Real> > integrator_;
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_model_;
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_adjoint_model_;
  std::string sensitivity_method_;
  const int g_index_;
  int Nt_;
  const bool onlyFinalTime_;
  const bool useTrapezoidalTimeIntegration_;
  bool supportsDgDxDot_;
  Real objectiveRecoveryValue_;
  bool useObjectiveRecoveryValue_;
  ROL::UpdateType updateType_;

  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;
  Teuchos::RCP<ROL_ObserverBase<Real>> observer_;

  Teuchos::RCP<Teuchos::ParameterList> tempus_params_;

}; // class ThyraProductME_ROL_DynamicObjective

template <typename Real>
ThyraProductME_ROL_DynamicObjective<Real>::
ThyraProductME_ROL_DynamicObjective(
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & model,
  const Teuchos::RCP<Tempus::Integrator<Real> >& integrator,
  const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & modelAdjoin,
  int g_index,
  Teuchos::ParameterList& piroParams,
  const int Nt,
  const bool onlyFinalTime,
  const bool useTrapezoidalTimeIntegration,
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> observer) :
  integrator_(integrator),
  thyra_model_(model),
  thyra_adjoint_model_(modelAdjoin),
  sensitivity_method_("Adjoint"),
  g_index_(g_index),
  Nt_(Nt),
  onlyFinalTime_(onlyFinalTime),
  useTrapezoidalTimeIntegration_(useTrapezoidalTimeIntegration),
  supportsDgDxDot_(false),
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel),
  observer_(observer),
  tempus_params_(Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList(piroParams.sublist("Tempus"))))
{
  typedef Thyra::ModelEvaluatorBase MEB;

  TEUCHOS_TEST_FOR_EXCEPTION(!onlyFinalTime_ && !useTrapezoidalTimeIntegration_, std::logic_error,
    std::endl << "Piro::ThyraProductME_ROL_DynamicObjective: time-integrated responses are "
    "only implemented with the trapezoidal rule." << std::endl);

  TEUCHOS_TEST_FOR_EXCEPTION(thyra_model_->get_g_space(g_index_)->dim() != 1, Teuchos::Exceptions::InvalidParameter,
    std::endl << "Piro::ThyraProductME_ROL_DynamicObjective: the objective response " << g_index_ <<
    " has dimension " << thyra_model_->get_g_space(g_index_)->dim() << ", but it must be a scalar." << std::endl);

  const MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

  TEUCHOS_TEST_FOR_EXCEPTION(
    !outArgs.supports(MEB::OUT_ARG_DgDx, g_index_).supports(MEB::DERIV_MV_GRADIENT_FORM),
    std::logic_error,
    std::endl << "Piro::ThyraProductME_ROL_DynamicObjective: DgDx of response " << g_index_ <<
    " must be supported in gradient form (DERIV_MV_GRADIENT_FORM)." << std::endl);

  const MEB::DerivativeSupport dgdx_dot_support = outArgs.supports(MEB::OUT_ARG_DgDx_dot, g_index_);
  supportsDgDxDot_ = dgdx_dot_support.supports(MEB::DERIV_MV_GRADIENT_FORM);
  if (!supportsDgDxDot_ && !dgdx_dot_support.none()) {
    *out_ << "WARNING: Piro::ThyraProductME_ROL_DynamicObjective: DgDx_dot of response " << g_index_
          << " is supported, but not in gradient form (DERIV_MV_GRADIENT_FORM). The dependence "
          << "of the response on x_dot is ignored when computing gradients." << std::endl;
  }
}

template <typename Real>
Real
ThyraProductME_ROL_DynamicObjective<Real>::
value( const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
              const ROL::Vector<Real> &p, const ROL::TimeStamp<Real> &timeStamp ) const
{
  using Teuchos::RCP;

  const Real t_old = timeStamp.t.front();
  const Real t_new = timeStamp.t.back();

  if(!isActive(timeStamp)) {
    if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
      *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::value final time of the time stamp " << t_new << " is not the final time." << std::endl;
    return 0;
  }
  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::value time stamp [" << t_old << ", " << t_new << "]" << std::endl;

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  const Real dt = t_new - t_old;

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  Real g = value( thyra_u_new.getVector(), u_dot, thyra_p.getVector(), t_new);

  if (!onlyFinalTime_) {
    g += value( thyra_u_old.getVector(), u_dot, thyra_p.getVector(), t_old);
    g *= dt/2.;
  }

  return g;
}

template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
gradient_uo( ROL::Vector<Real> &grad, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
              const ROL::Vector<Real> &p, const ROL::TimeStamp<Real> &timeStamp ) const
{
  using Teuchos::RCP;

  RCP<Thyra::VectorBase<Real> > dgduo = Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector();

  const Real t_old = timeStamp.t.front();
  const Real t_new = timeStamp.t.back();

  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_uo " << t_old << " " << t_new << " " << timeStamp.k << " " << Nt_ << std::endl;

  // In the final-time case u_old enters only through u_dot.
  if(!isActive(timeStamp) || (onlyFinalTime_ && !supportsDgDxDot_)) {
    Thyra::assign(dgduo.ptr(), Teuchos::ScalarTraits<Real>::zero());
    return;
  }

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  const Real dt = t_new - t_old;
  const Real one(1.0), half(0.5);

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  RCP<Thyra::VectorBase<Real> > dgdx_new = dgduo->clone_v();
  RCP<Thyra::VectorBase<Real> > dgdx_dot_new = dgduo->clone_v();
  evalStateDerivatives(dgdx_new, dgdx_dot_new, thyra_u_new.getVector(), u_dot, thyra_p.getVector(), t_new);

  if (onlyFinalTime_) {
    // G = g(u_new, u_dot)  ==>  dG/du_old = -1/dt dg/dx_dot(new)
    Thyra::V_StV(dgduo.ptr(), -one/dt, *dgdx_dot_new);
  }
  else {
    // G = dt/2 [g(u_new, u_dot) + g(u_old, u_dot)]
    // ==>  dG/du_old = dt/2 dg/dx(old) - 1/2 [dg/dx_dot(new) + dg/dx_dot(old)]
    RCP<Thyra::VectorBase<Real> > dgdx_old = dgduo->clone_v();
    RCP<Thyra::VectorBase<Real> > dgdx_dot_old = dgduo->clone_v();
    evalStateDerivatives(dgdx_old, dgdx_dot_old, thyra_u_old.getVector(), u_dot, thyra_p.getVector(), t_old);
    Thyra::V_StVpStV(dgduo.ptr(), dt/2., *dgdx_old, -half, *dgdx_dot_new);
    Thyra::Vp_StV(dgduo.ptr(), -half, *dgdx_dot_old);
  }
}

template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
gradient_un( ROL::Vector<Real> &grad, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
              const ROL::Vector<Real> &p, const ROL::TimeStamp<Real> &timeStamp ) const
{
  using Teuchos::RCP;

  RCP<Thyra::VectorBase<Real> > dgdun = Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector();

  const Real t_old = timeStamp.t.front();
  const Real t_new = timeStamp.t.back();

  if(!isActive(timeStamp)) {
    if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
      *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_un final time of the time stamp " << t_new << " is not the final time." << std::endl;
    Thyra::assign(dgdun.ptr(), Teuchos::ScalarTraits<Real>::zero());
    return;
  }
  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_un time stamp [" << t_old << ", " << t_new << "]" << std::endl;

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  const Real dt = t_new - t_old;
  const Real one(1.0), half(0.5);

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  RCP<Thyra::VectorBase<Real> > dgdx_new = dgdun->clone_v();
  RCP<Thyra::VectorBase<Real> > dgdx_dot_new = dgdun->clone_v();
  evalStateDerivatives(dgdx_new, dgdx_dot_new, thyra_u_new.getVector(), u_dot, thyra_p.getVector(), t_new);

  if (onlyFinalTime_) {
    // G = g(u_new, u_dot)  ==>  dG/du_new = dg/dx(new) + 1/dt dg/dx_dot(new)
    Thyra::V_StVpStV(dgdun.ptr(), one, *dgdx_new, one/dt, *dgdx_dot_new);
  }
  else {
    // G = dt/2 [g(u_new, u_dot) + g(u_old, u_dot)]
    // ==>  dG/du_new = dt/2 dg/dx(new) + 1/2 [dg/dx_dot(new) + dg/dx_dot(old)]
    RCP<Thyra::VectorBase<Real> > dgdx_old = dgdun->clone_v();
    RCP<Thyra::VectorBase<Real> > dgdx_dot_old = dgdun->clone_v();
    evalStateDerivatives(dgdx_old, dgdx_dot_old, thyra_u_old.getVector(), u_dot, thyra_p.getVector(), t_old);
    Thyra::V_StVpStV(dgdun.ptr(), dt/2., *dgdx_new, half, *dgdx_dot_new);
    Thyra::Vp_StV(dgdun.ptr(), half, *dgdx_dot_old);
  }
}

template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
gradient_z( ROL::Vector<Real> &grad, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
              const ROL::Vector<Real> &p, const ROL::TimeStamp<Real> &timeStamp ) const
{
  using Teuchos::RCP;

  const Real t_old = timeStamp.t.front();
  const Real t_new = timeStamp.t.back();

  if(!isActive(timeStamp)) {
    if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
      *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_z final time of the time stamp " << t_new << " is not the final time." << std::endl;
    Thyra::assign(Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector().ptr(), Teuchos::ScalarTraits<Real>::zero());
    return;
  }
  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_z time stamp [" << t_old << ", " << t_new << "]" << std::endl;

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  ROL::ThyraVector<Real>  & thyra_dgdp = dynamic_cast<ROL::ThyraVector<Real>&>(grad);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  const Real dt = t_new - t_old;

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  RCP<Thyra::VectorBase<Real> > dgdp = thyra_dgdp.getVector();
  gradient_z( dgdp, thyra_u_new.getVector(), u_dot, thyra_p.getVector(), t_new);

  if (!onlyFinalTime_) {
    // G = dt/2 [g(u_new, u_dot) + g(u_old, u_dot)]  ==>  dG/dz = dt/2 [dg/dp(new) + dg/dp(old)]
    RCP<Thyra::VectorBase<Real> > dgdp_old = thyra_dgdp.getVector()->clone_v();
    gradient_z( dgdp_old, thyra_u_old.getVector(), u_dot, thyra_p.getVector(), t_old);
    Thyra::Vp_V(dgdp.ptr(), *dgdp_old);
    Thyra::Vt_S(dgdp.ptr(), dt/2.);
  }
}

template <typename Real>
Real
ThyraProductME_ROL_DynamicObjective<Real>::
value( const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot,
              const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t) const
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  // Run tempus and compute response for specified parameter values
  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

  inArgs.set_p(0, p);
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(t);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(x_dot);

  Teuchos::RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));
  outArgs.set_g(g_index_, g);

  thyra_model_->evalModel(inArgs, outArgs);

  return ::Thyra::get_ele(*g,0);
}

template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
evalStateDerivatives( const Teuchos::RCP<Thyra::VectorBase<Real>> &dgdx,
                      const Teuchos::RCP<Thyra::VectorBase<Real>> &dgdx_dot,
                      const Teuchos::RCP<const Thyra::VectorBase<Real>> &x,
                      const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot,
                      const Teuchos::RCP<const Thyra::VectorBase<Real>> &p,
                      const Real &t) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

  inArgs.set_p(0,p);
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(t);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(x_dot);

  outArgs.set_DgDx(g_index_, MEB::DerivativeMultiVector<Real>(dgdx, MEB::DERIV_MV_GRADIENT_FORM));
  if (supportsDgDxDot_)
    outArgs.set_DgDx_dot(g_index_, MEB::DerivativeMultiVector<Real>(dgdx_dot, MEB::DERIV_MV_GRADIENT_FORM));

  thyra_model_->evalModel(inArgs, outArgs);

  if (!supportsDgDxDot_)
    Thyra::assign(dgdx_dot.ptr(), Teuchos::ScalarTraits<Real>::zero());
}


template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
gradient_z( Teuchos::RCP<Thyra::VectorBase<Real>> &grad, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot, 
              const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t) const
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

  inArgs.set_p(0,p);
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(t);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(x_dot);

  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));

  outArgs.set_g(g_index_, g);

  Teuchos::RCP<Thyra::ProductMultiVectorBase<Real> > prodvec_dgdp =
      Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(grad);
  if ( !g.is_null()) {
    if ( !prodvec_dgdp.is_null()) {
      Teuchos::RCP<const Piro::ProductModelEvaluator<Real>> model_PME = getProductModelEvaluator(thyra_model_);

      if ( !model_PME.is_null()) {
        Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> dgdp_op =
            Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(model_PME->create_DgDp_op(g_index_, 0, prodvec_dgdp));
        Thyra::ModelEvaluatorBase::Derivative<Real> dgdp_der(Teuchos::rcp_dynamic_cast<Thyra::LinearOpBase<Real>>(dgdp_op));
        outArgs.set_DgDp(g_index_, 0, dgdp_der);
      }
      else {
        ROL_TEST_FOR_EXCEPTION( true, std::logic_error, "Piro::ThyraProductME_ROL_DynamicObjective::gradient_z: dgdp is not supported for the used ModelEvaluator.");
      }
    }
    else {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
          outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, g_index_, 0);
      Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient;
      if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
        dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
      else if(dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
        dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
      else {
        ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
            "Piro::ThyraProductME_ROL_DynamicObjective::gradient_z: DgDp does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
      }
      outArgs.set_DgDp(g_index_, 0, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(grad, dgdp_orient));
    }
  }

  thyra_model_->evalModel(inArgs, outArgs);
}

} // namespace Piro

#endif
