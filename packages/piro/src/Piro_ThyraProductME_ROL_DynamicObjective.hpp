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

template <typename Real>
  Teuchos::RCP<Thyra::VectorBase<Real> > compute_u_dot(Teuchos::RCP<const Thyra::VectorBase<Real> > u_new, Teuchos::RCP<const Thyra::VectorBase<Real> > u_old, Real dt) {
    Real minus_one = -1.;
    Teuchos::RCP<Thyra::VectorBase<Real> > u_dot = u_new->clone_v();
    Thyra::Vp_V(u_dot.ptr(), *u_old, minus_one);
    Thyra::Vt_S(u_dot.ptr(), dt);
    return u_dot;
  }

/** \brief ThyraProductME_ROL_DynamicObjective
 * 
 * This class is used to be able to call Tempus from ROL in the context
 * of time integrated responses.
 * 
 * ROL needs to compute derivative with respect to u_old, u_new, and z.
 * However, Tempus allows to compute derivative with respect to x, x_dot, and p.
 * 
 * In this class member functions, we used the approximation:
 * \f[
 *      \dot{x} = \frac{u_{new}-u_{old}}{dt}
 * \f]
 * to evaluate the response value (if it depends on x_dot) and the derivatives.  
 *
 * The time integration over a ROL time stamp is done here using a trapezoidal
 * approach if the response does not depend only on the final time step.
 * 
 * The implemented trapezoidal approach visits both the first and the last time steps
 * of the ROL time stamp and evaluates the time integrand for both of these time values.
 * For the first one, x is set to u_old, x_dot is computed as above, and the time t of 
 * the first time step of the time stamp is used.
 * For the last one, x is set to u_new, x_dot is computed as above, and the time t of 
 * the last time step of the time stamp is used.
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

  Real value( const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot, 
              const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t) const;

  void gradient_uo( Teuchos::RCP<Thyra::VectorBase<Real>> &g, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot, 
                    const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t, const bool u_new, const Real &dt) const;

  void gradient_un( Teuchos::RCP<Thyra::VectorBase<Real>> &g, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot, 
                    const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t, const bool u_new, const Real &dt) const;

  void gradient_z( Teuchos::RCP<Thyra::VectorBase<Real>> &g, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot, 
                    const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t) const;

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
  const Teuchos::RCP<Tempus::Integrator<Real> > integrator_;
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_model_;
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_adjoint_model_;
  std::string sensitivity_method_;
  const int g_index_;
  int Nt_;
  const bool onlyFinalTime_;
  const bool useTrapezoidalTimeIntegration_;
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
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel),
  observer_(observer),
  tempus_params_(Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList(piroParams.sublist("Tempus"))))
{
  
}

template <typename Real>
Real
ThyraProductME_ROL_DynamicObjective<Real>::
value( const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new, 
              const ROL::Vector<Real> &p, const ROL::TimeStamp<Real> &timeStamp ) const
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  if(onlyFinalTime_ && (int) timeStamp.k != Nt_-1) {
    if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
      *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::value final time of the time stamp " << timeStamp.t[timeStamp.t.size()-1] << " is not the final time." << std::endl;
    return 0;
  }
  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::value final time of the time stamp " << timeStamp.t[timeStamp.t.size()-1] << " is the final time." << std::endl;

  TEUCHOS_TEST_FOR_EXCEPTION(!onlyFinalTime_ && !useTrapezoidalTimeIntegration_, std::logic_error,
    std::endl <<
    "Not implemented yet" << std::endl);

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  Real dt = timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0];

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  Real g = value( thyra_u_new.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[timeStamp.t.size()-1]);

  if (!onlyFinalTime_ && useTrapezoidalTimeIntegration_) {
    g += value( thyra_u_old.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[0]);
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
  typedef Thyra::ModelEvaluatorBase MEB;

  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_uo " << timeStamp.t[0] << " " << timeStamp.t[timeStamp.t.size()-1] << " " << timeStamp.k << " " << Nt_ << std::endl;
  
  if(onlyFinalTime_ || timeStamp.k <= 0) {
    Thyra::assign(Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector().ptr(), Teuchos::ScalarTraits<Real>::zero());
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(!onlyFinalTime_ && !useTrapezoidalTimeIntegration_, std::logic_error,
    std::endl <<
    "Not implemented yet" << std::endl);

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  ROL::ThyraVector<Real>  & thyra_dgdx = dynamic_cast<ROL::ThyraVector<Real>&>(grad);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  Real dt = timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0];

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  RCP<Thyra::VectorBase<Real> > dgdx = thyra_dgdx.getVector();
  gradient_uo( dgdx, thyra_u_new.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[timeStamp.t.size()-1], true, dt);

  if (!onlyFinalTime_ && useTrapezoidalTimeIntegration_) {
    RCP<Thyra::VectorBase<Real> > dgdx_old = thyra_dgdx.getVector()->clone_v();
    gradient_uo( dgdx_old, thyra_u_old.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[0], false, dt);
    Thyra::Vp_V(dgdx.ptr(), *dgdx_old);
    Thyra::Vt_S(dgdx.ptr(), dt/2.);
  }
}

template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
gradient_un( ROL::Vector<Real> &grad, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new, 
              const ROL::Vector<Real> &p, const ROL::TimeStamp<Real> &timeStamp ) const
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  if(onlyFinalTime_ && (int) timeStamp.k != Nt_-1) {
    if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
      *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_un final time of the time stamp " << timeStamp.t[timeStamp.t.size()-1] << " is not the final time." << std::endl;
    Thyra::assign(Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector().ptr(), Teuchos::ScalarTraits<Real>::zero());
    return;
  }
  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_un final time of the time stamp " << timeStamp.t[timeStamp.t.size()-1] << " is the final time." << std::endl;

  TEUCHOS_TEST_FOR_EXCEPTION(!onlyFinalTime_ && !useTrapezoidalTimeIntegration_, std::logic_error,
    std::endl <<
    "Not implemented yet" << std::endl);

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  ROL::ThyraVector<Real>  & thyra_dgdx = dynamic_cast<ROL::ThyraVector<Real>&>(grad);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  Real dt = timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0];

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  RCP<Thyra::VectorBase<Real> > dgdx = thyra_dgdx.getVector();
  gradient_un( dgdx, thyra_u_new.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[timeStamp.t.size()-1], true, dt);

  if (!onlyFinalTime_ && useTrapezoidalTimeIntegration_) {
    RCP<Thyra::VectorBase<Real> > dgdx_old = thyra_dgdx.getVector()->clone_v();
    gradient_un( dgdx_old, thyra_u_old.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[0], false, dt);
    Thyra::Vp_V(dgdx.ptr(), *dgdx_old);
    Thyra::Vt_S(dgdx.ptr(), dt/2.);
  }
}

template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
gradient_z( ROL::Vector<Real> &grad, const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new, 
              const ROL::Vector<Real> &p, const ROL::TimeStamp<Real> &timeStamp ) const
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  if(onlyFinalTime_ && (int) timeStamp.k != Nt_-1) {
    if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
      *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_z final time of the time stamp " << timeStamp.t[timeStamp.t.size()-1] << " is not the final time." << std::endl;
    Thyra::assign(Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector().ptr(), Teuchos::ScalarTraits<Real>::zero());
    return;
  }
  if(verbosityLevel_ >= Teuchos::VERB_EXTREME)
    *out_ << "Piro::ThyraProductME_ROL_DynamicObjective::gradient_z final time of the time stamp " << timeStamp.t[timeStamp.t.size()-1] << " is the final time." << std::endl;

  TEUCHOS_TEST_FOR_EXCEPTION(!onlyFinalTime_ && !useTrapezoidalTimeIntegration_, std::logic_error,
    std::endl <<
    "Not implemented yet" << std::endl);

  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);

  ROL::ThyraVector<Real>  & thyra_dgdp = dynamic_cast<ROL::ThyraVector<Real>&>(grad);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  Real dt = timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0];

  RCP<Thyra::VectorBase<Real> > u_dot = compute_u_dot(thyra_u_new.getVector(), thyra_u_old.getVector(), dt);

  RCP<Thyra::VectorBase<Real> > dgdp = thyra_dgdp.getVector();
  gradient_z( dgdp, thyra_u_new.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[timeStamp.t.size()-1]);

  if (!onlyFinalTime_ && useTrapezoidalTimeIntegration_) {
    RCP<Thyra::VectorBase<Real> > dgdp_old = thyra_dgdp.getVector()->clone_v();
    gradient_z( dgdp_old, thyra_u_old.getVector(), u_dot, thyra_p.getVector(), timeStamp.t[0]);
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
gradient_uo( Teuchos::RCP<Thyra::VectorBase<Real>> &grad, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot, 
              const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t, const bool u_new, const Real &dt) const
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

  RCP<Thyra::VectorBase<Real> > dgdx = grad->clone_v();

  if (!u_new) {
    const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, g_index_);
    Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient;
    if (dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
      dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
    else if(dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
      dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Piro::ThyraProductME_ROL_DynamicObjective::gradient_un: DgDx does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
    }

    outArgs.set_DgDx(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(dgdx, dgdx_orient));
  }

  bool use_dgdx_dot = true;

  const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_dot_support =
      outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx_dot, g_index_);
  Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_dot_orient;
  if (dgdx_dot_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
    dgdx_dot_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
  else if(dgdx_dot_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
    dgdx_dot_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
  else {
    use_dgdx_dot = false;
  }

  if (u_new && !use_dgdx_dot)
    return;

  RCP<Thyra::VectorBase<Real> > dgdx_dot = grad->clone_v();

  if (use_dgdx_dot)
    outArgs.set_DgDx_dot(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(dgdx_dot, dgdx_dot_orient));

  thyra_model_->evalModel(inArgs, outArgs);

  if (use_dgdx_dot) {
    if (u_new) {
      Thyra::V_StV(grad.ptr(), 1./dt, *dgdx_dot);
    }
    else {
      Thyra::V_StVpStV(grad.ptr(), 1., *dgdx, -1./dt, *dgdx_dot);
    }
  }
}

template <typename Real>
void
ThyraProductME_ROL_DynamicObjective<Real>::
gradient_un( Teuchos::RCP<Thyra::VectorBase<Real>> &grad, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x, const Teuchos::RCP<const Thyra::VectorBase<Real>> &x_dot, 
              const Teuchos::RCP<const Thyra::VectorBase<Real>> &p, const Real &t, const bool u_new, const Real &dt) const
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

  RCP<Thyra::VectorBase<Real> > dgdx = grad->clone_v();

  if (u_new) {
    const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, g_index_);
    Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient;
    if (dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
      dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
    else if(dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
      dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Piro::ThyraProductME_ROL_DynamicObjective::gradient_un: DgDx does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
    }

    outArgs.set_DgDx(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(dgdx, dgdx_orient));
  }

  bool use_dgdx_dot = true;

  const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_dot_support =
      outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx_dot, g_index_);
  Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_dot_orient;
  if (dgdx_dot_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
    dgdx_dot_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
  else if(dgdx_dot_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
    dgdx_dot_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
  else {
    use_dgdx_dot = false;
  }

  if (!u_new && !use_dgdx_dot)
    return;

  RCP<Thyra::VectorBase<Real> > dgdx_dot = grad->clone_v();

  if (use_dgdx_dot)
    outArgs.set_DgDx_dot(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(dgdx_dot, dgdx_dot_orient));

  thyra_model_->evalModel(inArgs, outArgs);

  if (use_dgdx_dot) {
    if (u_new) {
      Thyra::V_StVpStV(grad.ptr(), 1., *dgdx, 1./dt, *dgdx_dot);
    }
    else {
      Thyra::V_StV(grad.ptr(), -1./dt, *dgdx_dot);
    }
  }
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