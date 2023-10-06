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

#include "Tempus_TimeDerivative.hpp"
#include "Tempus_StepperBackwardEuler_decl.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "ROL_Objective.hpp"
#include "ROL_DynamicObjective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_ThyraVector.hpp"

namespace Piro {

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
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> observer) :
  integrator_(integrator),
  thyra_model_(model),
  thyra_adjoint_model_(modelAdjoin),
  sensitivity_method_("Adjoint"),
  g_index_(g_index),
  Nt_(Nt),
  onlyFinalTime_(onlyFinalTime),
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

  // Run tempus and compute response for specified parameter values
  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  inArgs.set_p(0, thyra_p.getVector());
  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));
  outArgs.set_g(g_index_, g);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  RCP<Thyra::VectorBase<Real> > u_dot = thyra_u_new.getVector()->clone_v();
  Teuchos::RCP<Tempus::TimeDerivative<Real> > timeDer = integrator_->getStepper()->getTimeDerivative(timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0], thyra_u_old.getVector());
  //  Teuchos::rcp(new Tempus::StepperBackwardEulerTimeDerivative<Real>(Real(1.0)/(timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0]),thyra_u_old.getVector()));
  timeDer->compute(thyra_u_new.getVector(), u_dot);


  inArgs.set_x(thyra_u_new.getVector());
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(timeStamp.t[timeStamp.t.size()-1]);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(u_dot);

  thyra_model_->evalModel(inArgs, outArgs);

  return ::Thyra::get_ele(*g,0);
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

  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  inArgs.set_p(0, thyra_p.getVector());
  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));

  ROL::ThyraVector<Real>  & thyra_dgdx = dynamic_cast<ROL::ThyraVector<Real>&>(grad);

  //const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
  //    outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, g_index_);
  //Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient;
  //if (dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
  //  dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
  //else if(dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
  //  dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
  //else {
  //  ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
  //      "Piro::ThyraProductME_ROL_DynamicObjective::gradient_un: DgDx does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
  //}
  //
  //outArgs.set_DgDx(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdx.getVector(), dgdx_orient));
  
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

  if (use_dgdx_dot) {
    outArgs.set_DgDx_dot(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdx.getVector(), dgdx_dot_orient));
  }
  else {
    Thyra::assign(Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad).getVector().ptr(), Teuchos::ScalarTraits<Real>::zero());
    return;
  }

  outArgs.set_g(g_index_, g);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  //Real dt = timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0];

  RCP<Thyra::VectorBase<Real> > u_dot = thyra_u_new.getVector()->clone_v();
  Teuchos::RCP<Tempus::TimeDerivative<Real> > timeDer = integrator_->getStepper()->getTimeDerivative(timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0], thyra_u_old.getVector());
  //  Teuchos::rcp(new Tempus::StepperBackwardEulerTimeDerivative<Real>(Real(1.0)/dt,thyra_u_old.getVector()));
  timeDer->compute(thyra_u_new.getVector(), u_dot);

  inArgs.set_x(thyra_u_new.getVector());
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(timeStamp.t[timeStamp.t.size()-1]);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(u_dot);

  thyra_model_->evalModel(inArgs, outArgs);

  Thyra::V_S(thyra_dgdx.getVector().ptr(), timeDer->get_DxDot_Dx_old());
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

  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  inArgs.set_p(0, thyra_p.getVector());
  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));

  ROL::ThyraVector<Real>  & thyra_dgdx = dynamic_cast<ROL::ThyraVector<Real>&>(grad);

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

  outArgs.set_DgDx(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdx.getVector(), dgdx_orient));

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

  RCP<Thyra::VectorBase<Real> > dgdx_dot = thyra_dgdx.getVector()->clone_v();

  if (use_dgdx_dot)
    outArgs.set_DgDx_dot(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(dgdx_dot, dgdx_dot_orient));

  outArgs.set_g(g_index_, g);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  //Real dt = timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0];

  RCP<Thyra::VectorBase<Real> > u_dot = thyra_u_new.getVector()->clone_v();
  Teuchos::RCP<Tempus::TimeDerivative<Real> > timeDer = integrator_->getStepper()->getTimeDerivative(timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0], thyra_u_old.getVector());
  //  Teuchos::rcp(new Tempus::StepperBackwardEulerTimeDerivative<Real>(Real(1.0)/dt,thyra_u_old.getVector()));
  timeDer->compute(thyra_u_new.getVector(), u_dot);

  inArgs.set_x(thyra_u_new.getVector());
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(timeStamp.t[timeStamp.t.size()-1]);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(u_dot);

  thyra_model_->evalModel(inArgs, outArgs);

  if (use_dgdx_dot)
    Thyra::V_StV(thyra_dgdx.getVector().ptr(), timeDer->get_DxDot_Dx_new(), *dgdx_dot);
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

  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  inArgs.set_p(0, thyra_p.getVector());
  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));

  ROL::ThyraVector<Real>  & thyra_dgdp = dynamic_cast<ROL::ThyraVector<Real>&>(grad);
  Teuchos::RCP<Thyra::ProductMultiVectorBase<Real> > prodvec_dgdp =
      Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(thyra_dgdp.getVector());
  if ( !thyra_dgdp.getVector().is_null()) {
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
      outArgs.set_DgDp(g_index_, 0, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdp.getVector(), dgdp_orient));
    }
  }

  outArgs.set_g(g_index_, g);

  const ROL::ThyraVector<Real>& thyra_u_new =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_new);

  const ROL::ThyraVector<Real>& thyra_u_old =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(u_old);

  RCP<Thyra::VectorBase<Real> > u_dot = thyra_u_new.getVector()->clone_v();
  Teuchos::RCP<Tempus::TimeDerivative<Real> > timeDer = integrator_->getStepper()->getTimeDerivative(timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0], thyra_u_old.getVector());
  //  Teuchos::rcp(new Tempus::StepperBackwardEulerTimeDerivative<Real>(Real(1.0)/(timeStamp.t[timeStamp.t.size()-1] - timeStamp.t[0]),thyra_u_old.getVector()));
  timeDer->compute(thyra_u_new.getVector(), u_dot);

  inArgs.set_x(thyra_u_new.getVector());
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(timeStamp.t[timeStamp.t.size()-1]);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(u_dot);

  thyra_model_->evalModel(inArgs, outArgs);
}

} // namespace Piro

#endif