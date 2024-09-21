// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_THYRAPRODUCTME_TEMPUS_FINALOBJECTIVE_HPP
#define PIRO_THYRAPRODUCTME_TEMPUS_FINALOBJECTIVE_HPP

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

#include "Piro_TempusIntegrator.hpp"

namespace Piro {

template <typename Real>
class ThyraProductME_TempusFinalObjective : public virtual ROL::Objective<Real> {
public:

  ThyraProductME_TempusFinalObjective(
    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & model,
    const Teuchos::RCP<Piro::TempusSolver<Real>> & piroTempusSolver,
    const Teuchos::RCP<Tempus::Integrator<Real> >& integrator,
    const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & modelAdjoin,
    int g_index,
    Teuchos::ParameterList& piroParams,
    const int Nt,
    Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH,
    Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null);

  virtual ~ThyraProductME_TempusFinalObjective() {}

  //! Compute value of objective
  Real value( const ROL::Vector<Real> &p, Real &tol ) override;

  void gradient( ROL::Vector<Real> &grad, const ROL::Vector<Real> &p, Real &tol ) override;

  //! Helper function to run tempus, computing responses and derivatives
  void run_tempus(ROL::Vector<Real>& r, const ROL::Vector<Real>& p) const;
  void run_tempus(const Thyra::ModelEvaluatorBase::InArgs<Real>&  inArgs,
                  const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs) const;
  
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
  Teuchos::RCP<Piro::TempusIntegrator<Real>> piroTempusIntegrator_;
  std::string sensitivity_method_;
  const int g_index_;
  int Nt_;
  Real objectiveRecoveryValue_;
  bool useObjectiveRecoveryValue_;
  ROL::UpdateType updateType_;

  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;
  Teuchos::RCP<ROL_ObserverBase<Real>> observer_;

  Teuchos::RCP<Teuchos::ParameterList> tempus_params_;
  Real time_final_;

}; // class ThyraProductME_TempusFinalObjective

template <typename Real>
ThyraProductME_TempusFinalObjective<Real>::
ThyraProductME_TempusFinalObjective(
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & model,
  const Teuchos::RCP<Piro::TempusSolver<Real>> & piroTempusSolver,
  const Teuchos::RCP<Tempus::Integrator<Real> >& integrator,
  const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & modelAdjoin,
  int g_index,
  Teuchos::ParameterList& piroParams,
  const int Nt,
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> observer) :
  integrator_(integrator),
  thyra_model_(model),
  thyra_adjoint_model_(modelAdjoin),
  piroTempusIntegrator_(piroTempusSolver->getNonconstPiroTempusIntegrator()),
  sensitivity_method_("Adjoint"),
  g_index_(g_index),
  Nt_(Nt),
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel),
  observer_(observer),
  tempus_params_(Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList(piroParams.sublist("Tempus")))),
  time_final_(piroParams.get<Real>("Final Time", 1.))
{
  
}

template <typename Real>
Real
ThyraProductME_TempusFinalObjective<Real>::
value( const ROL::Vector<Real> &p, Real &tol )
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
    *out_ << "Piro::ThyraProductME_TempusFinalObjective::value" << std::endl;

  // Run tempus and compute response for specified parameter values
  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  inArgs.set_p(0, thyra_p.getVector());
  RCP<Thyra::VectorBase<Real> > g =
    Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));
  outArgs.set_g(g_index_, g);
  run_tempus(inArgs, outArgs);

  return ::Thyra::get_ele(*g,0);
}

template <typename Real>
void
ThyraProductME_TempusFinalObjective<Real>::
gradient( ROL::Vector<Real> &grad, const ROL::Vector<Real> &p, Real &tol )
{
  *out_ << "Piro::ThyraProductME_TempusFinalObjective::gradient" << std::endl;

  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  // Run tempus and compute response gradient for specified parameter values
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
      Teuchos::RCP<const Piro::ProductModelEvaluator<Real>> model_PME = 
        Teuchos::rcp_dynamic_cast<const Piro::ProductModelEvaluator<Real>>(thyra_model_);
      if (model_PME.is_null()) {
        Teuchos::RCP<const Thyra::ModelEvaluatorDelegatorBase<Real>> model_MEDB =
          Teuchos::rcp_dynamic_cast<const Thyra::ModelEvaluatorDelegatorBase<Real>>(thyra_model_);
        if (!model_MEDB.is_null()) {
          model_PME = Teuchos::rcp_dynamic_cast<const Piro::ProductModelEvaluator<Real>>(model_MEDB->getUnderlyingModel());
        }
      }

      if ( !model_PME.is_null()) {
        Teko::BlockedLinearOp dgdp_op =
            Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(model_PME->create_DgDp_op(g_index_, 0, prodvec_dgdp));
        Thyra::ModelEvaluatorBase::Derivative<Real> dgdp_der(Teuchos::rcp_dynamic_cast<Thyra::LinearOpBase<Real>>(dgdp_op));
        outArgs.set_DgDp(g_index_, 0, dgdp_der);
      }
      else {
        ROL_TEST_FOR_EXCEPTION( true, std::logic_error, "Piro::ThyraProductME_TempusFinalObjective::gradient_z: dgdp is not supported for the used ModelEvaluator.");
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
            "Piro::ThyraProductME_TempusFinalObjective::gradient_z: DgDp does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
      }
      outArgs.set_DgDp(g_index_, 0, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdp.getVector(), dgdp_orient));
    }
  }

  outArgs.set_g(g_index_, g);
  run_tempus(inArgs, outArgs);
}

template <typename Real>
void
ThyraProductME_TempusFinalObjective<Real>::
run_tempus(ROL::Vector<Real>& r, const ROL::Vector<Real>& p) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
  const ROL::ThyraVector<Real>& thyra_p =
    Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p =
    Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
  inArgs.set_p(0, thyra_p.getVector());
  ROL::ThyraVector<Real>& thyra_r =
    Teuchos::dyn_cast<ROL::ThyraVector<Real> >(r);
  outArgs.set_g(g_index_, thyra_r.getVector());
  run_tempus(inArgs, outArgs);
}

template <typename Real>
void
ThyraProductME_TempusFinalObjective<Real>::
run_tempus(const Thyra::ModelEvaluatorBase::InArgs<Real>&  inArgs,
           const Thyra::ModelEvaluatorBase::OutArgs<Real>& outArgs) const
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultNominalBoundsOverrideModelEvaluator<Real> DNBOME;

  // Override nominal values in model to supplied inArgs
  RCP<DNBOME> wrapped_model = rcp(new DNBOME(thyra_model_, rcpFromRef(inArgs)));
  RCP<DNBOME> wrapped_adjoint_model = rcp(new DNBOME(thyra_adjoint_model_, rcpFromRef(inArgs)));


  Real t;
  RCP<const Thyra::VectorBase<Real> > x, x_dot;
  RCP<const Thyra::MultiVectorBase<double> > dxdp, dxdotdp;
  RCP<Thyra::VectorBase<Real> > g = outArgs.get_g(g_index_);
  const bool compute_dgdp = !outArgs.get_DgDp(g_index_, 0).isEmpty();


  RCP<Thyra::MultiVectorBase<Real> > dgdp_mv =
    outArgs.get_DgDp(g_index_, 0).getMultiVector();
  MEB::EDerivativeMultiVectorOrientation dgdp_orientation =
    outArgs.get_DgDp(g_index_, 0).getMultiVectorOrientation();
  Teko::BlockedLinearOp dgdp_op =
    Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(outArgs.get_DgDp(g_index_, 0).getLinearOp());

  piroTempusIntegrator_->clearSolutionHistory();

  // Create and run integrator
  if (compute_dgdp && sensitivity_method_ == "Adjoint") {
    RCP<Tempus::IntegratorAdjointSensitivity<Real> > integrator =
      Tempus::createIntegratorAdjointSensitivity<Real>(tempus_params_, wrapped_model, wrapped_adjoint_model);
    const bool integratorStatus = integrator->advanceTime(time_final_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorStatus, std::logic_error, "Integrator failed!");

    // Get final state
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
    if(!dgdp_mv.is_null())
      Thyra::assign(dgdp_mv.ptr(), *(integrator->getDgDp()));
    else {
      //dgdp_op->setBlock(0, 0, piroTempusIntegrator->getDgDp());
      dgdp_mv = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Real>>( Teuchos::rcp_dynamic_cast<Thyra::DefaultScaledAdjointLinearOp<Real>>(dgdp_op->getNonconstBlock(0, 0))->getNonconstOp() );
      Thyra::assign(dgdp_mv.ptr(), *(integrator->getDgDp()));
    }
  }
  else {
    RCP<Tempus::IntegratorBasic<Real> > integrator =
      Tempus::createIntegratorBasic<Real>(tempus_params_, wrapped_model);
    const bool integratorStatus = integrator->advanceTime(time_final_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorStatus, std::logic_error, "Integrator failed!");

    // Get final state
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
  }

  // Evaluate response at final state
  MEB::InArgs<Real> modelInArgs   = inArgs;
  MEB::OutArgs<Real> modelOutArgs = outArgs;
  modelInArgs.set_x(x);
  if (modelInArgs.supports(MEB::IN_ARG_x_dot)) modelInArgs.set_x_dot(x_dot);
  if (modelInArgs.supports(MEB::IN_ARG_t)) modelInArgs.set_t(t);
  RCP<Thyra::MultiVectorBase<Real> > dgdx, dgdxdot;
  if (compute_dgdp && sensitivity_method_ == "Adjoint") {
    // Clear dg/dp as an out arg since it was already computed by the adjoint
    // integrator
    modelOutArgs.set_DgDp(g_index_, 0,
                          MEB::Derivative<Real>());
  }

  thyra_model_->evalModel(modelInArgs, modelOutArgs);

  // dg/dp = dg/dp + dg/dx*dx/dp + dg/dx_dot*dx_dot/dp
  // We assume dg/dx, dg/dxdot are in gradient form while dxdp, dxdotdp are in
  // Jacobian form
  if (compute_dgdp && dgdx != Teuchos::null) {
    if (dgdp_orientation == MEB::DERIV_MV_JACOBIAN_FORM)
      dgdx->apply(Thyra::TRANS, *dxdp, dgdp_mv.ptr(), Real(1.0), Real(1.0));
    else
      dxdp->apply(Thyra::TRANS, *dgdx, dgdp_mv.ptr(), Real(1.0), Real(1.0));
  }
  if (compute_dgdp && dgdxdot != Teuchos::null) {
    if (dgdp_orientation == MEB::DERIV_MV_JACOBIAN_FORM)
      dgdxdot->apply(Thyra::TRANS, *dxdotdp, dgdp_mv.ptr(), Real(1.0), Real(1.0));
    else
      dxdotdp->apply(Thyra::TRANS, *dgdxdot, dgdp_mv.ptr(), Real(1.0), Real(1.0));
  }
}


} // namespace Piro

#endif