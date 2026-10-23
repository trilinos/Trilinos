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
#include "Piro_TempusSolver.hpp"
#include "Piro_ROL_ObserverBase.hpp"

namespace Piro {

/** \brief Reduced objective g(x(T), p) for transient problems, computed with Tempus.
 *
 * The value is computed by integrating the forward problem with Tempus up to the
 * final time T of the Tempus time step control. The gradient is computed with the
 * Tempus adjoint sensitivity integrator.
 *
 * The parameters are the product parameter vector (parameter 0) of the model
 * evaluator, so several parameter vectors are supported through
 * Piro::ProductModelEvaluator. The response must be scalar.
 *
 * If no adjoint model evaluator is provided, Tempus uses an implicit adjoint
 * model, which requires a linear solver that supports transposed solves.
 *
 * When a forward time integration fails, the "State Solve Converged" flag of
 * the "Optimization Status" sublist is set to false. In this case value()
 * returns the response at the last computed state if an "Objective Recovery
 * Value" is set in "Optimization Status" (see ThyraProductME_ROL_ObservedObjective,
 * which replaces it with the recovery value), and throws otherwise.
 * gradient() always throws.
 */
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

  void update( const ROL::Vector<Real> &x, ROL::UpdateType type, int iter = -1 ) override {
    (void) x;
    (void) type;
    (void) iter;
  }

  /** \brief State at the final time computed by the last forward integration,
   *  or null if the last integration was not performed with parameters p. */
  Teuchos::RCP<const Thyra::VectorBase<Real>> getFinalState(const ROL::Vector<Real>& p) const;

private:
  //! Integrate the forward problem (and the adjoint problem if dgdp is not null)
  //! with parameters p, and compute the response g at the final time.
  //! Returns false if the time integration failed.
  bool run_tempus(const Teuchos::RCP<const Thyra::VectorBase<Real>>& p,
                  const Teuchos::RCP<Thyra::VectorBase<Real>>& g,
                  const Teuchos::RCP<Thyra::VectorBase<Real>>& dgdp);

  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_model_;
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_adjoint_model_;
  Teuchos::RCP<Piro::TempusIntegrator<Real>> piroTempusIntegrator_;
  const int g_index_;

  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;

  Teuchos::RCP<Teuchos::ParameterList> tempus_params_;
  Real time_final_;

  //! Parameters and final state of the last successful forward integration
  Teuchos::RCP<Thyra::VectorBase<Real>> last_p_;
  Teuchos::RCP<Thyra::VectorBase<Real>> last_x_;

}; // class ThyraProductME_TempusFinalObjective

template <typename Real>
ThyraProductME_TempusFinalObjective<Real>::
ThyraProductME_TempusFinalObjective(
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & model,
  const Teuchos::RCP<Piro::TempusSolver<Real>> & piroTempusSolver,
  const Teuchos::RCP<Tempus::Integrator<Real> >& /* integrator */,
  const Teuchos::RCP<Tempus::Integrator<Real>> & /* adjoint_integrator */,
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & modelAdjoin,
  int g_index,
  Teuchos::ParameterList& piroParams,
  const int /* Nt */,
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> /* observer */) :
  thyra_model_(model),
  thyra_adjoint_model_(modelAdjoin),
  piroTempusIntegrator_(piroTempusSolver->getNonconstPiroTempusIntegrator()),
  g_index_(g_index),
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel),
  tempus_params_(Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList(piroParams.sublist("Tempus"))))
{
  TEUCHOS_TEST_FOR_EXCEPTION(thyra_model_->get_g_space(g_index_)->dim() != 1, Teuchos::Exceptions::InvalidParameter,
    std::endl << "Piro::ThyraProductME_TempusFinalObjective: the objective response " << g_index_ <<
    " has dimension " << thyra_model_->get_g_space(g_index_)->dim() << ", but it must be a scalar." << std::endl);

  // Integrate to the final time of the Tempus time step control used by the Piro::TempusSolver.
  time_final_ = piroTempusSolver->getFinalTime();

  // The adjoint integrator differentiates the response and parameter selected in the
  // Tempus "Sensitivities" sublist. Make sure they are the objective response and the
  // control (parameter 0 of the product model evaluator), otherwise value() and
  // gradient() would refer to different functions.
  Teuchos::ParameterList& sensPL = tempus_params_->sublist("Sensitivities");
  if (sensPL.isParameter("Response Function Index") && sensPL.get<int>("Response Function Index") != g_index_) {
    *out_ << "WARNING: Piro::ThyraProductME_TempusFinalObjective: overriding Tempus->Sensitivities->"
          << "\"Response Function Index\" = " << sensPL.get<int>("Response Function Index")
          << " with the ROL \"Response Vector Index\" = " << g_index_ << std::endl;
  }
  sensPL.set<int>("Response Function Index", g_index_);
  if (sensPL.isParameter("Sensitivity Parameter Index") && sensPL.get<int>("Sensitivity Parameter Index") != 0) {
    *out_ << "WARNING: Piro::ThyraProductME_TempusFinalObjective: overriding Tempus->Sensitivities->"
          << "\"Sensitivity Parameter Index\" = " << sensPL.get<int>("Sensitivity Parameter Index")
          << " with 0, the index of the optimization parameters in the product model evaluator" << std::endl;
  }
  sensPL.set<int>("Sensitivity Parameter Index", 0);
}

template <typename Real>
Real
ThyraProductME_TempusFinalObjective<Real>::
value( const ROL::Vector<Real> &p, Real &/* tol */ )
{
  if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
    *out_ << "Piro::ThyraProductME_TempusFinalObjective::value" << std::endl;

  const ROL::ThyraVector<Real>& thyra_p = Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  Teuchos::RCP<Thyra::VectorBase<Real> > g = Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));

  const bool success = run_tempus(thyra_p.getVector(), g, Teuchos::null);

  TEUCHOS_TEST_FOR_EXCEPTION(!success && !optParams_.isParameter("Objective Recovery Value"), std::runtime_error,
    std::endl << "Piro::ThyraProductME_TempusFinalObjective::value: the time integration failed. "
    "Set the ROL parameter \"Objective Recovery Value\" to continue the optimization in this case." << std::endl);

  return ::Thyra::get_ele(*g,0);
}

template <typename Real>
void
ThyraProductME_TempusFinalObjective<Real>::
gradient( ROL::Vector<Real> &grad, const ROL::Vector<Real> &p, Real &/* tol */ )
{
  if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
    *out_ << "Piro::ThyraProductME_TempusFinalObjective::gradient" << std::endl;

  const ROL::ThyraVector<Real>& thyra_p = Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  ROL::ThyraVector<Real>& thyra_dgdp = Teuchos::dyn_cast<ROL::ThyraVector<Real> >(grad);
  Teuchos::RCP<Thyra::VectorBase<Real> > g = Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));

  bool success = false;
  std::string tempusError;
  try {
    success = run_tempus(thyra_p.getVector(), g, thyra_dgdp.getVector());
  }
  catch (const std::exception& e) {
    // After a failed forward integration the Tempus adjoint integrator can throw
    optParams_.set("State Solve Converged", false);
    tempusError = e.what();
  }

  TEUCHOS_TEST_FOR_EXCEPTION(!success, std::runtime_error,
    std::endl << "Piro::ThyraProductME_TempusFinalObjective::gradient: the time integration failed, "
    "the gradient cannot be computed." << (tempusError.empty() ? "" : " Tempus error:\n") << tempusError << std::endl);
}

template <typename Real>
Teuchos::RCP<const Thyra::VectorBase<Real>>
ThyraProductME_TempusFinalObjective<Real>::
getFinalState(const ROL::Vector<Real>& p) const
{
  if (last_p_.is_null() || last_x_.is_null())
    return Teuchos::null;
  const ROL::ThyraVector<Real>& thyra_p = Teuchos::dyn_cast<const ROL::ThyraVector<Real> >(p);
  Teuchos::RCP<Thyra::VectorBase<Real>> diff = last_p_->clone_v();
  Thyra::Vp_StV(diff.ptr(), Real(-1.0), *thyra_p.getVector());
  return (Thyra::norm_inf(*diff) == Real(0.0)) ? last_x_ : Teuchos::null;
}

template <typename Real>
bool
ThyraProductME_TempusFinalObjective<Real>::
run_tempus(const Teuchos::RCP<const Thyra::VectorBase<Real>>& p,
           const Teuchos::RCP<Thyra::VectorBase<Real>>& g,
           const Teuchos::RCP<Thyra::VectorBase<Real>>& dgdp)
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultNominalBoundsOverrideModelEvaluator<Real> DNBOME;

  // Override the nominal parameters of the model with p
  MEB::InArgs<Real> inArgs = thyra_model_->getNominalValues();
  inArgs.set_p(0, p);
  const RCP<const MEB::InArgs<Real>> nominalValues = rcp(new MEB::InArgs<Real>(inArgs));
  RCP<DNBOME> wrapped_model = rcp(new DNBOME(thyra_model_, nominalValues));
  // If the adjoint model is null, Tempus uses an implicit adjoint model.
  RCP<DNBOME> wrapped_adjoint_model = thyra_adjoint_model_.is_null() ? Teuchos::null :
    rcp(new DNBOME(thyra_adjoint_model_, nominalValues));

  piroTempusIntegrator_->clearSolutionHistory();

  // Create and run integrator
  bool integratorStatus;
  Real t;
  RCP<const Thyra::VectorBase<Real> > x, x_dot;
  if (Teuchos::nonnull(dgdp)) {
    RCP<Tempus::IntegratorAdjointSensitivity<Real> > integrator =
      Tempus::createIntegratorAdjointSensitivity<Real>(tempus_params_, wrapped_model, wrapped_adjoint_model);
    integratorStatus = integrator->advanceTime(time_final_);
    if (integratorStatus) {
      // dg/dp is computed in gradient form, with one column per response component
      Thyra::assign(dgdp.ptr(), *(integrator->getDgDp()->col(0)));
    }
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
  }
  else {
    RCP<Tempus::IntegratorBasic<Real> > integrator =
      Tempus::createIntegratorBasic<Real>(tempus_params_, wrapped_model);
    integratorStatus = integrator->advanceTime(time_final_);
    t = integrator->getTime();
    x = integrator->getX();
    x_dot = integrator->getXDot();
  }

  optParams_.set("State Solve Converged", integratorStatus);
  if (!integratorStatus && verbosityLevel_ >= Teuchos::VERB_LOW)
    *out_ << "Piro::ThyraProductME_TempusFinalObjective: the time integration failed at time " << t << std::endl;

  // Evaluate the response at the final state
  MEB::InArgs<Real> modelInArgs = inArgs;
  MEB::OutArgs<Real> modelOutArgs = thyra_model_->createOutArgs();
  modelInArgs.set_x(x);
  if (modelInArgs.supports(MEB::IN_ARG_x_dot)) modelInArgs.set_x_dot(x_dot);
  if (modelInArgs.supports(MEB::IN_ARG_t)) modelInArgs.set_t(t);
  modelOutArgs.set_g(g_index_, g);
  thyra_model_->evalModel(modelInArgs, modelOutArgs);

  // Store the final state for the observer
  if (integratorStatus) {
    last_p_ = p->clone_v();
    last_x_ = x->clone_v();
  }
  else {
    last_p_ = Teuchos::null;
    last_x_ = Teuchos::null;
  }

  return integratorStatus;
}

} // namespace Piro

#endif
