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

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorForwardSensitivity.hpp"
#include "Tempus_IntegratorAdjointSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientForwardSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientAdjointSensitivity.hpp"
#include "Tempus_IntegratorObserverComposite.hpp"
#include "Tempus_IntegratorObserverNoOp.hpp"
#include "Tempus_AdjointSensitivityUtils.hpp"

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

/** \brief Tempus observer that computes a "Time Integral" or "Sum Over Times"
 *  response during the forward time integration.
 *
 * - "Time Integral": J = int_{t0}^{T} g(x(t),p,t) dt, computed with the
 *   trapezoidal rule over the accepted time steps.
 * - "Sum Over Times": J = sum_k g(x(t_k),p,t_k), where the time steps must end
 *   at the times t_k (up to the tolerance). The model is evaluated exactly at
 *   the times t_k.
 */
template <typename Real>
class TempusResponseAccumulator : public Tempus::IntegratorObserverNoOp<Real> {
public:
  TempusResponseAccumulator(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Real>>& model, const int g_index,
    const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs,
    const Tempus::AdjointResponseType responseType,
    const std::vector<Real>& responseTimes, const Real timeTol) :
    model_(model), g_index_(g_index), inArgs_(inArgs), responseType_(responseType),
    responseTimes_(responseTimes), timeTol_(timeTol),
    value_(0), t_prev_(0), g_prev_(0), next_(0) {}

  void observeStartIntegrator(const Tempus::Integrator<Real>& integrator) override {
    const Teuchos::RCP<const Tempus::SolutionState<Real>> state =
      integrator.getSolutionHistory()->getCurrentState();
    value_ = 0;
    next_ = 0;
    missed_.clear();
    t_prev_ = state->getTime();
    if (responseType_ == Tempus::AdjointResponseType::TimeIntegral)
      g_prev_ = evalResponse(*state, t_prev_);
    else
      addPointResponses(*state);
  }

  void observeEndTimeStep(const Tempus::Integrator<Real>& integrator) override {
    const Teuchos::RCP<const Tempus::SolutionState<Real>> state =
      integrator.getSolutionHistory()->getCurrentState();
    const Real t = state->getTime();
    if (t <= t_prev_) return;  // the time step failed
    if (responseType_ == Tempus::AdjointResponseType::TimeIntegral) {
      const Real g = evalResponse(*state, t);
      value_ += (t - t_prev_) / Real(2) * (g_prev_ + g);
      g_prev_ = g;
    }
    t_prev_ = t;
    if (responseType_ == Tempus::AdjointResponseType::SumOverTimes)
      addPointResponses(*state);
  }

  //! Value of the response computed by the last time integration
  Real getValue() const { return value_; }

  //! Response times that are not the time of a time step, including those that
  //! have not been reached
  std::vector<Real> getMissedTimes() const {
    std::vector<Real> missed = missed_;
    missed.insert(missed.end(), responseTimes_.begin() + next_, responseTimes_.end());
    return missed;
  }

private:
  //! Add the responses at the times t_k up to the time of the state
  void addPointResponses(const Tempus::SolutionState<Real>& state) {
    const Real t = state.getTime();
    for (; next_ < responseTimes_.size() && responseTimes_[next_] <= t + timeTol_; ++next_) {
      if (std::abs(responseTimes_[next_] - t) <= timeTol_)
        value_ += evalResponse(state, responseTimes_[next_]);
      else
        missed_.push_back(responseTimes_[next_]);
    }
  }

  Real evalResponse(const Tempus::SolutionState<Real>& state, const Real t) const {
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Real> inArgs = inArgs_;
    inArgs.set_x(state.getX());
    if (inArgs.supports(MEB::IN_ARG_x_dot) && Teuchos::nonnull(state.getXDot()))
      inArgs.set_x_dot(state.getXDot());
    if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(t);
    MEB::OutArgs<Real> outArgs = model_->createOutArgs();
    const Teuchos::RCP<Thyra::VectorBase<Real>> g = Thyra::createMember<Real>(model_->get_g_space(g_index_));
    outArgs.set_g(g_index_, g);
    model_->evalModel(inArgs, outArgs);
    return Thyra::get_ele(*g, 0);
  }

  const Teuchos::RCP<const Thyra::ModelEvaluator<Real>> model_;
  const int g_index_;
  const Thyra::ModelEvaluatorBase::InArgs<Real> inArgs_;
  const Tempus::AdjointResponseType responseType_;
  const std::vector<Real> responseTimes_;
  const Real timeTol_;
  Real value_, t_prev_, g_prev_;
  std::size_t next_;
  std::vector<Real> missed_;
};

/** \brief Reduced objective for transient problems, computed with Tempus.
 *
 * Three types of objectives are supported (see setResponseType()), with g the
 * scalar response of the model:
 * - "Final Time" (default): J(p) = g(x(T),p,T)
 * - "Time Integral": J(p) = int_{t0}^{T} g(x(t),p,t) dt
 * - "Sum Over Times": J(p) = sum_k g(x(t_k),p,t_k), e.g., for observations
 *   available at the times t_k.
 *
 * The value is computed by integrating the forward problem with Tempus up to the
 * final time T of the Tempus time step control (the time integral with the
 * trapezoidal rule over the time steps). The gradient is computed with the
 * Tempus adjoint sensitivity integrator, with the "Response Type" and
 * "Response Time List" of the Tempus "Sensitivities" sublist set accordingly.
 * For "Time Integral" and "Sum Over Times" the response must not depend on
 * x_dot. For "Sum Over Times" the time steps land on the times t_k with
 * "Variable" time steps, with "Constant" time steps the times t_k must be on
 * the time grid.
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

  /** \brief Set the type of objective: "Final Time" (default), "Time Integral" or
   *  "Sum Over Times", with responseTimes the times t_k of "Sum Over Times". */
  void setResponseType(const Tempus::AdjointResponseType responseType,
                       const std::vector<Real>& responseTimes = std::vector<Real>());

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
  //! with parameters p, and compute the response g.
  //! Returns false if the time integration failed.
  bool run_tempus(const Teuchos::RCP<const Thyra::VectorBase<Real>>& p,
                  const Teuchos::RCP<Thyra::VectorBase<Real>>& g,
                  const Teuchos::RCP<Thyra::VectorBase<Real>>& dgdp);

  //! Tolerance for comparing times (same as in Tempus::IntegratorAdjointSensitivity)
  Real timeTolerance() const {
    return Real(1.0e-9) * std::max(std::abs(time_final_ - time_initial_), std::abs(time_final_) + Real(1.0e-6));
  }

  //! Warn (once) if a final-time response depends on x_dot: the Tempus adjoint
  //! sensitivities neglect dg/dx_dot.
  void checkFinalTimeDgDxDot(const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs);

  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_model_;
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_adjoint_model_;
  Teuchos::RCP<Piro::TempusIntegrator<Real>> piroTempusIntegrator_;
  const int g_index_;

  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;

  Teuchos::RCP<Teuchos::ParameterList> tempus_params_;
  Real time_initial_;
  Real time_final_;

  Tempus::AdjointResponseType responseType_;
  std::vector<Real> responseTimes_;
  bool checkedDgDxDot_;

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
  tempus_params_(Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList(piroParams.sublist("Tempus")))),
  responseType_(Tempus::AdjointResponseType::FinalTime),
  checkedDgDxDot_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(thyra_model_->get_g_space(g_index_)->dim() != 1, Teuchos::Exceptions::InvalidParameter,
    std::endl << "Piro::ThyraProductME_TempusFinalObjective: the objective response " << g_index_ <<
    " has dimension " << thyra_model_->get_g_space(g_index_)->dim() << ", but it must be a scalar." << std::endl);

  // Integrate to the final time of the Tempus time step control used by the Piro::TempusSolver.
  time_initial_ = piroTempusSolver->getStartTime();
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

  // The adjoint integration evaluates the forward solution between the stored
  // time steps (e.g., at the stage times of Runge-Kutta methods) by interpolation.
  // A stepper of order q needs an interpolation of order q-1 to retain its order.
  const Teuchos::RCP<Tempus::Stepper<Real>> stepper = piroTempusIntegrator_->getStepper();
  const std::string integratorName = tempus_params_->get<std::string>("Integrator Name", "");
  int interpolationOrder = 0;
  if (!integratorName.empty() && tempus_params_->isSublist(integratorName) &&
      tempus_params_->sublist(integratorName).isSublist("Solution History") &&
      tempus_params_->sublist(integratorName).sublist("Solution History").isSublist("Interpolator"))
    interpolationOrder = tempus_params_->sublist(integratorName).sublist("Solution History")
      .sublist("Interpolator").get<int>("Order", 0);
  if (Teuchos::nonnull(stepper) && stepper->getOrder() - 1 > interpolationOrder) {
    *out_ << "WARNING: Piro::ThyraProductME_TempusFinalObjective: the Tempus stepper \""
          << stepper->getStepperType() << "\" has order " << stepper->getOrder()
          << ", but the interpolation of the forward solution used by the adjoint integration has order "
          << interpolationOrder << ". Set Tempus->" << integratorName << "->Solution History->Interpolator->"
          << "\"Order\" to at least " << stepper->getOrder() - 1 << " to retain the order of the stepper "
          << "in the gradient." << std::endl;
  }
}

template <typename Real>
void
ThyraProductME_TempusFinalObjective<Real>::
setResponseType(const Tempus::AdjointResponseType responseType, const std::vector<Real>& responseTimes)
{
  responseType_ = responseType;
  responseTimes_.clear();
  std::string typeName = "Final Time";
  if (responseType_ == Tempus::AdjointResponseType::TimeIntegral)
    typeName = "Time Integral";
  else if (responseType_ == Tempus::AdjointResponseType::SumOverTimes) {
    typeName = "Sum Over Times";
    responseTimes_ = responseTimes;
    std::sort(responseTimes_.begin(), responseTimes_.end());
    TEUCHOS_TEST_FOR_EXCEPTION(responseTimes_.empty(), Teuchos::Exceptions::InvalidParameter,
      std::endl << "Piro::ThyraProductME_TempusFinalObjective::setResponseType: the \"Sum Over Times\" "
      "response requires at least one response time." << std::endl);
    const Real tol = timeTolerance();
    TEUCHOS_TEST_FOR_EXCEPTION(responseTimes_.front() < time_initial_ - tol || responseTimes_.back() > time_final_ + tol,
      Teuchos::Exceptions::InvalidParameter,
      std::endl << "Piro::ThyraProductME_TempusFinalObjective::setResponseType: the response times [" <<
      responseTimes_.front() << ", " << responseTimes_.back() << "] must be in the time interval [" <<
      time_initial_ << ", " << time_final_ << "]." << std::endl);
  }

  Teuchos::ParameterList& sensPL = tempus_params_->sublist("Sensitivities");
  sensPL.set<std::string>("Response Type", typeName);
  std::ostringstream timeList;
  timeList.precision(17);
  for (std::size_t k = 0; k < responseTimes_.size(); ++k)
    timeList << (k > 0 ? ", " : "") << responseTimes_[k];
  sensPL.set<std::string>("Response Time List", timeList.str());
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
    std::endl << "Piro::ThyraProductME_TempusFinalObjective::gradient: the time integration or the adjoint "
    "sensitivity computation failed, the gradient cannot be computed." << (tempusError.empty() ? "" : " Tempus error:\n") <<
    tempusError << std::endl);
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

  // For "Time Integral" and "Sum Over Times" responses the value is computed by an
  // observer during the forward integration
  RCP<TempusResponseAccumulator<Real>> accumulator;
  if (responseType_ != Tempus::AdjointResponseType::FinalTime)
    accumulator = rcp(new TempusResponseAccumulator<Real>(thyra_model_, g_index_, inArgs,
      responseType_, responseTimes_, timeTolerance()));
  auto addAccumulator = [&accumulator](auto& integrator) {
    if (accumulator.is_null()) return;
    RCP<Tempus::IntegratorObserverComposite<Real>> observers = rcp(new Tempus::IntegratorObserverComposite<Real>);
    if (Teuchos::nonnull(integrator.getObserver()))
      observers->addObserver(integrator.getObserver());
    observers->addObserver(accumulator);
    integrator.setObserver(observers);
  };

  // Create and run integrator
  bool integratorStatus;
  Real t;
  RCP<const Thyra::VectorBase<Real> > x, x_dot;
  if (Teuchos::nonnull(dgdp)) {
    RCP<Tempus::IntegratorAdjointSensitivity<Real> > integrator =
      Tempus::createIntegratorAdjointSensitivity<Real>(tempus_params_, wrapped_model, wrapped_adjoint_model);
    addAccumulator(*integrator);
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
    // Same time steps as in the forward integration of the adjoint sensitivity integrator
    if (responseType_ == Tempus::AdjointResponseType::SumOverTimes)
      Tempus::AdjointSensitivityUtils::addResponseTimeEvents(*integrator->getNonConstTimeStepControl(), responseTimes_);
    addAccumulator(*integrator);
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
  modelInArgs.set_x(x);
  if (modelInArgs.supports(MEB::IN_ARG_x_dot)) modelInArgs.set_x_dot(x_dot);
  if (modelInArgs.supports(MEB::IN_ARG_t)) modelInArgs.set_t(t);
  if (accumulator.is_null()) {
    MEB::OutArgs<Real> modelOutArgs = thyra_model_->createOutArgs();
    modelOutArgs.set_g(g_index_, g);
    thyra_model_->evalModel(modelInArgs, modelOutArgs);
    if (integratorStatus && Teuchos::nonnull(dgdp))
      checkFinalTimeDgDxDot(modelInArgs);
  }
  else {
    if (integratorStatus) {
      const std::vector<Real> missed = accumulator->getMissedTimes();
      std::ostringstream missedList;
      for (const Real t_k : missed) missedList << " " << t_k;
      TEUCHOS_TEST_FOR_EXCEPTION(!missed.empty(), std::logic_error,
        std::endl << "Piro::ThyraProductME_TempusFinalObjective: no time step ends at the response times" <<
        missedList.str() << ". With 'Variable' time steps, the time steps land on the response times; with " <<
        "'Constant' time steps, the response times must be on the time grid." << std::endl);
    }
    Thyra::set_ele(0, accumulator->getValue(), g.ptr());
  }

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

template <typename Real>
void
ThyraProductME_TempusFinalObjective<Real>::
checkFinalTimeDgDxDot(const Thyra::ModelEvaluatorBase::InArgs<Real>& inArgs)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  if (checkedDgDxDot_) return;
  checkedDgDxDot_ = true;
  MEB::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
  if (!outArgs.supports(MEB::OUT_ARG_DgDx_dot, g_index_).supports(MEB::DERIV_MV_GRADIENT_FORM))
    return;
  const Teuchos::RCP<Thyra::MultiVectorBase<Real>> dgdx_dot = Thyra::createMembers(thyra_model_->get_x_space(), 1);
  outArgs.set_DgDx_dot(g_index_, MEB::Derivative<Real>(dgdx_dot, MEB::DERIV_MV_GRADIENT_FORM));
  thyra_model_->evalModel(inArgs, outArgs);
  const Real norm = Thyra::norm_inf(*dgdx_dot->col(0));
  if (norm != Real(0)) {
    *out_ << "WARNING: Piro::ThyraProductME_TempusFinalObjective: the response " << g_index_
          << " depends on x_dot (|dg/dx_dot| = " << norm << " at the final time), but the Tempus adjoint "
          << "sensitivities neglect this dependence, so the gradient is not consistent with the objective. "
          << "Use \"Tempus Driver\" = false, or a response that does not depend on x_dot." << std::endl;
  }
}

} // namespace Piro

#endif
