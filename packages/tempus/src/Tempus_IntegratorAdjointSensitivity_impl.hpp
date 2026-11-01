//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorAdjointSensitivity_impl_hpp
#define Tempus_IntegratorAdjointSensitivity_impl_hpp

#include <algorithm>
#include <cmath>

#include "Teuchos_ParameterList.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_ImplicitAdjointModelEvaluator.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

#include "Tempus_StepperImplicit_decl.hpp"

namespace Tempus {

template <class Scalar>
IntegratorAdjointSensitivity<Scalar>::IntegratorAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<IntegratorBasic<Scalar>>& state_integrator,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_model,
    const Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar>>&
        adjoint_aux_model,
    const Teuchos::RCP<IntegratorBasic<Scalar>>& adjoint_integrator,
    const Teuchos::RCP<SolutionHistory<Scalar>>& solutionHistory,
    const int p_index, const int g_index, const bool g_depends_on_p,
    const bool f_depends_on_p, const bool ic_depends_on_p,
    const bool mass_matrix_is_identity, const AdjointResponseType response_type,
    const std::vector<Scalar>& response_times,
    const Teuchos::RCP<const DAEStructure<Scalar>>& dae_structure,
    const Teuchos::RCP<const Teuchos::ParameterList>& dae_jump_solver_pl)
  : model_(model),
    state_integrator_(state_integrator),
    adjoint_model_(adjoint_model),
    adjoint_aux_model_(adjoint_aux_model),
    adjoint_integrator_(adjoint_integrator),
    solutionHistory_(solutionHistory),
    p_index_(p_index),
    g_index_(g_index),
    g_depends_on_p_(g_depends_on_p),
    f_depends_on_p_(f_depends_on_p),
    ic_depends_on_p_(ic_depends_on_p),
    mass_matrix_is_identity_(mass_matrix_is_identity),
    response_type_(response_type),
    response_times_(response_times),
    stepMode_(SensitivityStepMode::Forward),
    dae_structure_(dae_structure)
{
  std::sort(response_times_.begin(), response_times_.end());
  TEUCHOS_TEST_FOR_EXCEPTION(
      response_type_ == AdjointResponseType::SumOverTimes &&
          response_times_.empty(),
      std::logic_error,
      "Error - IntegratorAdjointSensitivity(): the \"Sum Over Times\" "
      "response requires a non-empty \"Response Time List\".\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
      dae_structure_ != Teuchos::null && mass_matrix_is_identity_,
      std::logic_error,
      "Error - IntegratorAdjointSensitivity(): \""
          << DAEStructure<Scalar>::algebraicMaskName()
          << "\" cannot be used with \"Mass Matrix Is Identity\" = true.\n");
  if (dae_structure_ != Teuchos::null)
    dae_jump_solver_ = Teuchos::rcp(new AdjointDAEJumpSolver<Scalar>(dae_structure_, dae_jump_solver_pl));
  TEUCHOS_TEST_FOR_EXCEPTION(
      getStepper()->getUseFSAL(), std::logic_error,
      "Error - IntegratorAdjointSensitivity(): Cannot use FSAL with\n"
      "        IntegratorAdjointSensitivity, because the state and adjoint\n"
      "        integrators require ModelEvaluator evaluation in the\n"
      "        constructor to make the initial conditions consistent.\n"
      "        For the adjoint integrator, this requires special construction\n"
      "        which has not been implemented yet.\n");
}

template <class Scalar>
IntegratorAdjointSensitivity<Scalar>::IntegratorAdjointSensitivity()
{
  state_integrator_   = createIntegratorBasic<Scalar>();
  adjoint_integrator_ = createIntegratorBasic<Scalar>();
  response_type_      = AdjointResponseType::FinalTime;
  stepMode_           = SensitivityStepMode::Forward;
}

template <class Scalar>
bool IntegratorAdjointSensitivity<Scalar>::advanceTime()
{
  const Scalar tfinal = state_integrator_->getTimeStepControl()->getFinalTime();
  return advanceTime(tfinal);
}

template <class Scalar>
bool IntegratorAdjointSensitivity<Scalar>::advanceTime(const Scalar timeFinal)
{
  TEMPUS_FUNC_TIME_MONITOR_DIFF(
      "Tempus::IntegratorAdjointSensitivity::advanceTime()", TEMPUS_AS_AT);

  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::createMember;
  using Thyra::createMembers;
  using Thyra::LinearOpBase;
  using Thyra::LinearOpWithSolveBase;
  using Thyra::LinearOpWithSolveFactoryBase;
  using Thyra::MultiVectorBase;
  using Thyra::PreconditionerBase;
  using Thyra::PreconditionerFactoryBase;
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  typedef Thyra::DefaultProductVector<Scalar> DPV;

  // Get initial state for later
  RCP<const SolutionHistory<Scalar>> state_solution_history =
      state_integrator_->getSolutionHistory();
  RCP<const SolutionState<Scalar>> initial_state = (*state_solution_history)[0];

  // Run state integrator and get solution
  bool state_status = true;
  {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::IntegratorAdjointSensitivity::advanceTime::state",
        TEMPUS_AS_AT_FWD);
    stepMode_    = SensitivityStepMode::Forward;
    state_status = state_integrator_->advanceTime(timeFinal);
  }

  // The adjoint problem is not solved if the forward integration failed
  if (!state_status) {
    dgdp_ = Teuchos::null;
    return false;
  }

  // For at least some time-stepping methods, the time of the last time step
  // may not be timeFinal (e.g., it may be greater by at most delta_t).
  // But since the adjoint model requires timeFinal in its formulation, reset
  // it to the achieved final time.
  adjoint_aux_model_->setFinalTime(state_integrator_->getTime());

  // Set solution history in adjoint stepper
  adjoint_aux_model_->setForwardSolutionHistory(state_solution_history);

  RCP<const VectorSpaceBase<Scalar>> g_space = model_->get_g_space(g_index_);
  RCP<const VectorSpaceBase<Scalar>> x_space = model_->get_x_space();
  RCP<const VectorSpaceBase<Scalar>> p_space = model_->get_p_space(p_index_);
  const int num_g                            = g_space->dim();
  MEB::InArgs<Scalar> inargs                 = model_->getNominalValues();
  MEB::OutArgs<Scalar> adj_outargs           = adjoint_model_->createOutArgs();

  // The adjoint is integrated in the reversed time tau = t_final + t_init - t
  // (see AdjointAuxSensitivityModelEvaluator).
  const Scalar t_init   = state_integrator_->getTimeStepControl()->getInitTime();
  const Scalar t_final  = state_integrator_->getTime();
  const Scalar time_tol = AdjointSensitivityUtils::timeTolerance(
      t_init, t_final, time_tolerance_pl_.get());

  // Times at which the response contributes to the adjoint with a jump
  //   y += df/dx_dot^{-T} * dg/dx^T,
  // in decreasing order, i.e., in the order of the adjoint integration.
  // For a "Time Integral" response there are no jumps, the response is a
  // source term of the adjoint equations.
  std::vector<Scalar> jump_times;
  if (response_type_ == AdjointResponseType::FinalTime)
    jump_times.push_back(t_final);
  else if (response_type_ == AdjointResponseType::SumOverTimes) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        response_times_.front() < t_init - time_tol || response_times_.back() > t_final + time_tol,
        std::logic_error,
        "Error - IntegratorAdjointSensitivity::advanceTime(): the response "
        "times ["
            << response_times_.front() << ", " << response_times_.back()
            << "] must be in the time interval [" << t_init << ", " << t_final
            << "].\n");
    jump_times.assign(response_times_.rbegin(), response_times_.rend());
  }

  // For an index-1 DAE, check that the rows of the algebraic equations of the
  // mass matrix are zero (see the class description)
  if (dae_structure_ != Teuchos::null) {
    RCP<const SolutionState<Scalar>> final_state = state_solution_history->getCurrentState();
    MEB::InArgs<Scalar> check_inargs             = model_->getNominalValues();
    check_inargs.set_t(final_state->getTime());
    check_inargs.set_x(final_state->getX());
    if (check_inargs.supports(MEB::IN_ARG_x_dot))
      check_inargs.set_x_dot(final_state->getXDot());
    dae_structure_->checkMassMatrix(*model_, check_inargs, "Tempus::IntegratorAdjointSensitivity::advanceTime");
  }

  // The adjoint state (y, z) is zero at t_final, before the jumps.
  RCP<DPV> adjoint_state = rcp_dynamic_cast<DPV>(Thyra::createMember(adjoint_aux_model_->get_x_space()));
  assign(Teuchos::ptr_implicit_cast<VectorBase<Scalar>>(adjoint_state.ptr()), Teuchos::ScalarTraits<Scalar>::zero());
  RCP<MultiVectorBase<Scalar>> adjoint_y_mv =
      rcp_dynamic_cast<DMVPV>(adjoint_state->getNonconstVectorBlock(0))->getNonconstMultiVector();

  // Sum of the dg/dp terms of the jump times, sum_k dg/dp(t_k)^T
  dgdp_ = createMembers(p_space, num_g);
  assign(dgdp_.ptr(), Scalar(0.0));

  RCP<MultiVectorBase<Scalar>> dgdx = createMembers(x_space, num_g);
  RCP<MultiVectorBase<Scalar>> dgdp_k;
  RCP<MultiVectorBase<Scalar>> jump;
  RCP<MultiVectorBase<Scalar>> dae_jump;
  RCP<MultiVectorBase<Scalar>> algebraic_jump;
  RCP<LinearOpWithSolveBase<Scalar>> W;

  // Run the adjoint integration, restarted at each jump time
  const bool restart_at_jumps =
      response_type_ == AdjointResponseType::SumOverTimes;
  adjoint_solution_history_ = Teuchos::null;
  bool sens_status          = true;
  {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::IntegratorAdjointSensitivity::advanceTime::adjoint",
        TEMPUS_AS_AT_ADJ);
    stepMode_ = SensitivityStepMode::Adjoint;
    const Scalar tau_init =
        adjoint_integrator_->getTimeStepControl()->getInitTime();
    Scalar tau        = tau_init;
    int num_segments  = 0;
    int adjoint_index = -1;  // time step index of the adjoint integration
    for (std::size_t k = 0; sens_status && k < jump_times.size(); ++k) {
      const Scalar t_k   = jump_times[k];
      const Scalar tau_k = t_final + t_init - t_k;
      if (tau_k > tau + time_tol) {
        if (restart_at_jumps && num_segments > 0)
          appendAdjointSolutionHistory(time_tol);
        sens_status =
            integrateAdjoint(tau, tau_k, *adjoint_state, adjoint_index);
        tau = tau_k;
        ++num_segments;
        if (!sens_status) break;
      }

      // Forward state at t_k
      RCP<const SolutionState<Scalar>> state;
      if (response_type_ == AdjointResponseType::FinalTime)
        state = state_solution_history->getCurrentState();
      else
        state = findForwardState(*state_solution_history, t_k, time_tol);
      inargs.set_t(response_type_ == AdjointResponseType::FinalTime
                       ? state->getTime()
                       : t_k);
      inargs.set_x(state->getX());
      if (inargs.supports(MEB::IN_ARG_x_dot)) inargs.set_x_dot(state->getXDot());

      // Compute dg/dx and the jump of the adjoint.
      // For explicit form, we are relying on the user to inform us the
      // the mass matrix is the identity.  It would be nice to be able to
      // determine somehow automatically that we are using an explicit stepper.
      AdjointSensitivityUtils::evalResponseGradientX(
          *model_, g_index_, inargs, *dgdx,
          response_type_ != AdjointResponseType::FinalTime,
          "Tempus::IntegratorAdjointSensitivity::advanceTime");
      if (mass_matrix_is_identity_)
        Thyra::update(Scalar(1.0), *dgdx, adjoint_y_mv.ptr());
      else if (dae_structure_ == Teuchos::null) {
        // ODE: jump df/dx_dot^{-T} * dg/dx^T.  The mass matrix is constant, so
        // its solver is reused for all jumps.
        if (W == Teuchos::null) W = createAdjointMassMatrixSolver(inargs);
        if (jump == Teuchos::null) jump = createMembers(x_space, num_g);
        // Initialize jump to zero before solve for linear solvers that
        // use what is passed in as the initial guess
        assign(jump.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
        W->solve(Thyra::NOTRANS, *dgdx, jump.ptr());
        Thyra::update(Scalar(1.0), *jump, adjoint_y_mv.ptr());
      }
      else {
        // Index-1 DAE: C^T z = dg/dx^T, with C = P_a df/dx + P_d df/dx_dot at
        // the state at t_k (see the class description).
        if (dae_jump == Teuchos::null)
          dae_jump = createMembers(model_->get_f_space(), num_g);
        if (!dae_jump_solver_->solve(*adjoint_model_, inargs, state->getTimeStep(), *dgdx, *dae_jump)) {
          sens_status = false;
          break;
        }
        // The algebraic components of z add -df/dp^T (P_a z) to the sensitivity
        if (f_depends_on_p_) {
          if (algebraic_jump == Teuchos::null)
            algebraic_jump = createMembers(model_->get_f_space(), num_g);
          assign(algebraic_jump.ptr(), *dae_jump);
          dae_structure_->zeroDifferentialComponents(*algebraic_jump);
          AdjointSensitivityUtils::applyDfDpTranspose(
              *model_, p_index_, inargs, *algebraic_jump, *dgdp_, Scalar(-1.0),
              Scalar(1.0), dfdp_op_, dfdp_mv_);
        }
        // The differential components of z are the jump of the adjoint (the
        // algebraic adjoint unknowns are determined by the adjoint equations)
        dae_structure_->zeroAlgebraicComponents(*dae_jump);
        Thyra::update(Scalar(1.0), *dae_jump, adjoint_y_mv.ptr());
      }

      // Add dg/dp(t_k)^T
      if (g_depends_on_p_) {
        if (dgdp_k == Teuchos::null) dgdp_k = createMembers(p_space, num_g);
        AdjointSensitivityUtils::evalResponseGradientP(*model_, g_index_, p_index_, inargs, *dgdp_k);
        Thyra::update(Scalar(1.0), *dgdp_k, dgdp_.ptr());
      }

      // Store the adjoint state after the jump in the adjoint solution history
      if (num_segments > 0)
        assign(adjoint_integrator_->getNonConstSolutionHistory()->getCurrentState()->getX().ptr(),
               *adjoint_state);
    }

    // Integrate to the initial time, unless the last jump is at the initial
    // time (or the integration failed).
    if (sens_status && (num_segments == 0 || tau < timeFinal - time_tol)) {
      if (restart_at_jumps && num_segments > 0)
        appendAdjointSolutionHistory(time_tol);
      sens_status =
          integrateAdjoint(tau, timeFinal, *adjoint_state, adjoint_index);
      ++num_segments;
    }
    if (restart_at_jumps && num_segments > 0)
      appendAdjointSolutionHistory(time_tol);
  }
  RCP<const SolutionHistory<Scalar>> adjoint_solution_history = getSensSolutionHistory();

  // Add in initial condition term = (dx/dp^T(0))*(df/dx_dot^T(0))*y(0)
  // If dxdp_init_ is null, assume it is zero
  if (ic_depends_on_p_ && dxdp_init_ != Teuchos::null) {
    RCP<const MultiVectorBase<Scalar>> adjoint_mv = adjoint_y_mv;
    if (mass_matrix_is_identity_)
      dxdp_init_->apply(Thyra::CONJTRANS, *adjoint_mv, dgdp_.ptr(), Scalar(1.0),
                        Scalar(1.0));
    else {
      inargs.set_t(initial_state->getTime());
      inargs.set_x(initial_state->getX());
      if (inargs.supports(MEB::IN_ARG_x_dot))
        inargs.set_x_dot(initial_state->getXDot());
      inargs.set_alpha(1.0);
      inargs.set_beta(0.0);
      RCP<LinearOpBase<Scalar>> W_op = adjoint_model_->create_W_op();
      adj_outargs.set_W_op(W_op);
      adjoint_model_->evalModel(inargs, adj_outargs);
      adj_outargs.set_W_op(Teuchos::null);
      RCP<MultiVectorBase<Scalar>> tmp = createMembers(x_space, num_g);
      W_op->apply(Thyra::NOTRANS, *adjoint_mv, tmp.ptr(), Scalar(1.0),
                  Scalar(0.0));
      dxdp_init_->apply(Thyra::CONJTRANS, *tmp, dgdp_.ptr(), Scalar(1.0),
                        Scalar(1.0));
    }
  }

  // Add in model parameter term = \int_0^T( (df/dp^T(t)*y(t) )dt which
  // is computed during the adjoint integration as an auxiliary integral
  // (2nd block of the solution vector).  For a time-integrated response, the
  // auxiliary integral also includes -\int_0^T dg/dp^T(t) dt.
  if (f_depends_on_p_ || (response_type_ == AdjointResponseType::TimeIntegral && g_depends_on_p_)) {
    RCP<const MultiVectorBase<Scalar>> z_mv =
        rcp_dynamic_cast<const DMVPV>(adjoint_state->getVectorBlock(1))
            ->getMultiVector();
    Thyra::V_VmV(dgdp_.ptr(), *dgdp_, *z_mv);
  }

  buildSolutionHistory(state_solution_history, adjoint_solution_history);

  return state_status && sens_status;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::createAdjointMassMatrixSolver(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs) const
{
  return AdjointSensitivityUtils::createAdjointLinearOpWithSolve(
      *adjoint_model_, inArgs, Scalar(1.0), Scalar(0.0));
}

template <class Scalar>
Teuchos::RCP<const SolutionState<Scalar>>
IntegratorAdjointSensitivity<Scalar>::findForwardState(
    const SolutionHistory<Scalar>& state_solution_history, const Scalar t, const Scalar tol) const
{
  Teuchos::RCP<const SolutionState<Scalar>> closest;
  for (int i = 0; i < state_solution_history.getNumStates(); ++i) {
    Teuchos::RCP<const SolutionState<Scalar>> state = state_solution_history[i];
    if (closest == Teuchos::null || std::abs(state->getTime() - t) < std::abs(closest->getTime() - t))
      closest = state;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      closest == Teuchos::null || std::abs(closest->getTime() - t) > tol,
      std::logic_error,
      "Error - IntegratorAdjointSensitivity: no stored time step of the "
      "forward integration ends at the response time "
          << t << " (closest stored time = "
          << (closest == Teuchos::null ? t : closest->getTime())
          << ").\n"
             "  With \"Step Type\" = \"Variable\", the time steps land on the "
             "response times.  With \"Step Type\" = \"Constant\", the response "
             "times must be on the time grid.  The forward solution history must "
             "store all time steps (\"Storage Type\" = \"Unlimited\").\n");
  return closest;
}

template <class Scalar>
bool IntegratorAdjointSensitivity<Scalar>::integrateAdjoint(
    const Scalar tau_begin, const Scalar tau_end,
    Thyra::VectorBase<Scalar>& adjoint_state, int& index)
{
  adjoint_integrator_->initializeSolutionHistory(tau_begin, adjoint_state.clone_v());
  // Continue the time step index of the previous segment, "Constant" time
  // steps compute the time from the index.
  if (index >= 0)
    adjoint_integrator_->getNonConstSolutionHistory()->getCurrentState()->setIndex(index);
  // Earlier segments may have reduced the final time of the time step control
  Teuchos::RCP<TimeStepControl<Scalar>> tsc = adjoint_integrator_->getNonConstTimeStepControl();
  if (tau_end > tsc->getFinalTime()) tsc->setFinalTime(tau_end);
  const bool status = adjoint_integrator_->advanceTime(tau_end);
  Thyra::assign(Teuchos::ptrFromRef(adjoint_state), *(adjoint_integrator_->getX()));
  index = adjoint_integrator_->getIndex();
  return status;
}

template <class Scalar>
void IntegratorAdjointSensitivity<Scalar>::appendAdjointSolutionHistory(
    const Scalar tol)
{
  Teuchos::RCP<const SolutionHistory<Scalar>> sh = adjoint_integrator_->getSolutionHistory();
  if (adjoint_solution_history_ == Teuchos::null) {
    adjoint_solution_history_ = Teuchos::rcp(new SolutionHistory<Scalar>());
    adjoint_solution_history_->setName("Adjoint Solution History");
    adjoint_solution_history_->setStorageType(STORAGE_TYPE_UNLIMITED);
  }
  for (int i = 0; i < sh->getNumStates(); ++i) {
    Teuchos::RCP<const SolutionState<Scalar>> state = (*sh)[i];
    // The first state of a segment is the last state of the previous segment
    if (adjoint_solution_history_->getNumStates() > 0 &&
        state->getTime() <= adjoint_solution_history_->maxTime() + tol)
      continue;
    adjoint_solution_history_->addState(state->clone());
  }
}

template <class Scalar>
Scalar IntegratorAdjointSensitivity<Scalar>::getTime() const
{
  // The combined solution history is empty if the forward integration failed
  if (solutionHistory_->getNumStates() == 0)
    return state_integrator_->getTime();
  return solutionHistory_->getCurrentTime();
}

template <class Scalar>
int IntegratorAdjointSensitivity<Scalar>::getIndex() const
{
  if (solutionHistory_->getNumStates() == 0)
    return state_integrator_->getIndex();
  return solutionHistory_->getCurrentIndex();
}

template <class Scalar>
Status IntegratorAdjointSensitivity<Scalar>::getStatus() const
{
  Status state_status = state_integrator_->getStatus();
  Status sens_status  = adjoint_integrator_->getStatus();
  if (state_status == FAILED || sens_status == FAILED) return FAILED;
  if (state_status == WORKING || sens_status == WORKING) return WORKING;
  return PASSED;
}

template <class Scalar>
void IntegratorAdjointSensitivity<Scalar>::setStatus(const Status st)
{
  state_integrator_->setStatus(st);
  adjoint_integrator_->setStatus(st);
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar>> IntegratorAdjointSensitivity<Scalar>::getStepper()
    const
{
  return state_integrator_->getStepper();
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getSolutionHistory() const
{
  return solutionHistory_;
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getStateSolutionHistory() const
{
  return state_integrator_->getSolutionHistory();
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getSensSolutionHistory() const
{
  // For the "Sum Over Times" response, the adjoint integration is restarted at
  // the response times and the segments are combined in
  // adjoint_solution_history_.
  if (adjoint_solution_history_ != Teuchos::null)
    return adjoint_solution_history_;
  return adjoint_integrator_->getSolutionHistory();
}

template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getNonConstSolutionHistory()
{
  return solutionHistory_;
}

template <class Scalar>
Teuchos::RCP<const TimeStepControl<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getTimeStepControl() const
{
  return state_integrator_->getTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getStateNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getSensNonConstTimeStepControl()
{
  return adjoint_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
void IntegratorAdjointSensitivity<Scalar>::initializeSolutionHistory(
    Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdot0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdotdot0,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxDp0,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> /* DxdotDp0 */,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> /* DxdotdotDp0 */)
{
  state_integrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0);
  dxdp_init_ = DxDp0;
}

template <class Scalar>
Teuchos::RCP<IntegratorObserver<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getObserver()
{
  return state_integrator_->getObserver();
}

template <class Scalar>
void IntegratorAdjointSensitivity<Scalar>::setObserver(
    Teuchos::RCP<IntegratorObserver<Scalar>> obs)
{
  state_integrator_->setObserver(obs);
  // ETP 1/12/22 Disabling passing of the observer to the adjoint
  // integrator to work around issues in Piro
  // adjoint_integrator_->setObserver(obs);
}

template <class Scalar>
void IntegratorAdjointSensitivity<Scalar>::initialize()
{
  state_integrator_->initialize();
  adjoint_integrator_->initialize();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getX() const
{
  return state_integrator_->getX();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getXDot() const
{
  return state_integrator_->getXDot();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getXDotDot() const
{
  return state_integrator_->getXDotDot();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getY() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultProductVector<Scalar> DPV;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  RCP<const DPV> pv     = rcp_dynamic_cast<const DPV>(adjoint_integrator_->getX());
  RCP<const DMVPV> mvpv = rcp_dynamic_cast<const DMVPV>(pv->getVectorBlock(0));
  return mvpv->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getYDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultProductVector<Scalar> DPV;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  RCP<const DPV> pv =
      rcp_dynamic_cast<const DPV>(adjoint_integrator_->getXDot());
  RCP<const DMVPV> mvpv = rcp_dynamic_cast<const DMVPV>(pv->getVectorBlock(0));
  return mvpv->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getYDotDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultProductVector<Scalar> DPV;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  RCP<const DPV> pv =
      rcp_dynamic_cast<const DPV>(adjoint_integrator_->getXDotDot());
  RCP<const DMVPV> mvpv = rcp_dynamic_cast<const DMVPV>(pv->getVectorBlock(0));
  return mvpv->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorAdjointSensitivity<Scalar>::getDgDp() const
{
  return dgdp_;
}

template <class Scalar>
std::string IntegratorAdjointSensitivity<Scalar>::description() const
{
  std::string name = "Tempus::IntegratorAdjointSensitivity";
  return (name);
}

template <class Scalar>
void IntegratorAdjointSensitivity<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << description() << "::describe" << std::endl;
  state_integrator_->describe(*l_out, verbLevel);
  adjoint_integrator_->describe(*l_out, verbLevel);
}

template <class Scalar>
SensitivityStepMode IntegratorAdjointSensitivity<Scalar>::getStepMode() const
{
  return stepMode_;
}

template <class Scalar>
Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar>>
IntegratorAdjointSensitivity<Scalar>::createAdjointModel(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_model,
    const Teuchos::RCP<Teuchos::ParameterList>& inputPL)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<ParameterList> spl = Teuchos::parameterList();
  if (inputPL != Teuchos::null) {
    *spl = inputPL->sublist("Sensitivities");
  }
  if (spl->isParameter("IC Depends on Parameters"))
    spl->remove("IC Depends on Parameters");
  spl->remove(DAEStructure<Scalar>::algebraicMaskName(), false);
  spl->remove(AdjointDAEJumpSolver<Scalar>::sublistName(), false);

  const Scalar tinit  = state_integrator_->getTimeStepControl()->getInitTime();
  const Scalar tfinal = state_integrator_->getTimeStepControl()->getFinalTime();
  return rcp(new AdjointAuxSensitivityModelEvaluator<Scalar>(
      model, adjoint_model, tinit, tfinal, spl));
}

template <class Scalar>
void IntegratorAdjointSensitivity<Scalar>::buildSolutionHistory(
    const Teuchos::RCP<const SolutionHistory<Scalar>>& state_solution_history,
    const Teuchos::RCP<const SolutionHistory<Scalar>>& adjoint_solution_history)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::createMembers;
  using Thyra::MultiVectorBase;
  using Thyra::multiVectorProductVector;
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;
  typedef Thyra::DefaultProductVectorSpace<Scalar> DPVS;
  typedef Thyra::DefaultProductVector<Scalar> DPV;

  RCP<const VectorSpaceBase<Scalar>> x_space = model_->get_x_space();
  RCP<const VectorSpaceBase<Scalar>> adjoint_space =
      rcp_dynamic_cast<const DPVS>(adjoint_aux_model_->get_x_space())
          ->getBlock(0);
  Teuchos::Array<RCP<const VectorSpaceBase<Scalar>>> spaces(2);
  spaces[0]                  = x_space;
  spaces[1]                  = adjoint_space;
  RCP<const DPVS> prod_space = Thyra::productVectorSpace(spaces());

  int num_states        = state_solution_history->getNumStates();
  const Scalar t_init   = state_integrator_->getTimeStepControl()->getInitTime();
  const Scalar t_final  = state_integrator_->getTime();
  const Scalar time_tol = AdjointSensitivityUtils::timeTolerance(
      t_init, t_final, time_tolerance_pl_.get());
  for (int i = 0; i < num_states; ++i) {
    RCP<const SolutionState<Scalar>> forward_state =
        (*state_solution_history)[i];
    // The adjoint time steps do not coincide with the forward ones when the
    // time step size varies (e.g., when the forward time steps land on the
    // times of a "Sum Over Times" response), use interpolation in that case.
    const Scalar adjoint_time = t_final + t_init - forward_state->getTime();
    RCP<const SolutionState<Scalar>> adjoint_state;
    {
      // Binary search for the first adjoint state with time >= adjoint_time-tol
      int lo = 0, hi = adjoint_solution_history->getNumStates();
      while (lo < hi) {
        const int mid = (lo + hi) / 2;
        if ((*adjoint_solution_history)[mid]->getTime() <
            adjoint_time - time_tol)
          lo = mid + 1;
        else
          hi = mid;
      }
      if (lo < adjoint_solution_history->getNumStates() &&
          std::abs((*adjoint_solution_history)[lo]->getTime() -
                   adjoint_time) <= time_tol)
        adjoint_state = (*adjoint_solution_history)[lo];
    }
    if (adjoint_state == Teuchos::null)
      adjoint_state = adjoint_solution_history->interpolateState(adjoint_time);

    // X
    RCP<DPV> x = Thyra::defaultProductVector(prod_space);
    RCP<const VectorBase<Scalar>> adjoint_x =
        rcp_dynamic_cast<const DPV>(adjoint_state->getX())->getVectorBlock(0);
    assign(x->getNonconstVectorBlock(0).ptr(), *(forward_state->getX()));
    assign(x->getNonconstVectorBlock(1).ptr(), *(adjoint_x));
    RCP<VectorBase<Scalar>> x_b = x;

    // X-Dot
    RCP<DPV> x_dot = Thyra::defaultProductVector(prod_space);
    RCP<const VectorBase<Scalar>> adjoint_x_dot =
        rcp_dynamic_cast<const DPV>(adjoint_state->getXDot())
            ->getVectorBlock(0);
    assign(x_dot->getNonconstVectorBlock(0).ptr(), *(forward_state->getXDot()));
    assign(x_dot->getNonconstVectorBlock(1).ptr(), *(adjoint_x_dot));
    RCP<VectorBase<Scalar>> x_dot_b = x_dot;

    // X-Dot-Dot
    RCP<DPV> x_dot_dot;
    if (forward_state->getXDotDot() != Teuchos::null) {
      x_dot_dot = Thyra::defaultProductVector(prod_space);
      RCP<const VectorBase<Scalar>> adjoint_x_dot_dot =
          rcp_dynamic_cast<const DPV>(adjoint_state->getXDotDot())
              ->getVectorBlock(0);
      assign(x_dot_dot->getNonconstVectorBlock(0).ptr(),
             *(forward_state->getXDotDot()));
      assign(x_dot_dot->getNonconstVectorBlock(1).ptr(), *(adjoint_x_dot_dot));
    }
    RCP<VectorBase<Scalar>> x_dot_dot_b = x_dot_dot;

    RCP<SolutionState<Scalar>> prod_state = forward_state->clone();
    prod_state->setX(x_b);
    prod_state->setXDot(x_dot_b);
    prod_state->setXDotDot(x_dot_dot_b);
    prod_state->setPhysicsState(Teuchos::null);
    solutionHistory_->addState(prod_state);
  }
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar>>
createIntegratorAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> inputPL,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_model)
{
  // set the parameters
  Teuchos::RCP<Teuchos::ParameterList> spl = Teuchos::parameterList();
  if (inputPL != Teuchos::null) *spl = inputPL->sublist("Sensitivities");
  // Tolerances for comparing times (see AdjointSensitivityUtils::timeTolerance)
  const Teuchos::RCP<const Teuchos::ParameterList> time_tolerance_pl =
      Teuchos::rcp(new Teuchos::ParameterList(*spl));

  int p_index          = spl->get<int>("Sensitivity Parameter Index", 0);
  int g_index          = spl->get<int>("Response Function Index", 0);
  bool g_depends_on_p  = spl->get<bool>("Response Depends on Parameters", true);
  bool f_depends_on_p  = spl->get<bool>("Residual Depends on Parameters", true);
  bool ic_depends_on_p = spl->get<bool>("IC Depends on Parameters", true);
  bool mass_matrix_is_identity =
      spl->get<bool>("Mass Matrix Is Identity", false);
  const Teuchos::RCP<const DAEStructure<Scalar>> dae_structure = getDAEStructure<Scalar>(*spl);
  Teuchos::RCP<const Teuchos::ParameterList> dae_jump_solver_pl;
  if (spl->isSublist(AdjointDAEJumpSolver<Scalar>::sublistName()))
    dae_jump_solver_pl = Teuchos::rcp(new Teuchos::ParameterList(spl->sublist(AdjointDAEJumpSolver<Scalar>::sublistName())));
  const AdjointResponseType response_type =
      AdjointSensitivityUtils::parseResponseType(
          spl->get<std::string>("Response Type", "Final Time"));
  const std::vector<Scalar> response_times =
      AdjointSensitivityUtils::parseTimeList<Scalar>(
          spl->get<std::string>("Response Time List", ""));
  TEUCHOS_TEST_FOR_EXCEPTION(
      response_type == AdjointResponseType::SumOverTimes &&
          response_times.empty(),
      std::logic_error,
      "Error - createIntegratorAdjointSensitivity(): the \"Sum Over Times\" "
      "response requires a non-empty \"Response Time List\".\n");

  auto state_integrator = createIntegratorBasic<Scalar>(inputPL, model);

  // The forward time steps must end at the response times
  if (response_type == AdjointResponseType::SumOverTimes)
    AdjointSensitivityUtils::addResponseTimeEvents(
        *state_integrator->getNonConstTimeStepControl(), response_times);

  // createAdjointModel
  if (spl->isParameter("IC Depends on Parameters"))
    spl->remove("IC Depends on Parameters");
  // The DAE structure is used by the integrator, not by the adjoint model
  spl->remove(DAEStructure<Scalar>::algebraicMaskName(), false);
  spl->remove(AdjointDAEJumpSolver<Scalar>::sublistName(), false);

  const Scalar tinit  = state_integrator->getTimeStepControl()->getInitTime();
  const Scalar tfinal = state_integrator->getTimeStepControl()->getFinalTime();
  // auto adjoint_model  = Teuchos::rcp(new
  // AdjointAuxSensitivityModelEvaluator<Scalar>(model, tfinal, spl));
  // TODO: where is the adjoint ME coming from?

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> adjt_model = adjoint_model;
  if (adjoint_model == Teuchos::null)
    adjt_model = Thyra::implicitAdjointModelEvaluator(model);

  auto adjoint_aux_model =
      Teuchos::rcp(new AdjointAuxSensitivityModelEvaluator<Scalar>(
          model, adjt_model, tinit, tfinal, spl));

  // Create combined solution histories combining the forward and adjoint
  // solutions.  We do not include the auxiliary part from the adjoint solution.
  auto integrator_name           = inputPL->get<std::string>("Integrator Name");
  auto integratorPL              = Teuchos::sublist(inputPL, integrator_name, true);
  auto shPL                      = Teuchos::sublist(integratorPL, "Solution History", true);
  auto combined_solution_History = createSolutionHistoryPL<Scalar>(shPL);

  // The adjoint integrator uses a copy of the parameters, without the objects
  // of the forward problem that the nonlinear solver of the forward stepper
  // stores in them (e.g., for NOX row sum scaling), which do not apply to the
  // adjoint problem (see removeIntegratorSolverObjects).
  auto adjointPL = Teuchos::rcp(new Teuchos::ParameterList(*inputPL));
  removeIntegratorSolverObjects(*adjointPL);
  auto adjoint_integrator =
      createIntegratorBasic<Scalar>(adjointPL, adjoint_aux_model);

  Teuchos::RCP<IntegratorAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorAdjointSensitivity<Scalar>(
          model, state_integrator, adjt_model, adjoint_aux_model,
          adjoint_integrator, combined_solution_History, p_index, g_index,
          g_depends_on_p, f_depends_on_p, ic_depends_on_p,
          mass_matrix_is_identity, response_type, response_times,
          dae_structure, dae_jump_solver_pl));
  integrator->setTimeToleranceParameters(time_tolerance_pl);

  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar>>
createIntegratorAdjointSensitivity()
{
  Teuchos::RCP<IntegratorAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorAdjointSensitivity<Scalar>());
  return (integrator);
}

}  // namespace Tempus
#endif  // Tempus_IntegratorAdjointSensitivity_impl_hpp
