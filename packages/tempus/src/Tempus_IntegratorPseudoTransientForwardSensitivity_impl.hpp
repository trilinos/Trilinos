//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorPseudoTransientForwardSensitivity_impl_hpp
#define Tempus_IntegratorPseudoTransientForwardSensitivity_impl_hpp

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "Tempus_WrapStaggeredFSAModelEvaluator.hpp"

namespace Tempus {

template <class Scalar>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
    IntegratorPseudoTransientForwardSensitivity(
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
        const Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>> &sens_model,
        const Teuchos::RCP<IntegratorBasic<Scalar>> &fwd_integrator,
        const Teuchos::RCP<IntegratorBasic<Scalar>> &sens_integrator,
        const bool reuse_solver, const bool force_W_update)
  : model_(model),
    sens_model_(sens_model),
    state_integrator_(fwd_integrator),
    sens_integrator_(sens_integrator),
    reuse_solver_(reuse_solver),
    force_W_update_(force_W_update),
    stepMode_(SensitivityStepMode::Forward)
{
}

template <class Scalar>
IntegratorPseudoTransientForwardSensitivity<
    Scalar>::IntegratorPseudoTransientForwardSensitivity()
  : reuse_solver_(false),
    force_W_update_(false),
    stepMode_(SensitivityStepMode::Forward)
{
  state_integrator_ = createIntegratorBasic<Scalar>();
  sens_integrator_  = createIntegratorBasic<Scalar>();
}

template <class Scalar>
bool IntegratorPseudoTransientForwardSensitivity<Scalar>::advanceTime()
{
  using Teuchos::RCP;
  using Thyra::VectorBase;

  // Run state integrator and get solution
  stepMode_         = SensitivityStepMode::Forward;
  bool state_status = state_integrator_->advanceTime();

  // Set solution in sensitivity ME
  sens_model_->setForwardSolutionState(state_integrator_->getCurrentState());

  // Reuse state solver if requested
  if (reuse_solver_ &&
      state_integrator_->getStepper()->getSolver() != Teuchos::null) {
    sens_model_->setSolver(state_integrator_->getStepper()->getSolver(),
                           force_W_update_);
  }

  // Run sensitivity integrator
  stepMode_        = SensitivityStepMode::Sensitivity;
  bool sens_status = sens_integrator_->advanceTime();

  buildSolutionHistory();

  return state_status && sens_status;
}

template <class Scalar>
bool IntegratorPseudoTransientForwardSensitivity<Scalar>::advanceTime(
    const Scalar timeFinal)
{
  TEMPUS_FUNC_TIME_MONITOR_DIFF(
      "Tempus::IntegratorPseudoTransientForwardSensitivity::advanceTime()",
      TEMPUS_PTFS_AT);

  using Teuchos::RCP;
  using Thyra::VectorBase;

  // Run state integrator and get solution
  bool state_status = true;
  {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::IntegratorPseudoTransientForwardSensitivity::advanceTime::"
        "state",
        TEMPUS_PTFS_AT_FWD);
    state_status = state_integrator_->advanceTime(timeFinal);
  }

  // Set solution in sensitivity ME
  sens_model_->setForwardSolutionState(state_integrator_->getCurrentState());

  // Reuse state solver if requested
  if (reuse_solver_ &&
      state_integrator_->getStepper()->getSolver() != Teuchos::null) {
    sens_model_->setSolver(state_integrator_->getStepper()->getSolver(),
                           force_W_update_);
  }

  // Run sensitivity integrator
  bool sens_status = true;
  {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::IntegratorPseudoTransientForwardSensitivity::advanceTime::"
        "sensitivity",
        TEMPUS_PTFS_AT_SEN);
    sens_status = sens_integrator_->advanceTime(timeFinal);
  }

  buildSolutionHistory();

  return state_status && sens_status;
}

template <class Scalar>
Scalar IntegratorPseudoTransientForwardSensitivity<Scalar>::getTime() const
{
  return solutionHistory_->getCurrentTime();
}

template <class Scalar>
int IntegratorPseudoTransientForwardSensitivity<Scalar>::getIndex() const
{
  return solutionHistory_->getCurrentIndex();
}

template <class Scalar>
Status IntegratorPseudoTransientForwardSensitivity<Scalar>::getStatus() const
{
  Status state_status = state_integrator_->getStatus();
  Status sens_status  = sens_integrator_->getStatus();
  if (state_status == FAILED || sens_status == FAILED) return FAILED;
  if (state_status == WORKING || sens_status == WORKING) return WORKING;
  return PASSED;
}

template <class Scalar>
void IntegratorPseudoTransientForwardSensitivity<Scalar>::setStatus(
    const Status st)
{
  state_integrator_->setStatus(st);
  sens_integrator_->setStatus(st);
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getStepper() const
{
  return state_integrator_->getStepper();
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getStateStepper() const
{
  return state_integrator_->getStepper();
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getSensStepper() const
{
  return sens_integrator_->getStepper();
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getSolutionHistory() const
{
  return solutionHistory_;
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getStateSolutionHistory()
    const
{
  return state_integrator_->getSolutionHistory();
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getSensSolutionHistory()
    const
{
  return sens_integrator_->getSolutionHistory();
}

template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar>>
IntegratorPseudoTransientForwardSensitivity<
    Scalar>::getNonConstSolutionHistory()
{
  return solutionHistory_;
}

template <class Scalar>
Teuchos::RCP<const TimeStepControl<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getTimeStepControl() const
{
  return state_integrator_->getTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorPseudoTransientForwardSensitivity<
    Scalar>::getNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorPseudoTransientForwardSensitivity<
    Scalar>::getStateNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorPseudoTransientForwardSensitivity<
    Scalar>::getSensNonConstTimeStepControl()
{
  return sens_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<IntegratorObserver<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getObserver()
{
  return state_integrator_->getObserver();
}

template <class Scalar>
void IntegratorPseudoTransientForwardSensitivity<Scalar>::setObserver(
    Teuchos::RCP<IntegratorObserver<Scalar>> obs)
{
  state_integrator_->setObserver(obs);
  sens_integrator_->setObserver(obs);
}

template <class Scalar>
void IntegratorPseudoTransientForwardSensitivity<Scalar>::
    initializeSolutionHistory(
        Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0,
        Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdot0,
        Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdotdot0,
        Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxDp0,
        Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxdotDp0,
        Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxdotdotDp0)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::createMember;
  using Thyra::VectorSpaceBase;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  //
  // Create and initialize product X, Xdot, Xdotdot

  RCP<const VectorSpaceBase<Scalar>> space = sens_model_->get_x_space();
  RCP<DMVPV> X                             = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdot                          = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdotdot                       = rcp_dynamic_cast<DMVPV>(createMember(space));
  const Scalar zero                        = Teuchos::ScalarTraits<Scalar>::zero();

  // x
  if (DxDp0 == Teuchos::null)
    assign(X->getNonconstMultiVector().ptr(), zero);
  else
    assign(X->getNonconstMultiVector().ptr(), *DxDp0);

  // xdot
  if (DxdotDp0 == Teuchos::null)
    assign(Xdot->getNonconstMultiVector().ptr(), zero);
  else
    assign(Xdot->getNonconstMultiVector().ptr(), *DxdotDp0);

  // xdotdot
  if (DxdotDp0 == Teuchos::null)
    assign(Xdotdot->getNonconstMultiVector().ptr(), zero);
  else
    assign(Xdotdot->getNonconstMultiVector().ptr(), *DxdotdotDp0);

  state_integrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0);
  sens_integrator_->initializeSolutionHistory(t0, X, Xdot, Xdotdot);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getX() const
{
  return state_integrator_->getX();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getDxDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> X = rcp_dynamic_cast<const DMVPV>(sens_integrator_->getX());
  return X->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getXDot() const
{
  return state_integrator_->getXDot();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getDXDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdot =
      rcp_dynamic_cast<const DMVPV>(sens_integrator_->getXDot());
  return Xdot->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getXDotDot() const
{
  return state_integrator_->getXDotDot();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getDXDotDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdotdot =
      rcp_dynamic_cast<const DMVPV>(sens_integrator_->getXDotDot());
  return Xdotdot->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getG() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // Compute g which is computed by response 1 of the
  // sensitivity model evaluator
  MEB::InArgs<Scalar> inargs   = sens_model_->getNominalValues();
  MEB::OutArgs<Scalar> outargs = sens_model_->createOutArgs();
  inargs.set_t(sens_integrator_->getTime());
  inargs.set_x(sens_integrator_->getX());
  if (inargs.supports(MEB::IN_ARG_x_dot))
    inargs.set_x_dot(sens_integrator_->getXDot());
  if (inargs.supports(MEB::IN_ARG_x_dot_dot))
    inargs.set_x_dot_dot(sens_integrator_->getXDotDot());

  Teuchos::RCP<Thyra::VectorBase<Scalar>> g =
      Thyra::createMember(sens_model_->get_g_space(1));
  outargs.set_g(1, g);

  sens_model_->evalModel(inargs, outargs);
  return g;
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientForwardSensitivity<Scalar>::getDgDp() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  // Compute final dg/dp  which is computed by response 0  of the
  // sensitivity model evaluator
  MEB::InArgs<Scalar> inargs   = sens_model_->getNominalValues();
  MEB::OutArgs<Scalar> outargs = sens_model_->createOutArgs();
  inargs.set_t(sens_integrator_->getTime());
  inargs.set_x(sens_integrator_->getX());
  if (inargs.supports(MEB::IN_ARG_x_dot))
    inargs.set_x_dot(sens_integrator_->getXDot());
  if (inargs.supports(MEB::IN_ARG_x_dot_dot))
    inargs.set_x_dot_dot(sens_integrator_->getXDotDot());

  Teuchos::RCP<Thyra::VectorBase<Scalar>> G =
      Thyra::createMember(sens_model_->get_g_space(0));
  Teuchos::RCP<DMVPV> dgdp = Teuchos::rcp_dynamic_cast<DMVPV>(G);
  outargs.set_g(0, G);

  sens_model_->evalModel(inargs, outargs);
  return dgdp->getMultiVector();
}

template <class Scalar>
std::string IntegratorPseudoTransientForwardSensitivity<Scalar>::description()
    const
{
  std::string name = "Tempus::IntegratorPseudoTransientForwardSensitivity";
  return (name);
}

template <class Scalar>
void IntegratorPseudoTransientForwardSensitivity<Scalar>::describe(
    Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << description() << "::describe" << std::endl;
  state_integrator_->describe(*l_out, verbLevel);
  sens_integrator_->describe(*l_out, verbLevel);
}

template <class Scalar>
SensitivityStepMode
IntegratorPseudoTransientForwardSensitivity<Scalar>::getStepMode() const
{
  return stepMode_;
}

template <class Scalar>
void IntegratorPseudoTransientForwardSensitivity<Scalar>::buildSolutionHistory()
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
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  // TODO: get the solution history PL or create it

  // Create combined solution histories, first for the states with zero
  // sensitivities and then for the sensitivities with frozen states
  RCP<ParameterList> shPL;
  // Teuchos::sublist(state_integrator_->getIntegratorParameterList(), "Solution
  // History", true);
  solutionHistory_ = createSolutionHistoryPL<Scalar>(shPL);

  const int num_param = rcp_dynamic_cast<const DMVPV>(sens_integrator_->getX())
                            ->getMultiVector()
                            ->domain()
                            ->dim();
  RCP<const VectorSpaceBase<Scalar>> x_space = model_->get_x_space();
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar>> prod_space =
      Thyra::multiVectorProductVectorSpace(x_space, num_param + 1);
  const Teuchos::Range1D rng(1, num_param);
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  RCP<const SolutionHistory<Scalar>> state_solution_history =
      state_integrator_->getSolutionHistory();
  int num_states = state_solution_history->getNumStates();
  for (int i = 0; i < num_states; ++i) {
    RCP<const SolutionState<Scalar>> state = (*state_solution_history)[i];

    // X
    RCP<MultiVectorBase<Scalar>> x_mv = createMembers(x_space, num_param + 1);
    assign(x_mv->col(0).ptr(), *(state->getX()));
    assign(x_mv->subView(rng).ptr(), zero);
    RCP<VectorBase<Scalar>> x = multiVectorProductVector(prod_space, x_mv);

    // X-Dot
    RCP<VectorBase<Scalar>> x_dot;
    if (state->getXDot() != Teuchos::null) {
      RCP<MultiVectorBase<Scalar>> x_dot_mv =
          createMembers(x_space, num_param + 1);
      assign(x_dot_mv->col(0).ptr(), *(state->getXDot()));
      assign(x_dot_mv->subView(rng).ptr(), zero);
      x_dot = multiVectorProductVector(prod_space, x_dot_mv);
    }

    // X-Dot-Dot
    RCP<VectorBase<Scalar>> x_dot_dot;
    if (state->getXDotDot() != Teuchos::null) {
      RCP<MultiVectorBase<Scalar>> x_dot_dot_mv =
          createMembers(x_space, num_param + 1);
      assign(x_dot_dot_mv->col(0).ptr(), *(state->getXDotDot()));
      assign(x_dot_dot_mv->subView(rng).ptr(), zero);
      x_dot_dot = multiVectorProductVector(prod_space, x_dot_dot_mv);
    }

    RCP<SolutionState<Scalar>> prod_state = state->clone();
    prod_state->setX(x);
    prod_state->setXDot(x_dot);
    prod_state->setXDotDot(x_dot_dot);
    solutionHistory_->addState(prod_state);
  }

  RCP<const VectorBase<Scalar>> frozen_x =
      state_solution_history->getCurrentState()->getX();
  RCP<const VectorBase<Scalar>> frozen_x_dot =
      state_solution_history->getCurrentState()->getXDot();
  RCP<const VectorBase<Scalar>> frozen_x_dot_dot =
      state_solution_history->getCurrentState()->getXDotDot();
  RCP<const SolutionHistory<Scalar>> sens_solution_history =
      sens_integrator_->getSolutionHistory();
  num_states = sens_solution_history->getNumStates();
  for (int i = 0; i < num_states; ++i) {
    RCP<const SolutionState<Scalar>> state = (*sens_solution_history)[i];

    // X
    RCP<MultiVectorBase<Scalar>> x_mv = createMembers(x_space, num_param + 1);
    RCP<const MultiVectorBase<Scalar>> dxdp =
        rcp_dynamic_cast<const DMVPV>(state->getX())->getMultiVector();
    assign(x_mv->col(0).ptr(), *(frozen_x));
    assign(x_mv->subView(rng).ptr(), *dxdp);
    RCP<VectorBase<Scalar>> x = multiVectorProductVector(prod_space, x_mv);

    // X-Dot
    RCP<VectorBase<Scalar>> x_dot;
    if (state->getXDot() != Teuchos::null) {
      RCP<MultiVectorBase<Scalar>> x_dot_mv =
          createMembers(x_space, num_param + 1);
      RCP<const MultiVectorBase<Scalar>> dxdotdp =
          rcp_dynamic_cast<const DMVPV>(state->getXDot())->getMultiVector();
      assign(x_dot_mv->col(0).ptr(), *(frozen_x_dot));
      assign(x_dot_mv->subView(rng).ptr(), *dxdotdp);
      x_dot = multiVectorProductVector(prod_space, x_dot_mv);
    }

    // X-Dot-Dot
    RCP<VectorBase<Scalar>> x_dot_dot;
    if (state->getXDotDot() != Teuchos::null) {
      RCP<MultiVectorBase<Scalar>> x_dot_dot_mv =
          createMembers(x_space, num_param + 1);
      RCP<const MultiVectorBase<Scalar>> dxdotdotdp =
          rcp_dynamic_cast<const DMVPV>(state->getXDotDot())->getMultiVector();
      assign(x_dot_dot_mv->col(0).ptr(), *(frozen_x_dot_dot));
      assign(x_dot_dot_mv->subView(rng).ptr(), *dxdotdotdp);
      x_dot_dot = multiVectorProductVector(prod_space, x_dot_dot_mv);
    }

    RCP<SolutionState<Scalar>> prod_state = state->clone();
    prod_state->setX(x);
    prod_state->setXDot(x_dot);
    prod_state->setXDotDot(x_dot_dot);
    solutionHistory_->addState(prod_state, false);
  }
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
createIntegratorPseudoTransientForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &sens_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &sens_solve_model)
{
  auto fwd_integrator = createIntegratorBasic<Scalar>(pList, model);
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>> sens_model;
  Teuchos::RCP<IntegratorBasic<Scalar>> sens_integrator;

  {
    Teuchos::RCP<Teuchos::ParameterList> pl =
        Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::RCP<const Teuchos::ParameterList> integrator_pl =
        fwd_integrator->getValidParameters();
    Teuchos::RCP<const Teuchos::ParameterList> sensitivity_pl =
        StaggeredForwardSensitivityModelEvaluator<Scalar>::getValidParameters();
    pl->setParameters(*integrator_pl);
    pl->sublist("Sensitivities").setParameters(*sensitivity_pl);
    pl->sublist("Sensitivities").set("Reuse State Linear Solver", false);
    pl->sublist("Sensitivities").set("Force W Update", false);
    pl->sublist("Sensitivities").set("Cache Matrices", false);
    pList->setParametersNotAlreadySet(*pl);
  }

  bool reuse_solver =
      pList->sublist("Sensitivities").get("Reuse State Linear Solver", false);
  bool force_W_update =
      pList->sublist("Sensitivities").get("Force W Update", false);
  bool cache_matrices =
      pList->sublist("Sensitivities").get("Cache Matrices", false);

  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList != Teuchos::null) {
      *pl = pList->sublist("Sensitivities");
    }
    pl->remove("Reuse State Linear Solver");
    pl->remove("Force W Update");
    pl->remove("Cache Matrices");
    sens_model = wrapStaggeredFSAModelEvaluator(
        model, sens_residual_model, sens_solve_model, cache_matrices, pl);
    sens_integrator = createIntegratorBasic<Scalar>(pList, sens_model);
  }

  Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
      integrator = Teuchos::rcp(
          new Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>(
              model, sens_model, fwd_integrator, sens_integrator, reuse_solver,
              force_W_update));

  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
createIntegratorPseudoTransientForwardSensitivity()
{
  Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
      integrator = Teuchos::rcp(
          new Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>());
  return (integrator);
}

}  // namespace Tempus
#endif  // Tempus_IntegratorPseudoTransientForwardSensitivity_impl_hpp
