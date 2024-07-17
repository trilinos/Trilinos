//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorPseudoTransientAdjointSensitivity_impl_hpp
#define Tempus_IntegratorPseudoTransientAdjointSensitivity_impl_hpp

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_ImplicitAdjointModelEvaluator.hpp"

namespace Tempus {

template <class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
    IntegratorPseudoTransientAdjointSensitivity(
        Teuchos::RCP<Teuchos::ParameterList> inputPL,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>&
            adjoint_residual_model,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_solve_model)
{
  model_                  = model;
  adjoint_residual_model_ = adjoint_residual_model;
  adjoint_solve_model_    = adjoint_solve_model;
  state_integrator_       = createIntegratorBasic<Scalar>(inputPL, model_);
  sens_model_             = createSensitivityModel(model_, adjoint_residual_model_,
                                                   adjoint_solve_model_, inputPL);
  sens_integrator_        = createIntegratorBasic<Scalar>(inputPL, sens_model_);
  stepMode_               = SensitivityStepMode::Forward;
  do_forward_integration_ = true;
  do_adjoint_integration_ = true;
}

template <class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
    IntegratorPseudoTransientAdjointSensitivity(
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>&
            adjoint_residual_model,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_solve_model,
        std::string stepperType)
{
  model_                  = model;
  adjoint_residual_model_ = adjoint_residual_model;
  adjoint_solve_model_    = adjoint_solve_model;
  state_integrator_       = createIntegratorBasic<Scalar>(model_, stepperType);
  sens_model_             = createSensitivityModel(model_, adjoint_residual_model_,
                                                   adjoint_solve_model_, Teuchos::null);
  sens_integrator_        = createIntegratorBasic<Scalar>(sens_model_, stepperType);
  stepMode_               = SensitivityStepMode::Forward;
  do_forward_integration_ = true;
  do_adjoint_integration_ = true;
}

template <class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
    IntegratorPseudoTransientAdjointSensitivity(
        Teuchos::RCP<Teuchos::ParameterList> inputPL,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_model)
  : IntegratorPseudoTransientAdjointSensitivity(inputPL, model, adjoint_model,
                                                adjoint_model)
{
}

template <class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
    IntegratorPseudoTransientAdjointSensitivity(
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_model,
        std::string stepperType)
  : IntegratorPseudoTransientAdjointSensitivity(model, adjoint_model,
                                                adjoint_model, stepperType)
{
}

template <class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
    IntegratorPseudoTransientAdjointSensitivity(
        Teuchos::RCP<Teuchos::ParameterList> inputPL,
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model)
  : IntegratorPseudoTransientAdjointSensitivity(
        inputPL, model, Thyra::implicitAdjointModelEvaluator(model))
{
}

template <class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
    IntegratorPseudoTransientAdjointSensitivity(
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
        std::string stepperType)
  : IntegratorPseudoTransientAdjointSensitivity(
        model, Thyra::implicitAdjointModelEvaluator(model), stepperType)
{
}

template <class Scalar>
IntegratorPseudoTransientAdjointSensitivity<
    Scalar>::IntegratorPseudoTransientAdjointSensitivity()
{
  state_integrator_ = createIntegratorBasic<Scalar>();
  sens_integrator_  = createIntegratorBasic<Scalar>();
  stepMode_         = SensitivityStepMode::Forward;
}

template <class Scalar>
bool IntegratorPseudoTransientAdjointSensitivity<Scalar>::advanceTime()
{
  const Scalar tfinal = state_integrator_->getTimeStepControl()->getFinalTime();
  return advanceTime(tfinal);
}

template <class Scalar>
bool IntegratorPseudoTransientAdjointSensitivity<Scalar>::advanceTime(
    const Scalar timeFinal)
{
  TEMPUS_FUNC_TIME_MONITOR_DIFF(
      "Tempus::IntegratorPseudoTransientAdjointSensitivity::advanceTime()",
      TEMPUS_PTAS_AT);

  using Teuchos::RCP;
  using Thyra::VectorBase;
  typedef Thyra::ModelEvaluatorBase MEB;

  bool state_status = true;
  if (do_forward_integration_) {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::IntegratorPseudoTransientAdjointSensitivity::advanceTime::"
        "state",
        TEMPUS_PTAS_AT_FWD);

    // Run state integrator and get solution
    stepMode_    = SensitivityStepMode::Forward;
    state_status = state_integrator_->advanceTime(timeFinal);
  }

  bool sens_status = true;
  if (do_adjoint_integration_) {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::IntegratorPseudoTransientAdjointSensitivity::advanceTime::"
        "adjoint",
        TEMPUS_PTAS_AT_ADJ);

    // For at least some time-stepping methods, the time of the last time step
    // may not be timeFinal (e.g., it may be greater by at most delta_t).
    // But since the adjoint model requires timeFinal in its formulation, reset
    // it to the achieved final time.
    sens_model_->setFinalTime(state_integrator_->getTime());

    // Set solution in sensitivity ME
    sens_model_->setForwardSolutionHistory(
        state_integrator_->getSolutionHistory());

    // Run sensitivity integrator
    stepMode_   = SensitivityStepMode::Adjoint;
    sens_status = sens_integrator_->advanceTime(timeFinal);

    // Compute final dg/dp, g which is computed by response 0, 1 of the adjoint
    // model evaluator
    MEB::InArgs<Scalar> inargs   = sens_model_->getNominalValues();
    MEB::OutArgs<Scalar> outargs = sens_model_->createOutArgs();
    inargs.set_t(sens_integrator_->getTime());
    inargs.set_x(sens_integrator_->getX());
    if (inargs.supports(MEB::IN_ARG_x_dot))
      inargs.set_x_dot(sens_integrator_->getXDot());
    if (inargs.supports(MEB::IN_ARG_x_dot_dot))
      inargs.set_x_dot_dot(sens_integrator_->getXDotDot());
    RCP<VectorBase<Scalar>> G = dgdp_;
    if (G == Teuchos::null) {
      G     = Thyra::createMember(sens_model_->get_g_space(0));
      dgdp_ = Teuchos::rcp_dynamic_cast<DMVPV>(G);
    }
    if (g_ == Teuchos::null)
      g_ = Thyra::createMember(sens_model_->get_g_space(1));
    outargs.set_g(0, G);
    outargs.set_g(1, g_);
    sens_model_->evalModel(inargs, outargs);

    buildSolutionHistory();
  }

  return state_status && sens_status;
}

template <class Scalar>
Scalar IntegratorPseudoTransientAdjointSensitivity<Scalar>::getTime() const
{
  return solutionHistory_->getCurrentTime();
}

template <class Scalar>
int IntegratorPseudoTransientAdjointSensitivity<Scalar>::getIndex() const
{
  return solutionHistory_->getCurrentIndex();
}

template <class Scalar>
Status IntegratorPseudoTransientAdjointSensitivity<Scalar>::getStatus() const
{
  Status state_status = state_integrator_->getStatus();
  Status sens_status  = sens_integrator_->getStatus();
  if (state_status == FAILED || sens_status == FAILED) return FAILED;
  if (state_status == WORKING || sens_status == WORKING) return WORKING;
  return PASSED;
}

template <class Scalar>
void IntegratorPseudoTransientAdjointSensitivity<Scalar>::setStatus(
    const Status st)
{
  state_integrator_->setStatus(st);
  sens_integrator_->setStatus(st);
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getStepper() const
{
  return state_integrator_->getStepper();
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getStateStepper() const
{
  return state_integrator_->getStepper();
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getSensStepper() const
{
  return sens_integrator_->getStepper();
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getSolutionHistory() const
{
  return solutionHistory_;
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getStateSolutionHistory()
    const
{
  return state_integrator_->getSolutionHistory();
}

template <class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getSensSolutionHistory()
    const
{
  return sens_integrator_->getSolutionHistory();
}

template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<
    Scalar>::getNonConstSolutionHistory()
{
  return solutionHistory_;
}

template <class Scalar>
Teuchos::RCP<const TimeStepControl<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getTimeStepControl() const
{
  return state_integrator_->getTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<
    Scalar>::getNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<
    Scalar>::getStateNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<
    Scalar>::getSensNonConstTimeStepControl()
{
  return sens_integrator_->getNonConstTimeStepControl();
}

template <class Scalar>
Teuchos::RCP<IntegratorObserver<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getObserver()
{
  return state_integrator_->getObserver();
}

template <class Scalar>
void IntegratorPseudoTransientAdjointSensitivity<Scalar>::setObserver(
    Teuchos::RCP<IntegratorObserver<Scalar>> obs)
{
  state_integrator_->setObserver(obs);
  sens_integrator_->setObserver(obs);
}

template <class Scalar>
void IntegratorPseudoTransientAdjointSensitivity<Scalar>::
    initializeSolutionHistory(
        Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0,
        Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdot0,
        Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdotdot0,
        Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> y0,
        Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> ydot0,
        Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> ydotdot0)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::createMember;
  using Thyra::VectorSpaceBase;

  //
  // Create and initialize product X, Xdot, Xdotdot

  RCP<const VectorSpaceBase<Scalar>> space = sens_model_->get_x_space();
  RCP<DMVPV> Y                             = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Ydot                          = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Ydotdot                       = rcp_dynamic_cast<DMVPV>(createMember(space));
  const Scalar zero                        = Teuchos::ScalarTraits<Scalar>::zero();

  // y
  if (y0 == Teuchos::null)
    assign(Y->getNonconstMultiVector().ptr(), zero);
  else
    assign(Y->getNonconstMultiVector().ptr(), *y0);

  // ydot
  if (ydot0 == Teuchos::null)
    assign(Ydot->getNonconstMultiVector().ptr(), zero);
  else
    assign(Ydot->getNonconstMultiVector().ptr(), *ydot0);

  // ydotdot
  if (ydotdot0 == Teuchos::null)
    assign(Ydotdot->getNonconstMultiVector().ptr(), zero);
  else
    assign(Ydotdot->getNonconstMultiVector().ptr(), *ydotdot0);

  state_integrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0);
  sens_integrator_->initializeSolutionHistory(t0, Y, Ydot, Ydotdot);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getX() const
{
  return state_integrator_->getX();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getXDot() const
{
  return state_integrator_->getXDot();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getXDotDot() const
{
  return state_integrator_->getXDotDot();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getY() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  RCP<const DMVPV> mvpv =
      rcp_dynamic_cast<const DMVPV>(sens_integrator_->getX());
  return mvpv->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getYDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  RCP<const DMVPV> mvpv =
      rcp_dynamic_cast<const DMVPV>(sens_integrator_->getXDot());
  return mvpv->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getYDotDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  RCP<const DMVPV> mvpv =
      rcp_dynamic_cast<const DMVPV>(sens_integrator_->getXDotDot());
  return mvpv->getMultiVector();
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getG() const
{
  return g_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getDgDp() const
{
  return dgdp_->getMultiVector();
}

template <class Scalar>
std::string IntegratorPseudoTransientAdjointSensitivity<Scalar>::description()
    const
{
  std::string name = "Tempus::IntegratorPseudoTransientAdjointSensitivity";
  return (name);
}

template <class Scalar>
void IntegratorPseudoTransientAdjointSensitivity<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << description() << "::describe" << std::endl;
  state_integrator_->describe(*l_out, verbLevel);
  sens_integrator_->describe(*l_out, verbLevel);
}

template <class Scalar>
void IntegratorPseudoTransientAdjointSensitivity<Scalar>::setParameterList(
    const Teuchos::RCP<Teuchos::ParameterList>& inputPL)
{
  // IntegratorBasic is no longer a Teuchos::ParameterListAcceptor.
  // Since setting the ParameterList is essentially a complete reset,
  // we will rebuild from scratch and reuse the ModelEvaluator.
  auto model = Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar>>(
      state_integrator_->getStepper()->getModel());
  auto tmp_state_integrator = createIntegratorBasic<Scalar>(inputPL, model);
  state_integrator_->copy(tmp_state_integrator);

  model = Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar>>(
      sens_integrator_->getStepper()->getModel());
  auto tmp_sens_integrator = createIntegratorBasic<Scalar>(inputPL, model);
  sens_integrator_->copy(tmp_sens_integrator);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::unsetParameterList()
{
  // IntegratorBasic is no longer a Teuchos::ParameterListAcceptor.
  // We will treat unsetting the ParameterList as a reset to default
  // settings, and reuse the ModelEvaluator.
  auto tmp_state_integrator = createIntegratorBasic<Scalar>();
  auto model                = state_integrator_->getStepper()->getModel();
  tmp_state_integrator->setModel(model);
  state_integrator_->copy(tmp_state_integrator);

  auto tmp_sens_integrator = createIntegratorBasic<Scalar>();
  model                    = sens_integrator_->getStepper()->getModel();
  tmp_sens_integrator->setModel(model);
  sens_integrator_->copy(tmp_sens_integrator);

  auto pl = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      sens_integrator_->getValidParameters());
  return pl;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::RCP<const Teuchos::ParameterList> integrator_pl =
      state_integrator_->getValidParameters();
  Teuchos::RCP<const Teuchos::ParameterList> sensitivity_pl =
      AdjointSensitivityModelEvaluator<Scalar>::getValidParameters();
  pl->setParameters(*integrator_pl);
  pl->sublist("Sensitivities").setParameters(*sensitivity_pl);

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getNonconstParameterList()
{
  auto pl = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      state_integrator_->getValidParameters());
  return pl;
}

template <class Scalar>
SensitivityStepMode
IntegratorPseudoTransientAdjointSensitivity<Scalar>::getStepMode() const
{
  return stepMode_;
}

template <class Scalar>
Teuchos::RCP<AdjointSensitivityModelEvaluator<Scalar>>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::createSensitivityModel(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_solve_model,
    const Teuchos::RCP<Teuchos::ParameterList>& inputPL)
{
  using Teuchos::rcp;

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  if (inputPL != Teuchos::null) {
    *pl = inputPL->sublist("Sensitivities");
  }
  const Scalar tinit  = state_integrator_->getTimeStepControl()->getInitTime();
  const Scalar tfinal = state_integrator_->getTimeStepControl()->getFinalTime();
  return rcp(new AdjointSensitivityModelEvaluator<Scalar>(
      model, adjoint_residual_model, adjoint_solve_model_, tinit, tfinal, true,
      pl));
}

template <class Scalar>
void IntegratorPseudoTransientAdjointSensitivity<Scalar>::buildSolutionHistory()
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::defaultProductVector;
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;
  typedef Thyra::DefaultProductVector<Scalar> DPV;
  typedef Thyra::DefaultProductVectorSpace<Scalar> DPVS;

  // Create combined solution histories, first for the states with zero
  // adjoint and then for the adjoint with frozen states
  auto shPL = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      state_integrator_->getSolutionHistory()->getValidParameters());
  solutionHistory_ = createSolutionHistoryPL<Scalar>(shPL);

  RCP<const VectorSpaceBase<Scalar>> x_space       = model_->get_x_space();
  RCP<const VectorSpaceBase<Scalar>> adjoint_space = sens_model_->get_x_space();
  Teuchos::Array<RCP<const VectorSpaceBase<Scalar>>> spaces(2);
  spaces[0]                  = x_space;
  spaces[1]                  = adjoint_space;
  RCP<const DPVS> prod_space = Thyra::productVectorSpace(spaces());
  const Scalar zero          = Teuchos::ScalarTraits<Scalar>::zero();

  RCP<const SolutionHistory<Scalar>> state_solution_history =
      state_integrator_->getSolutionHistory();
  int num_states = state_solution_history->getNumStates();
  for (int i = 0; i < num_states; ++i) {
    RCP<const SolutionState<Scalar>> state = (*state_solution_history)[i];

    // X
    RCP<DPV> x = defaultProductVector(prod_space);
    assign(x->getNonconstVectorBlock(0).ptr(), *(state->getX()));
    assign(x->getNonconstVectorBlock(1).ptr(), zero);
    RCP<VectorBase<Scalar>> x_b = x;

    // X-Dot
    RCP<VectorBase<Scalar>> x_dot_b;
    if (state->getXDot() != Teuchos::null) {
      RCP<DPV> x_dot = defaultProductVector(prod_space);
      assign(x_dot->getNonconstVectorBlock(0).ptr(), *(state->getXDot()));
      assign(x_dot->getNonconstVectorBlock(1).ptr(), zero);
      x_dot_b = x_dot;
    }

    // X-Dot-Dot
    RCP<VectorBase<Scalar>> x_dot_dot_b;
    if (state->getXDotDot() != Teuchos::null) {
      RCP<DPV> x_dot_dot = defaultProductVector(prod_space);
      assign(x_dot_dot->getNonconstVectorBlock(0).ptr(),
             *(state->getXDotDot()));
      assign(x_dot_dot->getNonconstVectorBlock(1).ptr(), zero);
      x_dot_dot_b = x_dot_dot;
    }

    RCP<SolutionState<Scalar>> prod_state = state->clone();
    prod_state->setX(x_b);
    prod_state->setXDot(x_dot_b);
    prod_state->setXDotDot(x_dot_dot_b);
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
    RCP<DPV> x = defaultProductVector(prod_space);
    assign(x->getNonconstVectorBlock(0).ptr(), *frozen_x);
    assign(x->getNonconstVectorBlock(1).ptr(), *(state->getX()));
    RCP<VectorBase<Scalar>> x_b = x;

    // X-Dot
    RCP<VectorBase<Scalar>> x_dot_b;
    if (state->getXDot() != Teuchos::null) {
      RCP<DPV> x_dot = defaultProductVector(prod_space);
      assign(x_dot->getNonconstVectorBlock(0).ptr(), *frozen_x_dot);
      assign(x_dot->getNonconstVectorBlock(1).ptr(), *(state->getXDot()));
      x_dot_b = x_dot;
    }

    // X-Dot-Dot
    RCP<VectorBase<Scalar>> x_dot_dot_b;
    if (state->getXDotDot() != Teuchos::null) {
      RCP<DPV> x_dot_dot = defaultProductVector(prod_space);
      assign(x_dot_dot->getNonconstVectorBlock(0).ptr(), *frozen_x_dot_dot);
      assign(x_dot_dot->getNonconstVectorBlock(1).ptr(),
             *(state->getXDotDot()));
      x_dot_dot_b = x_dot_dot;
    }

    RCP<SolutionState<Scalar>> prod_state = state->clone();
    prod_state->setX(x_b);
    prod_state->setXDot(x_dot_b);
    prod_state->setXDotDot(x_dot_dot_b);
    solutionHistory_->addState(prod_state, false);
  }
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>>
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model)
{
  Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorPseudoTransientAdjointSensitivity<Scalar>(
          pList, model));
  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>>
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    std::string stepperType)
{
  Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorPseudoTransientAdjointSensitivity<Scalar>(
          model, stepperType));
  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>>
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_model)
{
  Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorPseudoTransientAdjointSensitivity<Scalar>(
          pList, model, adjoint_model));
  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>>
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_model,
    std::string stepperType)
{
  Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorPseudoTransientAdjointSensitivity<Scalar>(
          model, adjoint_model, stepperType));
  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>>
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_solve_model)
{
  Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorPseudoTransientAdjointSensitivity<Scalar>(
          pList, model, adjoint_residual_model, adjoint_solve_model));
  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>>
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& adjoint_solve_model,
    std::string stepperType)
{
  Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorPseudoTransientAdjointSensitivity<Scalar>(
          model, adjoint_residual_model, adjoint_solve_model, stepperType));
  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>>
integratorPseudoTransientAdjointSensitivity()
{
  Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorPseudoTransientAdjointSensitivity<Scalar>());
  return (integrator);
}

}  // namespace Tempus
#endif  // Tempus_IntegratorPseudoTransientAdjointSensitivity_impl_hpp
