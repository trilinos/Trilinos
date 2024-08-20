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

#include "Teuchos_ParameterList.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_ImplicitAdjointModelEvaluator.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

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
    const bool mass_matrix_is_identity)
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
    stepMode_(SensitivityStepMode::Forward)
{
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

  // For at least some time-stepping methods, the time of the last time step
  // may not be timeFinal (e.g., it may be greater by at most delta_t).
  // But since the adjoint model requires timeFinal in its formulation, reset
  // it to the achieved final time.
  adjoint_aux_model_->setFinalTime(state_integrator_->getTime());

  // Set solution history in adjoint stepper
  adjoint_aux_model_->setForwardSolutionHistory(state_solution_history);

  // Compute dg/dx
  RCP<const VectorSpaceBase<Scalar>> g_space = model_->get_g_space(g_index_);
  RCP<const VectorSpaceBase<Scalar>> x_space = model_->get_x_space();
  const int num_g                            = g_space->dim();
  RCP<MultiVectorBase<Scalar>> dgdx          = createMembers(x_space, num_g);
  MEB::InArgs<Scalar> inargs                 = model_->getNominalValues();
  RCP<const SolutionState<Scalar>> state =
      state_solution_history->getCurrentState();
  inargs.set_t(state->getTime());
  inargs.set_x(state->getX());
  inargs.set_x_dot(state->getXDot());
  MEB::OutArgs<Scalar> outargs     = model_->createOutArgs();
  MEB::OutArgs<Scalar> adj_outargs = adjoint_model_->createOutArgs();
  outargs.set_DgDx(g_index_,
                   MEB::Derivative<Scalar>(dgdx, MEB::DERIV_MV_GRADIENT_FORM));
  model_->evalModel(inargs, outargs);
  outargs.set_DgDx(g_index_, MEB::Derivative<Scalar>());

  // Compute ICs == [ (df/dx_dot)^{-T} (dg/dx)^T; 0 ]
  // For explicit form, we are relying on the user to inform us the
  // the mass matrix is the identity.  It would be nice to be able to determine
  // somehow automatically that we are using an explicit stepper.
  RCP<DPV> adjoint_init = rcp_dynamic_cast<DPV>(
      Thyra::createMember(adjoint_aux_model_->get_x_space()));
  RCP<MultiVectorBase<Scalar>> adjoint_init_mv =
      rcp_dynamic_cast<DMVPV>(adjoint_init->getNonconstVectorBlock(0))
          ->getNonconstMultiVector();
  assign(adjoint_init->getNonconstVectorBlock(1).ptr(),
         Teuchos::ScalarTraits<Scalar>::zero());
  if (mass_matrix_is_identity_)
    assign(adjoint_init_mv.ptr(), *dgdx);
  else {
    inargs.set_alpha(1.0);
    inargs.set_beta(0.0);
    RCP<LinearOpWithSolveBase<Scalar>> W;
    if (adj_outargs.supports(MEB::OUT_ARG_W)) {
      // Model supports W
      W = adjoint_model_->create_W();
      adj_outargs.set_W(W);
      adjoint_model_->evalModel(inargs, adj_outargs);
      adj_outargs.set_W(Teuchos::null);
    }
    else {
      // Otherwise model must support a W_op and W factory
      RCP<const LinearOpWithSolveFactoryBase<Scalar>> lowsfb =
          adjoint_model_->get_W_factory();
      TEUCHOS_TEST_FOR_EXCEPTION(lowsfb == Teuchos::null, std::logic_error,
                                 "Adjoint ME must support W out-arg or provide "
                                 "a W_factory for non-identity mass matrix");

      // Compute W_op (and W_prec if supported)
      RCP<LinearOpBase<Scalar>> W_op = adjoint_model_->create_W_op();
      adj_outargs.set_W_op(W_op);
      RCP<PreconditionerFactoryBase<Scalar>> prec_factory =
          lowsfb->getPreconditionerFactory();
      RCP<PreconditionerBase<Scalar>> W_prec;
      if (prec_factory != Teuchos::null)
        W_prec = prec_factory->createPrec();
      else if (adj_outargs.supports(MEB::OUT_ARG_W_prec)) {
        W_prec = adjoint_model_->create_W_prec();
        adj_outargs.set_W_prec(W_prec);
      }
      adjoint_model_->evalModel(inargs, adj_outargs);
      adj_outargs.set_W_op(Teuchos::null);
      if (adj_outargs.supports(MEB::OUT_ARG_W_prec))
        adj_outargs.set_W_prec(Teuchos::null);

      // Create and initialize W
      W = lowsfb->createOp();
      if (W_prec != Teuchos::null) {
        if (prec_factory != Teuchos::null)
          prec_factory->initializePrec(
              Thyra::defaultLinearOpSource<Scalar>(W_op), W_prec.get());
        Thyra::initializePreconditionedOp<Scalar>(*lowsfb, W_op, W_prec,
                                                  W.ptr());
      }
      else
        Thyra::initializeOp<Scalar>(*lowsfb, W_op, W.ptr());
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
        W == Teuchos::null, std::logic_error,
        "A null W has been encountered in "
        "Tempus::IntegratorAdjointSensitivity::advanceTime!\n");
    // Initialize adjoint_init_mv to zero before solve for linear solvers that
    // use what is passed in as the initial guess
    assign(adjoint_init_mv.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
    W->solve(Thyra::NOTRANS, *dgdx, adjoint_init_mv.ptr());
  }

  // Run sensitivity integrator and get solution
  bool sens_status = true;
  {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::IntegratorAdjointSensitivity::advanceTime::adjoint",
        TEMPUS_AS_AT_ADJ);
    stepMode_ = SensitivityStepMode::Adjoint;
    const Scalar tinit =
        adjoint_integrator_->getTimeStepControl()->getInitTime();
    adjoint_integrator_->initializeSolutionHistory(tinit, adjoint_init);
    sens_status = adjoint_integrator_->advanceTime(timeFinal);
  }
  RCP<const SolutionHistory<Scalar>> adjoint_solution_history =
      adjoint_integrator_->getSolutionHistory();

  // Compute dg/dp at final time T
  RCP<const VectorSpaceBase<Scalar>> p_space = model_->get_p_space(p_index_);
  dgdp_                                      = createMembers(p_space, num_g);
  if (g_depends_on_p_) {
    MEB::DerivativeSupport dgdp_support =
        outargs.supports(MEB::OUT_ARG_DgDp, g_index_, p_index_);
    if (dgdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
      outargs.set_DgDp(
          g_index_, p_index_,
          MEB::Derivative<Scalar>(dgdp_, MEB::DERIV_MV_GRADIENT_FORM));
      model_->evalModel(inargs, outargs);
    }
    else if (dgdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
      const int num_p                         = p_space->dim();
      RCP<MultiVectorBase<Scalar>> dgdp_trans = createMembers(g_space, num_p);
      outargs.set_DgDp(
          g_index_, p_index_,
          MEB::Derivative<Scalar>(dgdp_trans, MEB::DERIV_MV_JACOBIAN_FORM));
      model_->evalModel(inargs, outargs);
      Thyra::DetachedMultiVectorView<Scalar> dgdp_view(*dgdp_);
      Thyra::DetachedMultiVectorView<Scalar> dgdp_trans_view(*dgdp_trans);
      for (int i = 0; i < num_p; ++i)
        for (int j = 0; j < num_g; ++j) dgdp_view(i, j) = dgdp_trans_view(j, i);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "Invalid dg/dp support");
    outargs.set_DgDp(g_index_, p_index_, MEB::Derivative<Scalar>());
  }
  else
    assign(dgdp_.ptr(), Scalar(0.0));

  // Add in initial condition term = (dx/dp^T(0))*(df/dx_dot^T(0))*y(0)
  // If dxdp_init_ is null, assume it is zero
  if (ic_depends_on_p_ && dxdp_init_ != Teuchos::null) {
    RCP<const SolutionState<Scalar>> adjoint_state =
        adjoint_solution_history->getCurrentState();
    RCP<const VectorBase<Scalar>> adjoint_x =
        rcp_dynamic_cast<const DPV>(adjoint_state->getX())->getVectorBlock(0);
    RCP<const MultiVectorBase<Scalar>> adjoint_mv =
        rcp_dynamic_cast<const DMVPV>(adjoint_x)->getMultiVector();
    if (mass_matrix_is_identity_)
      dxdp_init_->apply(Thyra::CONJTRANS, *adjoint_mv, dgdp_.ptr(), Scalar(1.0),
                        Scalar(1.0));
    else {
      inargs.set_t(initial_state->getTime());
      inargs.set_x(initial_state->getX());
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
  // (2nd block of the solution vector)
  if (f_depends_on_p_) {
    RCP<const SolutionState<Scalar>> adjoint_state =
        adjoint_solution_history->getCurrentState();
    RCP<const VectorBase<Scalar>> z =
        rcp_dynamic_cast<const DPV>(adjoint_state->getX())->getVectorBlock(1);
    RCP<const MultiVectorBase<Scalar>> z_mv =
        rcp_dynamic_cast<const DMVPV>(z)->getMultiVector();
    Thyra::V_VmV(dgdp_.ptr(), *dgdp_, *z_mv);
  }

  buildSolutionHistory(state_solution_history, adjoint_solution_history);

  return state_status && sens_status;
}

template <class Scalar>
Scalar IntegratorAdjointSensitivity<Scalar>::getTime() const
{
  return solutionHistory_->getCurrentTime();
}

template <class Scalar>
int IntegratorAdjointSensitivity<Scalar>::getIndex() const
{
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
  if (spl->isParameter("Response Depends on Parameters"))
    spl->remove("Response Depends on Parameters");
  if (spl->isParameter("Residual Depends on Parameters"))
    spl->remove("Residual Depends on Parameters");
  if (spl->isParameter("IC Depends on Parameters"))
    spl->remove("IC Depends on Parameters");

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

  int num_states       = state_solution_history->getNumStates();
  const Scalar t_init  = state_integrator_->getTimeStepControl()->getInitTime();
  const Scalar t_final = state_integrator_->getTime();
  for (int i = 0; i < num_states; ++i) {
    RCP<const SolutionState<Scalar>> forward_state =
        (*state_solution_history)[i];
    RCP<const SolutionState<Scalar>> adjoint_state =
        adjoint_solution_history->findState(t_final + t_init -
                                            forward_state->getTime());

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

  int p_index          = spl->get<int>("Sensitivity Parameter Index", 0);
  int g_index          = spl->get<int>("Response Function Index", 0);
  bool g_depends_on_p  = spl->get<bool>("Response Depends on Parameters", true);
  bool f_depends_on_p  = spl->get<bool>("Residual Depends on Parameters", true);
  bool ic_depends_on_p = spl->get<bool>("IC Depends on Parameters", true);
  bool mass_matrix_is_identity =
      spl->get<bool>("Mass Matrix Is Identity", false);

  auto state_integrator = createIntegratorBasic<Scalar>(inputPL, model);

  // createAdjointModel
  if (spl->isParameter("Response Depends on Parameters"))
    spl->remove("Response Depends on Parameters");
  if (spl->isParameter("Residual Depends on Parameters"))
    spl->remove("Residual Depends on Parameters");
  if (spl->isParameter("IC Depends on Parameters"))
    spl->remove("IC Depends on Parameters");

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

  auto adjoint_integrator =
      createIntegratorBasic<Scalar>(inputPL, adjoint_aux_model);

  Teuchos::RCP<IntegratorAdjointSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorAdjointSensitivity<Scalar>(
          model, state_integrator, adjt_model, adjoint_aux_model,
          adjoint_integrator, combined_solution_History, p_index, g_index,
          g_depends_on_p, f_depends_on_p, ic_depends_on_p,
          mass_matrix_is_identity));

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
