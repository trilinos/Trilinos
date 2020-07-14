// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorAdjointSensitivity_impl_hpp
#define Tempus_IntegratorAdjointSensitivity_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "NOX_Thyra.H"

namespace Tempus {

template<class Scalar>
IntegratorAdjointSensitivity<Scalar>::
IntegratorAdjointSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                inputPL,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model)
{
  model_ = model;
  setParameterList(inputPL);
  state_integrator_ = integratorBasic<Scalar>(inputPL, model_);

  TEUCHOS_TEST_FOR_EXCEPTION( getStepper()->getUseFSAL(), std::logic_error,
    "Error - IntegratorAdjointSensitivity(): Cannot use FSAL with\n"
    "        IntegratorAdjointSensitivity, because the state and adjoint\n"
    "        integrators require ModelEvaluator evaluation in the\n"
    "        constructor to make the initial conditions consistent.\n"
    "        For the adjoint integrator, this requires special construction\n"
    "        which has not been implemented yet.\n");

  adjoint_model_ = createAdjointModel(model_, inputPL);
  adjoint_integrator_ = integratorBasic<Scalar>(inputPL, adjoint_model_);
}

template<class Scalar>
IntegratorAdjointSensitivity<Scalar>::
IntegratorAdjointSensitivity()
{
  state_integrator_ = integratorBasic<Scalar>();
  adjoint_integrator_ = integratorBasic<Scalar>();
}

template<class Scalar>
bool
IntegratorAdjointSensitivity<Scalar>::
advanceTime()
{
  const Scalar tfinal =
    state_integrator_->getTimeStepControl()->getFinalTime();
  return advanceTime(tfinal);
}

template<class Scalar>
bool
IntegratorAdjointSensitivity<Scalar>::
advanceTime(const Scalar timeFinal)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;
  using Thyra::MultiVectorBase;
  using Thyra::LinearOpBase;
  using Thyra::LinearOpWithSolveBase;
  using Thyra::createMember;
  using Thyra::createMembers;
  using Thyra::assign;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  typedef Thyra::DefaultProductVector<Scalar> DPV;

  // Get initial state for later
  RCP<const SolutionHistory<Scalar> > state_solution_history =
    state_integrator_->getSolutionHistory();
  RCP<const SolutionState<Scalar> > initial_state =
    (*state_solution_history)[0];

  // Run state integrator and get solution
  bool state_status = state_integrator_->advanceTime(timeFinal);

  // Set solution history in adjoint stepper
  adjoint_model_->setForwardSolutionHistory(state_solution_history);

  // Compute dg/dx
  RCP<const VectorSpaceBase<Scalar> > g_space = model_->get_g_space(g_index_);
  RCP<const VectorSpaceBase<Scalar> > x_space = model_->get_x_space();
  const int num_g = g_space->dim();
  RCP<MultiVectorBase<Scalar> > dgdx = createMembers(x_space, num_g);
  MEB::InArgs<Scalar> inargs = model_->getNominalValues();
  RCP<const SolutionState<Scalar> > state =
    state_solution_history->getCurrentState();
  inargs.set_t(state->getTime());
  inargs.set_x(state->getX());
  inargs.set_x_dot(state->getXDot());
  MEB::OutArgs<Scalar> outargs = model_->createOutArgs();
  outargs.set_DgDx(g_index_,
                   MEB::Derivative<Scalar>(dgdx, MEB::DERIV_MV_GRADIENT_FORM));
  model_->evalModel(inargs, outargs);
  outargs.set_DgDx(g_index_, MEB::Derivative<Scalar>());

  // Compute ICs == [ (df/dx_dot)^{-T} (dg/dx)^T; 0 ]
  // For explicit form, we are relying on the user to inform us the
  // the mass matrix is the identity.  It would be nice to be able to determine
  // somehow automatically that we are using an explicit stepper.
  RCP<DPV> adjoint_init =
    rcp_dynamic_cast<DPV>(Thyra::createMember(adjoint_model_->get_x_space()));
  RCP<MultiVectorBase<Scalar> > adjoint_init_mv =
    rcp_dynamic_cast<DMVPV>(adjoint_init->getNonconstVectorBlock(0))->getNonconstMultiVector();
  assign(adjoint_init->getNonconstVectorBlock(1).ptr(),
         Teuchos::ScalarTraits<Scalar>::zero());
  if (mass_matrix_is_identity_)
    assign(adjoint_init_mv.ptr(), *dgdx);
  else {
    inargs.set_alpha(1.0);
    inargs.set_beta(0.0);
    RCP<LinearOpWithSolveBase<Scalar> > W = model_->create_W();
    outargs.set_W(W);
    model_->evalModel(inargs, outargs);
    W->solve(Thyra::CONJTRANS, *dgdx, adjoint_init_mv.ptr());
    outargs.set_W(Teuchos::null);
  }

  // Run sensitivity integrator and get solution
  adjoint_integrator_->initializeSolutionHistory(Scalar(0.0), adjoint_init);
  bool sens_status = adjoint_integrator_->advanceTime(timeFinal);
  RCP<const SolutionHistory<Scalar> > adjoint_solution_history =
    adjoint_integrator_->getSolutionHistory();

  // Compute dg/dp at final time T
  RCP<const VectorSpaceBase<Scalar> > p_space = model_->get_p_space(p_index_);
  dgdp_ = createMembers(p_space, num_g);
  if (g_depends_on_p_) {
    MEB::DerivativeSupport dgdp_support =
      outargs.supports(MEB::OUT_ARG_DgDp, g_index_, p_index_);
    if (dgdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
      outargs.set_DgDp(g_index_, p_index_,
                       MEB::Derivative<Scalar>(dgdp_,
                                               MEB::DERIV_MV_GRADIENT_FORM));
      model_->evalModel(inargs, outargs);
    }
    else if (dgdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
      const int num_p = p_space->dim();
      RCP<MultiVectorBase<Scalar> > dgdp_trans =
        createMembers(g_space, num_p);
      outargs.set_DgDp(g_index_, p_index_,
                       MEB::Derivative<Scalar>(dgdp_trans,
                                               MEB::DERIV_MV_JACOBIAN_FORM));
      model_->evalModel(inargs, outargs);
      Thyra::DetachedMultiVectorView<Scalar> dgdp_view(*dgdp_);
      Thyra::DetachedMultiVectorView<Scalar> dgdp_trans_view(*dgdp_trans);
      for (int i=0; i<num_p; ++i)
        for (int j=0; j<num_g; ++j)
          dgdp_view(i,j) = dgdp_trans_view(j,i);
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
    RCP<const SolutionState<Scalar> > adjoint_state =
      adjoint_solution_history->getCurrentState();
    RCP<const VectorBase<Scalar> > adjoint_x =
      rcp_dynamic_cast<const DPV>(adjoint_state->getX())->getVectorBlock(0);
    RCP<const MultiVectorBase<Scalar> > adjoint_mv =
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
      RCP<LinearOpBase<Scalar> > W_op = model_->create_W_op();
      outargs.set_W_op(W_op);
      model_->evalModel(inargs, outargs);
      outargs.set_W_op(Teuchos::null);
      RCP<MultiVectorBase<Scalar> > tmp = createMembers(x_space, num_g);
      W_op->apply(Thyra::CONJTRANS, *adjoint_mv, tmp.ptr(), Scalar(1.0),
                  Scalar(0.0));
      dxdp_init_->apply(Thyra::CONJTRANS, *tmp, dgdp_.ptr(), Scalar(1.0),
                        Scalar(1.0));
    }
  }

  // Add in model parameter term = \int_0^T( (df/dp^T(t)*y(t) )dt which
  // is computed during the adjoint integration as an auxiliary integral
  // (2nd block of the solution vector)
  if (f_depends_on_p_) {
    RCP<const SolutionState<Scalar> > adjoint_state =
      adjoint_solution_history->getCurrentState();
    RCP<const VectorBase<Scalar> > z =
      rcp_dynamic_cast<const DPV>(adjoint_state->getX())->getVectorBlock(1);
    RCP<const MultiVectorBase<Scalar> > z_mv =
      rcp_dynamic_cast<const DMVPV>(z)->getMultiVector();
    Thyra::V_VmV(dgdp_.ptr(), *dgdp_, *z_mv);
  }

  buildSolutionHistory(state_solution_history, adjoint_solution_history);

  return state_status && sens_status;
}

template<class Scalar>
Scalar
IntegratorAdjointSensitivity<Scalar>::
getTime() const
{
  return solutionHistory_->getCurrentTime();
}

template<class Scalar>
int
IntegratorAdjointSensitivity<Scalar>::
getIndex() const
{
  return solutionHistory_->getCurrentIndex();
}

template<class Scalar>
Status
IntegratorAdjointSensitivity<Scalar>::
getStatus() const
{
  Status state_status = state_integrator_->getStatus();
  Status sens_status = adjoint_integrator_->getStatus();
  if (state_status == FAILED || sens_status == FAILED)
    return FAILED;
  if (state_status == WORKING || sens_status == WORKING)
    return WORKING;
  return PASSED;
}

template<class Scalar>
Teuchos::RCP<Stepper<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getStepper() const
{
  return state_integrator_->getStepper();
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorAdjointSensitivity<Scalar>::
getTempusParameterList()
{
  return state_integrator_->getTempusParameterList();
}

template<class Scalar>
void
IntegratorAdjointSensitivity<Scalar>::
setTempusParameterList(Teuchos::RCP<Teuchos::ParameterList> pl)
{
  state_integrator_->setTempusParameterList(pl);
  adjoint_integrator_->setTempusParameterList(pl);
}

template<class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getSolutionHistory() const
{
  return solutionHistory_;
}

template<class Scalar>
Teuchos::RCP<const TimeStepControl<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getTimeStepControl() const
{
  return state_integrator_->getTimeStepControl();
}

template<class Scalar>
Teuchos::RCP<TimeStepControl<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template<class Scalar>
void IntegratorAdjointSensitivity<Scalar>::
initializeSolutionHistory(Scalar t0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxDp0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > /* DxdotDp0 */,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > /* DxdotdotDp0 */)
{
  state_integrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0);
  dxdp_init_ = DxDp0;
}

template<class Scalar>
void IntegratorAdjointSensitivity<Scalar>::
setObserver(Teuchos::RCP<IntegratorObserver<Scalar> > obs)
{
  state_integrator_->setObserver(obs);
  // Currently not setting observer on adjoint integrator because it isn't
  // clear what we want to do with it
  //adjoint_integrator_->setObserver(obs);
}

template<class Scalar>
void IntegratorAdjointSensitivity<Scalar>::
initialize()
{
  state_integrator_->initialize();
  adjoint_integrator_->initialize();
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getX() const
{
  return state_integrator_->getX();
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getXdot() const
{
  return state_integrator_->getXdot();
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getXdotdot() const
{
  return state_integrator_->getXdotdot();
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
getDgDp() const
{
  return dgdp_;
}

template<class Scalar>
std::string
IntegratorAdjointSensitivity<Scalar>::
description() const
{
  std::string name = "Tempus::IntegratorAdjointSensitivity";
  return(name);
}

template<class Scalar>
void
IntegratorAdjointSensitivity<Scalar>::
describe(
  Teuchos::FancyOStream          &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  out << description() << "::describe" << std::endl;
  state_integrator_->describe(out, verbLevel);
  adjoint_integrator_->describe(out, verbLevel);
}

template<class Scalar>
void
IntegratorAdjointSensitivity<Scalar>::
setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & inputPL)
{
  if (state_integrator_ != Teuchos::null)
    state_integrator_->setParameterList(inputPL);
  if (adjoint_integrator_ != Teuchos::null)
    adjoint_integrator_->setParameterList(inputPL);
  Teuchos::ParameterList& spl = inputPL->sublist("Sensitivities");
  p_index_ = spl.get<int>("Sensitivity Parameter Index", 0);
  g_index_ = spl.get<int>("Response Function Index", 0);
  g_depends_on_p_ = spl.get<bool>("Response Depends on Parameters", true);
  f_depends_on_p_ = spl.get<bool>("Residual Depends on Parameters", true);
  ic_depends_on_p_ = spl.get<bool>("IC Depends on Parameters", true);
  mass_matrix_is_identity_ = spl.get<bool>("Mass Matrix Is Identity", false);
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorAdjointSensitivity<Scalar>::
unsetParameterList()
{
  state_integrator_->unsetParameterList();
  return adjoint_integrator_->unsetParameterList();
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorAdjointSensitivity<Scalar>::
getValidParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  RCP<ParameterList> pl = parameterList();
  RCP<const ParameterList> integrator_pl =
    state_integrator_->getValidParameters();
  RCP<const ParameterList> sensitivity_pl =
    AdjointAuxSensitivityModelEvaluator<Scalar>::getValidParameters();
  pl->setParameters(*integrator_pl);
  ParameterList& spl = pl->sublist("Sensitivities");
  spl.setParameters(*sensitivity_pl);
  spl.set<bool>("Response Depends on Parameters", true);
  spl.set<bool>("Residual Depends on Parameters", true);
  spl.set<bool>("IC Depends on Parameters", true);

  return pl;
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorAdjointSensitivity<Scalar>::
getNonconstParameterList()
{
  return state_integrator_->getNonconstParameterList();
}

template <class Scalar>
Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar> >
IntegratorAdjointSensitivity<Scalar>::
createAdjointModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const Teuchos::RCP<Teuchos::ParameterList>& inputPL)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  RCP<ParameterList> spl = Teuchos::parameterList();
  if (inputPL != Teuchos::null) {
    *spl = inputPL->sublist("Sensitivities");
  }
  spl->remove("Response Depends on Parameters");
  spl->remove("Residual Depends on Parameters");
  spl->remove("IC Depends on Parameters");
  const Scalar tfinal = state_integrator_->getTimeStepControl()->getFinalTime();
  return rcp(new AdjointAuxSensitivityModelEvaluator<Scalar>(
               model, tfinal, spl));
}

template<class Scalar>
void
IntegratorAdjointSensitivity<Scalar>::
buildSolutionHistory(
  const Teuchos::RCP<const SolutionHistory<Scalar> >& state_solution_history,
  const Teuchos::RCP<const SolutionHistory<Scalar> >& adjoint_solution_history)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ParameterList;
  using Thyra::VectorBase;
  using Thyra::MultiVectorBase;
  using Thyra::VectorSpaceBase;
  using Thyra::createMembers;
  using Thyra::multiVectorProductVector;
  using Thyra::assign;
  typedef Thyra::DefaultProductVectorSpace<Scalar> DPVS;
  typedef Thyra::DefaultProductVector<Scalar> DPV;

  // Create combined solution histories combining the forward and adjoint
  // solutions.  We do not include the auxiliary part from the adjoint solution.
  RCP<ParameterList> shPL =
    Teuchos::sublist(state_integrator_->getIntegratorParameterList(),
                     "Solution History", true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

  RCP<const VectorSpaceBase<Scalar> > x_space = model_->get_x_space();
  RCP<const VectorSpaceBase<Scalar> > adjoint_space =
    rcp_dynamic_cast<const DPVS>(adjoint_model_->get_x_space())->getBlock(0);
  Teuchos::Array< RCP<const VectorSpaceBase<Scalar> > > spaces(2);
  spaces[0] = x_space;
  spaces[1] = adjoint_space;
  RCP<const DPVS > prod_space = Thyra::productVectorSpace(spaces());

  int num_states = state_solution_history->getNumStates();
  const Scalar t_final = state_integrator_->getTime();
  for (int i=0; i<num_states; ++i) {
    RCP<const SolutionState<Scalar> > forward_state =
      (*state_solution_history)[i];
    RCP<const SolutionState<Scalar> > adjoint_state =
      adjoint_solution_history->findState(t_final-forward_state->getTime());

    // X
    RCP<DPV> x = Thyra::defaultProductVector(prod_space);
    RCP<const VectorBase<Scalar> > adjoint_x =
      rcp_dynamic_cast<const DPV>(adjoint_state->getX())->getVectorBlock(0);
    assign(x->getNonconstVectorBlock(0).ptr(), *(forward_state->getX()));
    assign(x->getNonconstVectorBlock(1).ptr(), *(adjoint_x));
    RCP<VectorBase<Scalar> > x_b = x;

    // X-Dot
    RCP<DPV> x_dot = Thyra::defaultProductVector(prod_space);
    RCP<const VectorBase<Scalar> > adjoint_x_dot =
      rcp_dynamic_cast<const DPV>(adjoint_state->getXDot())->getVectorBlock(0);
    assign(x_dot->getNonconstVectorBlock(0).ptr(), *(forward_state->getXDot()));
    assign(x_dot->getNonconstVectorBlock(1).ptr(), *(adjoint_x_dot));
    RCP<VectorBase<Scalar> > x_dot_b = x_dot;

    // X-Dot-Dot
    RCP<DPV> x_dot_dot;
    if (forward_state->getXDotDot() != Teuchos::null) {
       x_dot_dot = Thyra::defaultProductVector(prod_space);
       RCP<const VectorBase<Scalar> > adjoint_x_dot_dot =
         rcp_dynamic_cast<const DPV>(
           adjoint_state->getXDotDot())->getVectorBlock(0);
       assign(x_dot_dot->getNonconstVectorBlock(0).ptr(),
              *(forward_state->getXDotDot()));
       assign(x_dot_dot->getNonconstVectorBlock(1).ptr(),
              *(adjoint_x_dot_dot));
    }
    RCP<VectorBase<Scalar> > x_dot_dot_b = x_dot_dot;

    RCP<SolutionState<Scalar> > prod_state = forward_state->clone();
    prod_state->setX(x_b);
    prod_state->setXDot(x_dot_b);
    prod_state->setXDotDot(x_dot_dot_b);
    prod_state->setPhysicsState(Teuchos::null);
    solutionHistory_->addState(prod_state);
  }
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar> >
integratorAdjointSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model)
{
  Teuchos::RCP<IntegratorAdjointSensitivity<Scalar> > integrator =
    Teuchos::rcp(new IntegratorAdjointSensitivity<Scalar>(pList, model));
  return(integrator);
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar> >
integratorAdjointSensitivity()
{
  Teuchos::RCP<IntegratorAdjointSensitivity<Scalar> > integrator =
    Teuchos::rcp(new IntegratorAdjointSensitivity<Scalar>());
  return(integrator);
}

} // namespace Tempus
#endif // Tempus_IntegratorAdjointSensitivity_impl_hpp
