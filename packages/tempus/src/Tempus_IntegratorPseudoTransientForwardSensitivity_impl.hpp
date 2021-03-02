// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorPseudoTransientForwardSensitivity_impl_hpp
#define Tempus_IntegratorPseudoTransientForwardSensitivity_impl_hpp

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "Tempus_WrapStaggeredFSAModelEvaluator.hpp"


namespace Tempus {

template<class Scalar>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
IntegratorPseudoTransientForwardSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                inputPL,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model) :
  reuse_solver_(false)
{
  model_ = model;
  sens_model_ = createSensitivityModel(model_, inputPL);
  state_integrator_ = integratorBasic<Scalar>(inputPL, model_);
  sens_integrator_ = integratorBasic<Scalar>(inputPL, sens_model_);
}

template<class Scalar>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
IntegratorPseudoTransientForwardSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType) :
  reuse_solver_(false),
  force_W_update_(false)
{
  model_ = model;
  sens_model_ = createSensitivityModel(model, Teuchos::null);
  state_integrator_ = integratorBasic<Scalar>(model_, stepperType);
  sens_integrator_ = integratorBasic<Scalar>(sens_model_, stepperType);
}

template<class Scalar>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
IntegratorPseudoTransientForwardSensitivity() :
  reuse_solver_(false),
  force_W_update_(false)
{
  state_integrator_ = integratorBasic<Scalar>();
  sens_integrator_ = integratorBasic<Scalar>();
}

template<class Scalar>
bool
IntegratorPseudoTransientForwardSensitivity<Scalar>::
advanceTime()
{
  using Teuchos::RCP;
  using Thyra::VectorBase;

  // Run state integrator and get solution
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
  bool sens_status = sens_integrator_->advanceTime();

  buildSolutionHistory();

  return state_status && sens_status;
}

template<class Scalar>
bool
IntegratorPseudoTransientForwardSensitivity<Scalar>::
advanceTime(const Scalar timeFinal)
{
  using Teuchos::RCP;
  using Thyra::VectorBase;

  // Run state integrator and get solution
  bool state_status = state_integrator_->advanceTime(timeFinal);

  // Set solution in sensitivity ME
  sens_model_->setForwardSolutionState(state_integrator_->getCurrentState());

  // Reuse state solver if requested
  if (reuse_solver_ &&
      state_integrator_->getStepper()->getSolver() != Teuchos::null) {
    sens_model_->setSolver(state_integrator_->getStepper()->getSolver(),
                           force_W_update_);
  }

  // Run sensitivity integrator
  bool sens_status = sens_integrator_->advanceTime(timeFinal);

  buildSolutionHistory();

  return state_status && sens_status;
}

template<class Scalar>
Scalar
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getTime() const
{
  return solutionHistory_->getCurrentTime();
}

template<class Scalar>
int
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getIndex() const
{
  return solutionHistory_->getCurrentIndex();
}

template<class Scalar>
Status
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getStatus() const
{
  Status state_status = state_integrator_->getStatus();
  Status sens_status = sens_integrator_->getStatus();
  if (state_status == FAILED || sens_status == FAILED)
    return FAILED;
  if (state_status == WORKING || sens_status == WORKING)
    return WORKING;
  return PASSED;
}

template<class Scalar>
Teuchos::RCP<Stepper<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getStepper() const
{
  return state_integrator_->getStepper();
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getTempusParameterList()
{
  return state_integrator_->getTempusParameterList();
}

template<class Scalar>
void
IntegratorPseudoTransientForwardSensitivity<Scalar>::
setTempusParameterList(Teuchos::RCP<Teuchos::ParameterList> pl)
{
  state_integrator_->setTempusParameterList(pl);
  sens_integrator_->setTempusParameterList(pl);
}

template<class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getSolutionHistory() const
{
  return solutionHistory_;
}

template<class Scalar>
Teuchos::RCP<const TimeStepControl<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getTimeStepControl() const
{
  return state_integrator_->getTimeStepControl();
}

template<class Scalar>
Teuchos::RCP<TimeStepControl<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getNonConstTimeStepControl()
{
  return state_integrator_->getNonConstTimeStepControl();
}

template<class Scalar>
void IntegratorPseudoTransientForwardSensitivity<Scalar>::
initializeSolutionHistory(Scalar t0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxDp0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotDp0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotdotDp0)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorSpaceBase;
  using Thyra::assign;
  using Thyra::createMember;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  //
  // Create and initialize product X, Xdot, Xdotdot

  RCP< const VectorSpaceBase<Scalar> > space = sens_model_->get_x_space();
  RCP<DMVPV> X       = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdot    = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdotdot = rcp_dynamic_cast<DMVPV>(createMember(space));
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

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

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getX() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> X =
    rcp_dynamic_cast<const DMVPV>(solutionHistory_->getCurrentState()->getX());
  return X->getMultiVector()->col(0);
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getDxDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> X =
    rcp_dynamic_cast<const DMVPV>(solutionHistory_->getCurrentState()->getX());
  const int num_param = X->getMultiVector()->domain()->dim()-1;
  const Teuchos::Range1D rng(1,num_param);
  return X->getMultiVector()->subView(rng);
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getXDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdot =
    rcp_dynamic_cast<const DMVPV>(solutionHistory_->getCurrentState()->getXDot());
  return Xdot->getMultiVector()->col(0);
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getDXDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdot =
    rcp_dynamic_cast<const DMVPV>(solutionHistory_->getCurrentState()->getXDot());
  const int num_param = Xdot->getMultiVector()->domain()->dim()-1;
  const Teuchos::Range1D rng(1,num_param);
  return Xdot->getMultiVector()->subView(rng);
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getXDotDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdotdot =
    rcp_dynamic_cast<const DMVPV>(solutionHistory_->getCurrentState()->getXDotDot());
  return Xdotdot->getMultiVector()->col(0);
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getDXDotDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdotdot =
    rcp_dynamic_cast<const DMVPV>(solutionHistory_->getCurrentState()->getXDotDot());
  const int num_param = Xdotdot->getMultiVector()->domain()->dim()-1;
  const Teuchos::Range1D rng(1,num_param);
  return Xdotdot->getMultiVector()->subView(rng);
}

template<class Scalar>
std::string
IntegratorPseudoTransientForwardSensitivity<Scalar>::
description() const
{
  std::string name = "Tempus::IntegratorPseudoTransientForwardSensitivity";
  return(name);
}

template<class Scalar>
void
IntegratorPseudoTransientForwardSensitivity<Scalar>::
describe(
  Teuchos::FancyOStream          &in_out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  auto out = Teuchos::fancyOStream( in_out.getOStream() );
  out->setOutputToRootOnly(0);
  *out << description() << "::describe" << std::endl;
  state_integrator_->describe(in_out, verbLevel);
  sens_integrator_->describe(in_out, verbLevel);
}

template<class Scalar>
void
IntegratorPseudoTransientForwardSensitivity<Scalar>::
setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & inputPL)
{
  state_integrator_->setParameterList(inputPL);
  sens_integrator_->setParameterList(inputPL);
  reuse_solver_ =
    inputPL->sublist("Sensitivities").get("Reuse State Linear Solver", false);
  force_W_update_ =
    inputPL->sublist("Sensitivities").get("Force W Update", false);
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
unsetParameterList()
{
  state_integrator_->unsetParameterList();
  return sens_integrator_->unsetParameterList();
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::RCP<const Teuchos::ParameterList> integrator_pl =
    state_integrator_->getValidParameters();
  Teuchos::RCP<const Teuchos::ParameterList> sensitivity_pl =
    StaggeredForwardSensitivityModelEvaluator<Scalar>::getValidParameters();
  pl->setParameters(*integrator_pl);
  pl->sublist("Sensitivities").setParameters(*sensitivity_pl);
  pl->sublist("Sensitivities").set("Reuse State Linear Solver", false);
  pl->sublist("Sensitivities").set("Force W Update", false);

  return pl;
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientForwardSensitivity<Scalar>::
getNonconstParameterList()
{
  return state_integrator_->getNonconstParameterList();
}

template <class Scalar>
Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> >
IntegratorPseudoTransientForwardSensitivity<Scalar>::
createSensitivityModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const Teuchos::RCP<Teuchos::ParameterList>& inputPL)
{
  using Teuchos::rcp;

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  if (inputPL != Teuchos::null) {
    *pl = inputPL->sublist("Sensitivities");
  }
  reuse_solver_ = pl->get("Reuse State Linear Solver", false);
  force_W_update_ = pl->get("Force W Update", true);
  pl->remove("Reuse State Linear Solver");
  pl->remove("Force W Update");
  return wrapStaggeredFSAModelEvaluator(model, pl);
}

template<class Scalar>
void
IntegratorPseudoTransientForwardSensitivity<Scalar>::
buildSolutionHistory()
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
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  // Create combined solution histories, first for the states with zero
  // sensitivities and then for the sensitivities with frozen states
  RCP<ParameterList> shPL =
    Teuchos::sublist(state_integrator_->getIntegratorParameterList(),
                     "Solution History", true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

  const int num_param =
    rcp_dynamic_cast<const DMVPV>(sens_integrator_->getX())->getMultiVector()->domain()->dim();
  RCP<const VectorSpaceBase<Scalar> > x_space = model_->get_x_space();
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > prod_space =
    Thyra::multiVectorProductVectorSpace(x_space, num_param+1);
  const Teuchos::Range1D rng(1,num_param);
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  RCP<const SolutionHistory<Scalar> > state_solution_history =
    state_integrator_->getSolutionHistory();
  int num_states = state_solution_history->getNumStates();
  for (int i=0; i<num_states; ++i) {
    RCP<const SolutionState<Scalar> > state = (*state_solution_history)[i];

    // X
    RCP< MultiVectorBase<Scalar> > x_mv =
      createMembers(x_space, num_param+1);
    assign(x_mv->col(0).ptr(), *(state->getX()));
    assign(x_mv->subView(rng).ptr(), zero);
    RCP<VectorBase<Scalar> > x = multiVectorProductVector(prod_space, x_mv);

    // X-Dot
    RCP<VectorBase<Scalar> > x_dot;
    if (state->getXDot() != Teuchos::null) {
      RCP< MultiVectorBase<Scalar> > x_dot_mv =
        createMembers(x_space, num_param+1);
      assign(x_dot_mv->col(0).ptr(), *(state->getXDot()));
      assign(x_dot_mv->subView(rng).ptr(), zero);
      x_dot = multiVectorProductVector(prod_space, x_dot_mv);
    }

    // X-Dot-Dot
    RCP<VectorBase<Scalar> > x_dot_dot;
    if (state->getXDotDot() != Teuchos::null) {
      RCP< MultiVectorBase<Scalar> > x_dot_dot_mv =
        createMembers(x_space, num_param+1);
      assign(x_dot_dot_mv->col(0).ptr(), *(state->getXDotDot()));
      assign(x_dot_dot_mv->subView(rng).ptr(), zero);
      x_dot_dot = multiVectorProductVector(prod_space, x_dot_dot_mv);
    }

    RCP<SolutionState<Scalar> > prod_state = state->clone();
    prod_state->setX(x);
    prod_state->setXDot(x_dot);
    prod_state->setXDotDot(x_dot_dot);
    solutionHistory_->addState(prod_state);
  }

  RCP<const VectorBase<Scalar> > frozen_x =
    state_solution_history->getCurrentState()->getX();
  RCP<const VectorBase<Scalar> > frozen_x_dot =
    state_solution_history->getCurrentState()->getXDot();
  RCP<const VectorBase<Scalar> > frozen_x_dot_dot =
    state_solution_history->getCurrentState()->getXDotDot();
  RCP<const SolutionHistory<Scalar> > sens_solution_history =
    sens_integrator_->getSolutionHistory();
  num_states = sens_solution_history->getNumStates();
  for (int i=0; i<num_states; ++i) {
    RCP<const SolutionState<Scalar> > state = (*sens_solution_history)[i];

    // X
    RCP< MultiVectorBase<Scalar> > x_mv =
      createMembers(x_space, num_param+1);
    RCP<const MultiVectorBase<Scalar> > dxdp =
      rcp_dynamic_cast<const DMVPV>(state->getX())->getMultiVector();
    assign(x_mv->col(0).ptr(), *(frozen_x));
    assign(x_mv->subView(rng).ptr(), *dxdp);
    RCP<VectorBase<Scalar> > x = multiVectorProductVector(prod_space, x_mv);

    // X-Dot
    RCP<VectorBase<Scalar> > x_dot;
    if (state->getXDot() != Teuchos::null) {
      RCP< MultiVectorBase<Scalar> > x_dot_mv =
        createMembers(x_space, num_param+1);
      RCP<const MultiVectorBase<Scalar> > dxdotdp =
        rcp_dynamic_cast<const DMVPV>(state->getXDot())->getMultiVector();
      assign(x_dot_mv->col(0).ptr(), *(frozen_x_dot));
      assign(x_dot_mv->subView(rng).ptr(), *dxdotdp);
      x_dot = multiVectorProductVector(prod_space, x_dot_mv);
    }

    // X-Dot-Dot
    RCP<VectorBase<Scalar> > x_dot_dot;
    if (state->getXDotDot() != Teuchos::null) {
      RCP< MultiVectorBase<Scalar> > x_dot_dot_mv =
        createMembers(x_space, num_param+1);
      RCP<const MultiVectorBase<Scalar> > dxdotdotdp =
        rcp_dynamic_cast<const DMVPV>(state->getXDotDot())->getMultiVector();
      assign(x_dot_dot_mv->col(0).ptr(), *(frozen_x_dot_dot));
      assign(x_dot_dot_mv->subView(rng).ptr(), *dxdotdotdp);
      x_dot_dot = multiVectorProductVector(prod_space, x_dot_dot_mv);
    }

    RCP<SolutionState<Scalar> > prod_state = state->clone();
    prod_state->setX(x);
    prod_state->setXDot(x_dot);
    prod_state->setXDotDot(x_dot_dot);
    solutionHistory_->addState(prod_state);
  }
}

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar> >
integratorPseudoTransientForwardSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model)
{
  Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>(pList, model));
  return(integrator);
}

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar> >
integratorPseudoTransientForwardSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model,
  std::string stepperType)
{
  Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>(model, stepperType));
  return(integrator);
}

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar> >
integratorPseudoTransientForwardSensitivity()
{
  Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>());
  return(integrator);
}

} // namespace Tempus
#endif // Tempus_IntegratorPseudoTransientForwardSensitivity_impl_hpp
