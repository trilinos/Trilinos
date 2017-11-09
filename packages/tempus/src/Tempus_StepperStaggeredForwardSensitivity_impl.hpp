// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperStaggeredForwardSensitivity_impl_hpp
#define Tempus_StepperStaggeredForwardSensitivity_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_StaggeredForwardSensitivityModelEvaluator.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "NOX_Thyra.H"

namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;

// StepperStaggeredForwardSensitivity definitions:
template<class Scalar>
StepperStaggeredForwardSensitivity<Scalar>::
StepperStaggeredForwardSensitivity(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Teuchos::ParameterList>& pList,
  const Teuchos::RCP<Teuchos::ParameterList>& sens_pList)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Set all the input parameters and call initialize
  this->setParams(pList, sens_pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Create forward sensitivity model evaluator wrapper
  Teuchos::RCP<Teuchos::ParameterList> spl = Teuchos::parameterList();
  *spl = *sensPL_;
  spl->remove("Reuse State Linear Solver");
  spl->remove("Force W Update");
  fsa_model_ =
    rcp(new StaggeredForwardSensitivityModelEvaluator<Scalar>(
          appModel, spl));

  // Create state and sensitivity steppers
  RCP<StepperFactory<Scalar> > sf =Teuchos::rcp(new StepperFactory<Scalar>());
  if (stateStepper_ == Teuchos::null)
    stateStepper_ = sf->createStepper(appModel, stepperPL_);
  else
    stateStepper_->setModel(appModel);
  if (sensitivityStepper_ == Teuchos::null)
    sensitivityStepper_ = sf->createStepper(fsa_model_, stepperPL_);
  else
    sensitivityStepper_->setModel(fsa_model_);
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
setSolver(std::string solverName)
{
  stateStepper_->setSolver(solverName);
  sensitivityStepper_->setSolver(solverName);
}

template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
StepperStaggeredForwardSensitivity<Scalar>::
getModel()
{
  return stateStepper_->getModel();
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  stateStepper_->setSolver(solverPL);
  sensitivityStepper_->setSolver(solverPL);
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  stateStepper_->setSolver(solver);
  sensitivityStepper_->setSolver(solver);
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
initialize()
{
  this->setSolver();
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::MultiVectorBase;
  using Thyra::assign;
  using Thyra::createMember;
  using Thyra::multiVectorProductVector;
  using Thyra::multiVectorProductVectorSpace;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;

  // Initialize state, sensitivity solution histories if necessary.
  // We only need to split the solution history into state and sensitivity
  // components for the first step, otherwise the state and sensitivity
  // histories are updated from the previous step.
  if (stateSolutionHistory_ == Teuchos::null) {
    RCP<Teuchos::ParameterList> shPL =
      solutionHistory->getNonconstParameterList();

    // Get product X, XDot, XDotDot
    RCP<SolutionState<Scalar> > state = solutionHistory->getCurrentState();
    RCP<DMVPV> X, XDot, XDotDot;
    X = rcp_dynamic_cast<DMVPV>(state->getX());
    XDot = rcp_dynamic_cast<DMVPV>(state->getXDot());
    if (state->getXDotDot() != Teuchos::null)
      XDotDot = rcp_dynamic_cast<DMVPV>(state->getXDotDot());

    // Pull out state components (has to be non-const because of SolutionState
    // constructor)
    RCP<VectorBase<Scalar> > x, xdot, xdotdot;
    x = X->getNonconstMultiVector()->col(0);
    xdot = XDot->getNonconstMultiVector()->col(0);
    if (XDotDot != Teuchos::null)
      xdotdot = XDotDot->getNonconstMultiVector()->col(0);

    // Create state solution history
    RCP<SolutionState<Scalar> > state_state =
      rcp(new SolutionState<Scalar>(state->getMetaData()->clone(),
                                    x, xdot, xdotdot,
                                    state->getStepperState()->clone()));
    stateSolutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));
    stateSolutionHistory_->addState(state_state);

    const int num_param = X->getMultiVector()->domain()->dim()-1;
    const Teuchos::Range1D rng(1,num_param);

    // Pull out sensitivity components
    RCP<MultiVectorBase<Scalar> > dxdp, dxdotdp, dxdotdotdp;
    dxdp = X->getNonconstMultiVector()->subView(rng);
    dxdotdp = XDot->getNonconstMultiVector()->subView(rng);
    if (XDotDot != Teuchos::null)
      dxdotdotdp = XDotDot->getNonconstMultiVector()->subView(rng);

    // Create sensitivity product vectors
    RCP<DMVPV> dxdp_vec, dxdotdp_vec, dxdotdotdp_vec;
    RCP<const DMVPVS> prod_space =
      multiVectorProductVectorSpace(X->getMultiVector()->range(), num_param);
    dxdp_vec = multiVectorProductVector(prod_space, dxdp);
    dxdotdp_vec = multiVectorProductVector(prod_space, dxdotdp);
    if (XDotDot != Teuchos::null)
      dxdotdotdp_vec = multiVectorProductVector(prod_space, dxdotdotdp);

    // Create sensitivity solution history
    RCP<SolutionState<Scalar> > sens_state =
      rcp(new SolutionState<Scalar>(state->getMetaData()->clone(),
                                    dxdp_vec, dxdotdp_vec, dxdotdotdp_vec,
                                    state->getStepperState()->clone()));
    sensSolutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));
    sensSolutionHistory_->addState(sens_state);
  }

  // Take step for state equations
  stateSolutionHistory_->initWorkingState();
  stateSolutionHistory_->getWorkingState()->getMetaData()->copy(
    solutionHistory->getWorkingState()->getMetaData());
  stateStepper_->takeStep(stateSolutionHistory_);

  // Get current state and set in sensitivity model evaluator
   RCP<SolutionState<Scalar> > state =
     stateSolutionHistory_->getWorkingState();
  RCP<const VectorBase<Scalar> > x = state->getX();
  RCP<const VectorBase<Scalar> > xdot = state->getXDot();
  RCP<const VectorBase<Scalar> > xdotdot = state->getXDotDot();
  fsa_model_->setModelX(x);
  fsa_model_->setModelXDot(xdot);
  fsa_model_->setModelXDotDot(xdotdot);

  // Reuse state solver if requested
  if (reuse_solver_ && stateStepper_->getSolver() != Teuchos::null) {
    RCP<Thyra::NOXNonlinearSolver> nox_solver =
      rcp_dynamic_cast<Thyra::NOXNonlinearSolver>(stateStepper_->getSolver());
    fsa_model_->setLO(nox_solver->get_nonconst_W_op(force_W_update_));
    fsa_model_->setPO(nox_solver->get_nonconst_prec_op());
  }

  // Take step in senstivity equations
  sensSolutionHistory_->initWorkingState();
  sensSolutionHistory_->getWorkingState()->getMetaData()->copy(
    solutionHistory->getWorkingState()->getMetaData());
  sensitivityStepper_->takeStep(sensSolutionHistory_);

  // Get current sensitivities
  RCP<const MultiVectorBase<Scalar> > dxdp, dxdotdp, dxdotdotdp;
  RCP<SolutionState<Scalar> > sens_state =
    sensSolutionHistory_->getWorkingState();
  dxdp = rcp_dynamic_cast<const DMVPV>(
    sens_state->getX())->getMultiVector();
  dxdotdp = rcp_dynamic_cast<const DMVPV>(
    sens_state->getXDot())->getMultiVector();
  if (sens_state->getXDotDot() != Teuchos::null)
    dxdotdotdp = rcp_dynamic_cast<const DMVPV>(
      sens_state->getXDotDot())->getMultiVector();

  const int num_param = dxdp->domain()->dim();
  const Teuchos::Range1D rng(1,num_param);

  // Form product solutions
  RCP<DMVPV> X, XDot, XDotDot;
  RCP<const Thyra::VectorSpaceBase<Scalar> > prod_space =
    multiVectorProductVectorSpace(dxdp->range(), num_param+1);
  X = rcp_dynamic_cast<DMVPV>(createMember(prod_space));
  assign(X->getNonconstMultiVector()->col(0).ptr(), *x);
  assign(X->getNonconstMultiVector()->subView(rng).ptr(), *dxdp);
  XDot = rcp_dynamic_cast<DMVPV>(createMember(prod_space));
  assign(XDot->getNonconstMultiVector()->col(0).ptr(), *xdot);
  assign(XDot->getNonconstMultiVector()->subView(rng).ptr(), *dxdotdp);
  if (xdotdot != Teuchos::null) {
    XDotDot = rcp_dynamic_cast<DMVPV>(createMember(prod_space));
    assign(XDotDot->getNonconstMultiVector()->col(0).ptr(), *xdot);
    assign(XDotDot->getNonconstMultiVector()->subView(rng).ptr(), *dxdotdotdp);
  }

  // Add step to solution history
  RCP<SolutionState<Scalar> > prod_state = solutionHistory->getWorkingState();
  assign(prod_state->getX().ptr(), *X);
  assign(prod_state->getXDot().ptr(), *XDot);
  if (XDotDot != Teuchos::null)
    assign(prod_state->getXDotDot().ptr(), *XDotDot);
  prod_state->setOrder(std::min(state->getOrder(), sens_state->getOrder()));

  // Determine whether step passed or failed
  bool passed = true;
  if (state->getStepperStatus() == Status::FAILED ||
      state->getSolutionStatus() == Status::FAILED ||
      sens_state->getStepperStatus() == Status::FAILED ||
      sens_state->getSolutionStatus() == Status::FAILED)
    passed = false;

  if (passed) {
    prod_state->getStepperState()->stepperStatus_ = Status::PASSED;
    stateSolutionHistory_->promoteWorkingState();
    sensSolutionHistory_->promoteWorkingState();
  }
  else
     prod_state->getStepperState()->stepperStatus_ = Status::FAILED;
}


template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperStaggeredForwardSensitivity<Scalar>::
getDefaultStepperState()
{
  // ETP:  Note, maybe this should be a special stepper state that combines
  // states from both state and sensitivity steppers?
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperStaggeredForwardSensitivity<Scalar>::
description() const
{
  std::string name = "StepperStaggeredForwardSensitivity";
  return(name);
}


template<class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl;
  stateStepper_->describe(out, verbLevel);
  out << std::endl;
  sensitivityStepper_->describe(out, verbLevel);
}


template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  if (pList == Teuchos::null)
    stepperPL_ = this->getDefaultParameters();
  else
    stepperPL_ = pList;
  // Can not validate because of optional Parameters (e.g., Solver Name).
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperStaggeredForwardSensitivity<Scalar>::
getValidParameters() const
{
  return stateStepper_->getValidParameters();
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperStaggeredForwardSensitivity<Scalar>::
getDefaultParameters() const
{
  return stateStepper_->getDefaultParameters();
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperStaggeredForwardSensitivity<Scalar>::
getNonconstParameterList()
{
  return stepperPL_;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperStaggeredForwardSensitivity<Scalar>::
unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return temp_plist;
}


template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::
setParams(
  Teuchos::RCP<Teuchos::ParameterList> const& pList,
  Teuchos::RCP<Teuchos::ParameterList> const& spList)
{
  if (pList == Teuchos::null)
    stepperPL_ = this->getDefaultParameters();
  else
    stepperPL_ = pList;

  if (spList == Teuchos::null)
    sensPL_ = Teuchos::parameterList();
  else
    sensPL_ = spList;

  reuse_solver_ = sensPL_->get("Reuse State Linear Solver", false);
  force_W_update_ = sensPL_->get("Force W Update", true);

  // Can not validate because of optional Parameters (e.g., Solver Name).
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
}


template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StepperStaggeredForwardSensitivity<Scalar>::
get_x_space() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;

  RCP<const Thyra::VectorSpaceBase<Scalar> > x_space =
    stateStepper_->getModel()->get_x_space();
  RCP<const DMVPVS> dxdp_space =
    rcp_dynamic_cast<const DMVPVS>(sensitivityStepper_->getModel()->get_x_space());
  const int num_param = dxdp_space->numBlocks();
  RCP<const Thyra::VectorSpaceBase<Scalar> > prod_space =
    multiVectorProductVectorSpace(x_space, num_param+1);
  return prod_space;
}

} // namespace Tempus

#endif // Tempus_StepperStaggeredForwardSensitivity_impl_hpp
