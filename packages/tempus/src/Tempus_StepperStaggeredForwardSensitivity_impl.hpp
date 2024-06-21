//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperStaggeredForwardSensitivity_impl_hpp
#define Tempus_StepperStaggeredForwardSensitivity_impl_hpp

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapStaggeredFSAModelEvaluator.hpp"
#include "Tempus_WrapCombinedFSAModelEvaluator.hpp"

namespace Tempus {

template <class Scalar>
StepperStaggeredForwardSensitivity<Scalar>::StepperStaggeredForwardSensitivity()
  : stepMode_(SensitivityStepMode::Forward)
{
  this->setStepperName("StaggeredForwardSensitivity");
  this->setStepperType("StaggeredForwardSensitivity");
  this->setParams(Teuchos::null, Teuchos::null);
}

template <class Scalar>
StepperStaggeredForwardSensitivity<Scalar>::StepperStaggeredForwardSensitivity(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
        sens_residual_model,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& sens_solve_model,
    const Teuchos::RCP<Teuchos::ParameterList>& pList,
    const Teuchos::RCP<Teuchos::ParameterList>& sens_pList)
  : stepMode_(SensitivityStepMode::Forward)
{
  // Set all the input parameters and call initialize
  this->setStepperName("StaggeredForwardSensitivity");
  this->setStepperType("StaggeredForwardSensitivity");
  this->setParams(pList, sens_pList);
  this->setModel(appModel, sens_residual_model, sens_solve_model);
  this->initialize();
}

template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
        sens_residual_model,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& sens_solve_model)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Create forward sensitivity model evaluator wrapper
  Teuchos::RCP<Teuchos::ParameterList> spl = Teuchos::parameterList();
  *spl                                     = *sensPL_;
  spl->remove("Reuse State Linear Solver");
  spl->remove("Force W Update");
  fsa_model_ = wrapStaggeredFSAModelEvaluator(appModel, sens_residual_model,
                                              sens_solve_model, false, spl);

  // Create combined FSA ME which serves as "the" ME for this stepper,
  // so that getModel() has a ME consistent the FSA problem (including both
  // state and sensitivity components), e.g., the integrator may call
  // getModel()->getNominalValues(), which needs to be consistent.
  combined_fsa_model_ = wrapCombinedFSAModelEvaluator(
      appModel, sens_residual_model, sens_solve_model, spl);

  // Create state and sensitivity steppers
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
  if (stateStepper_ == Teuchos::null)
    stateStepper_ = sf->createStepper(stepperPL_, appModel);
  else
    stateStepper_->setModel(appModel);

  if (sensitivityStepper_ == Teuchos::null)
    sensitivityStepper_ = sf->createStepper(stepperPL_, fsa_model_);
  else
    sensitivityStepper_->setModel(fsa_model_);

  this->isInitialized_ = false;
}

template <class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
StepperStaggeredForwardSensitivity<Scalar>::getModel() const
{
  return combined_fsa_model_;
}

template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::setSolver(
    Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  stateStepper_->setSolver(solver);
  sensitivityStepper_->setSolver(solver);

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::createMember;
  using Thyra::MultiVectorBase;
  using Thyra::multiVectorProductVector;
  using Thyra::multiVectorProductVectorSpace;
  using Thyra::VectorBase;
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
    X = rcp_dynamic_cast<DMVPV>(state->getX(), true);

    XDot = rcp_dynamic_cast<DMVPV>(state->getXDot(), true);
    if (state->getXDotDot() != Teuchos::null)
      XDotDot = rcp_dynamic_cast<DMVPV>(state->getXDotDot(), true);

    // Pull out state components (has to be non-const because of SolutionState
    // constructor)
    RCP<VectorBase<Scalar> > x, xdot, xdotdot;
    x    = X->getNonconstMultiVector()->col(0);
    xdot = XDot->getNonconstMultiVector()->col(0);
    if (XDotDot != Teuchos::null)
      xdotdot = XDotDot->getNonconstMultiVector()->col(0);

    // Create state solution history
    RCP<SolutionState<Scalar> > state_state = state->clone();
    state_state->setX(x);
    state_state->setXDot(xdot);
    state_state->setXDotDot(xdotdot);
    stateSolutionHistory_ = createSolutionHistoryPL<Scalar>(shPL);
    stateSolutionHistory_->addState(state_state);

    const int num_param = X->getMultiVector()->domain()->dim() - 1;
    TEUCHOS_ASSERT(num_param > 0);
    const Teuchos::Range1D rng(1, num_param);

    // Pull out sensitivity components
    RCP<MultiVectorBase<Scalar> > dxdp, dxdotdp, dxdotdotdp;
    dxdp    = X->getNonconstMultiVector()->subView(rng);
    dxdotdp = XDot->getNonconstMultiVector()->subView(rng);
    if (XDotDot != Teuchos::null)
      dxdotdotdp = XDotDot->getNonconstMultiVector()->subView(rng);

    // Create sensitivity product vectors
    RCP<VectorBase<Scalar> > dxdp_vec, dxdotdp_vec, dxdotdotdp_vec;
    RCP<const DMVPVS> prod_space =
        multiVectorProductVectorSpace(X->getMultiVector()->range(), num_param);
    dxdp_vec    = multiVectorProductVector(prod_space, dxdp);
    dxdotdp_vec = multiVectorProductVector(prod_space, dxdotdp);
    if (XDotDot != Teuchos::null)
      dxdotdotdp_vec = multiVectorProductVector(prod_space, dxdotdotdp);

    // Create sensitivity solution history
    RCP<SolutionState<Scalar> > sens_state = state->clone();
    sens_state->setX(dxdp_vec);
    sens_state->setXDot(dxdotdp_vec);
    sens_state->setXDotDot(dxdotdotdp_vec);
    sensSolutionHistory_ = createSolutionHistoryPL<Scalar>(shPL);
    sensSolutionHistory_->addState(sens_state);
  }

  // Get our working state
  RCP<SolutionState<Scalar> > prod_state = solutionHistory->getWorkingState();
  RCP<DMVPV> X, XDot, XDotDot;
  X    = rcp_dynamic_cast<DMVPV>(prod_state->getX(), true);
  XDot = rcp_dynamic_cast<DMVPV>(prod_state->getXDot(), true);
  if (prod_state->getXDotDot() != Teuchos::null)
    XDotDot = rcp_dynamic_cast<DMVPV>(prod_state->getXDotDot(), true);

  // Take step for state equations
  stepMode_ = SensitivityStepMode::Forward;
  stateSolutionHistory_->initWorkingState();
  RCP<SolutionState<Scalar> > state = stateSolutionHistory_->getWorkingState();
  state->getMetaData()->copy(prod_state->getMetaData());
  stateStepper_->takeStep(stateSolutionHistory_);

  // Set state components of product state
  assign(X->getNonconstMultiVector()->col(0).ptr(), *(state->getX()));
  assign(XDot->getNonconstMultiVector()->col(0).ptr(), *(state->getXDot()));
  if (XDotDot != Teuchos::null)
    assign(XDotDot->getNonconstMultiVector()->col(0).ptr(),
           *(state->getXDotDot()));
  prod_state->setOrder(state->getOrder());

  // If step passed promote the state, otherwise fail and stop
  if (state->getSolutionStatus() == Status::FAILED) {
    prod_state->setSolutionStatus(Status::FAILED);
    return;
  }
  stateSolutionHistory_->promoteWorkingState();

  // Get forward state in sensitivity model evaluator
  fsa_model_->setForwardSolutionHistory(stateSolutionHistory_);
  if (reuse_solver_ && stateStepper_->getSolver() != Teuchos::null)
    fsa_model_->setSolver(stateStepper_->getSolver(), force_W_update_);

  // Take step in sensitivity equations
  stepMode_ = SensitivityStepMode::Sensitivity;
  sensSolutionHistory_->initWorkingState();
  RCP<SolutionState<Scalar> > sens_state =
      sensSolutionHistory_->getWorkingState();
  sens_state->getMetaData()->copy(prod_state->getMetaData());
  sensitivityStepper_->takeStep(sensSolutionHistory_);

  // Set sensitivity components of product state
  RCP<const MultiVectorBase<Scalar> > dxdp =
      rcp_dynamic_cast<const DMVPV>(sens_state->getX(), true)->getMultiVector();
  const int num_param = dxdp->domain()->dim();
  const Teuchos::Range1D rng(1, num_param);
  assign(X->getNonconstMultiVector()->subView(rng).ptr(), *dxdp);
  RCP<const MultiVectorBase<Scalar> > dxdotdp =
      rcp_dynamic_cast<const DMVPV>(sens_state->getXDot(), true)
          ->getMultiVector();
  assign(XDot->getNonconstMultiVector()->subView(rng).ptr(), *dxdotdp);
  if (sens_state->getXDotDot() != Teuchos::null) {
    RCP<const MultiVectorBase<Scalar> > dxdotdotdp =
        rcp_dynamic_cast<const DMVPV>(sens_state->getXDotDot(), true)
            ->getMultiVector();
    assign(XDotDot->getNonconstMultiVector()->subView(rng).ptr(), *dxdotdotdp);
  }
  prod_state->setOrder(std::min(state->getOrder(), sens_state->getOrder()));

  // If step passed promote the state, otherwise fail and stop
  if (sens_state->getSolutionStatus() == Status::FAILED) {
    prod_state->setSolutionStatus(Status::FAILED);
  }
  else {
    sensSolutionHistory_->promoteWorkingState();
    prod_state->setSolutionStatus(Status::PASSED);
  }
}

template <class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperStaggeredForwardSensitivity<Scalar>::getDefaultStepperState()
{
  // ETP:  Note, maybe this should be a special stepper state that combines
  // states from both state and sensitivity steppers?
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(description()));
  return stepperState;
}

template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);

  out << "--- StepperStaggeredForwardSensitivity ---\n";
  out << "  stateStepper_         = " << stateStepper_ << std::endl;
  out << "  sensitivityStepper_   = " << sensitivityStepper_ << std::endl;
  out << "  combined_fsa_model_   = " << combined_fsa_model_ << std::endl;
  out << "  fsa_model_            = " << fsa_model_ << std::endl;
  out << "  stateSolutionHistory_ = " << stateSolutionHistory_ << std::endl;
  out << "  sensSolutionHistory_  = " << sensSolutionHistory_ << std::endl;
  out << "------------------------------------------" << std::endl;

  out << description() << "::describe:" << std::endl;
  stateStepper_->describe(out, verbLevel);
  out << std::endl;
  sensitivityStepper_->describe(out, verbLevel);
}

template <class Scalar>
bool StepperStaggeredForwardSensitivity<Scalar>::isValidSetup(
    Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;

  if (stateStepper_ == Teuchos::null) {
    isValidSetup = false;
    out << "The state stepper is not set!\n";
  }

  if (sensitivityStepper_ == Teuchos::null) {
    isValidSetup = false;
    out << "The sensitivity stepper is not set!\n";
  }

  if (combined_fsa_model_ == Teuchos::null) {
    isValidSetup = false;
    out << "The combined FSA model is not set!\n";
  }

  if (fsa_model_ == Teuchos::null) {
    isValidSetup = false;
    out << "The FSA model is not set!\n";
  }

  return isValidSetup;
}

template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (this->stepperPL_ == Teuchos::null)
      this->stepperPL_ = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
          this->getValidParameters());
  }
  else {
    this->stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters (e.g., Solver Name).
  // stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperStaggeredForwardSensitivity<Scalar>::getValidParameters() const
{
  return stateStepper_->getValidParameters();
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperStaggeredForwardSensitivity<Scalar>::getNonconstParameterList()
{
  return stepperPL_;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperStaggeredForwardSensitivity<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_                                      = Teuchos::null;
  return temp_plist;
}

template <class Scalar>
void StepperStaggeredForwardSensitivity<Scalar>::setParams(
    Teuchos::RCP<Teuchos::ParameterList> const& pList,
    Teuchos::RCP<Teuchos::ParameterList> const& spList)
{
  if (pList == Teuchos::null)
    stepperPL_ = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
        this->getValidParameters());
  else
    stepperPL_ = pList;

  if (spList == Teuchos::null)
    sensPL_ = Teuchos::parameterList();
  else
    sensPL_ = spList;

  reuse_solver_   = sensPL_->get("Reuse State Linear Solver", false);
  force_W_update_ = sensPL_->get("Force W Update", true);

  // Can not validate because of optional Parameters (e.g., Solver Name).
  // stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StepperStaggeredForwardSensitivity<Scalar>::get_x_space() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;

  RCP<const Thyra::VectorSpaceBase<Scalar> > x_space =
      stateStepper_->getModel()->get_x_space();
  RCP<const DMVPVS> dxdp_space = rcp_dynamic_cast<const DMVPVS>(
      sensitivityStepper_->getModel()->get_x_space(), true);
  const int num_param = dxdp_space->numBlocks();
  RCP<const Thyra::VectorSpaceBase<Scalar> > prod_space =
      multiVectorProductVectorSpace(x_space, num_param + 1);
  return prod_space;
}

}  // namespace Tempus

#endif  // Tempus_StepperStaggeredForwardSensitivity_impl_hpp
