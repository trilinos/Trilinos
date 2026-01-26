//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI3_impl_hpp
#define Tempus_StepperEPI3_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_StepperEPI3ModifierDefault.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"

namespace Tempus {

template <class Scalar>
StepperEPI3<Scalar>::StepperEPI3()
{
  this->setStepperName("EPI3");
  this->setStepperType("EPI3");
  this->setUseFSAL(true);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setAppAction(Teuchos::null);
}

template <class Scalar>
StepperEPI3<Scalar>::StepperEPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPI3AppAction<Scalar> >&
        stepperEPI3AppAction)
{
  this->setStepperName("EPI3");
  this->setStepperType("EPI3");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);

  this->setAppAction(stepperEPI3AppAction);
  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperEPI3<Scalar>::setAppAction(
    Teuchos::RCP<StepperEPI3AppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperEPI3AppAction_ =
        Teuchos::rcp(new StepperEPI3ModifierDefault<Scalar>());
  }
  else {
    stepperEPI3AppAction_ = appAction;
  }
}

template<class Scalar>
void StepperEPI3<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperExplicit<Scalar>::setModel(appModel);

  auto phif = Teuchos::rcp(new PhiEvaluatorFactory<Scalar>());
  phiEvaluator_ = phif->createPhiEvaluator("PFD", appModel);

  TEUCHOS_TEST_FOR_EXCEPTION(
  phiEvaluator_.is_null(), std::logic_error,
  "phiEvaluator_ is null in StepperEPI3::setModel");

  phiEvaluator_->setModel(appModel);
  phiEvaluator_->initialize();

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperEPI3<Scalar>::setInitialConditions(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDot
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(initialState->getX()->clone_v());
  else
    this->setStepperXDot(initialState->getXDot());

  StepperExplicit<Scalar>::setInitialConditions(solutionHistory);
}

template <class Scalar>
void StepperEPI3<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();
  phiEvaluator_->checkInitialized();
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperEPI3::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperEPI3<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for EPI3.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\n Try setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");

    RCP<StepperEPI3<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEPI3AppAction_->execute(
        solutionHistory, thisStepper,
        StepperEPI3AppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > currentState =
        solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    if (currentState->getXDot() != Teuchos::null)
      this->setStepperXDot(currentState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();
    const Scalar dt                      = workingState->getTimeStep();

    // if (!(this->getUseFSAL()) || workingState->getNConsecutiveFailures() != 0) {
      // Need to compute XDotOld.
      stepperEPI3AppAction_->execute(
          solutionHistory, thisStepper,
          StepperEPI3AppAction<
              Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);

      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, currentState->getX(),
                                currentState->getTime(), p);

      RCP<const Thyra::VectorBase<Scalar> > xDot_Exponential = phiEvaluator_->buildATildeMatrix(2, dt, xDot);

      // For UseFSAL=false, x and xDot are now sync'ed or consistent
      // at the same time level for the currentState.
      currentState->setIsSynced(true);
    // }

    // EPI3 update, x^n = x^{n-1} + dt^n * xDot^{n-1}
    // Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
    //                *(currentState->getX()), dt, *(xDot));
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
                   *(currentState->getX()), dt, *(xDot_Exponential));

    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    xDot = this->getStepperXDot();

    if (this->getUseFSAL()) {
      // Get consistent xDot^n.
      stepperEPI3AppAction_->execute(
          solutionHistory, thisStepper,
          StepperEPI3AppAction<
              Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);

      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, workingState->getX(),
                                workingState->getTime(), p);

      // For UseFSAL=true, x and xDot are now sync'ed or consistent
      // for the workingState.
      workingState->setIsSynced(true);
    }
    else {
      assign(xDot.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      workingState->setIsSynced(false);
    }

    workingState->setSolutionStatus(Status::PASSED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperEPI3AppAction_->execute(
        solutionHistory, thisStepper,
        StepperEPI3AppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}

/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template <class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperEPI3<Scalar>::getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperEPI3<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << std::endl;
  Stepper<Scalar>::describe(*l_out, verbLevel);
  StepperExplicit<Scalar>::describe(*l_out, verbLevel);
  *l_out << "  stepperEPI3AppAction_ = " << stepperEPI3AppAction_ << std::endl
         << "----------------------------" << std::endl;
}

template <class Scalar>
bool StepperEPI3<Scalar>::isValidSetup(Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);

  bool isValidSetup = true;

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;
  if (!StepperExplicit<Scalar>::isValidSetup(out)) isValidSetup = false;
  if (stepperEPI3AppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The EPI3 AppAction is not set!\n";
  }
  return isValidSetup;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperEPI3<Scalar> > createStepperEPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperEPI3<Scalar>());
  stepper->setStepperExplicitValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperEPI3_impl_hpp
