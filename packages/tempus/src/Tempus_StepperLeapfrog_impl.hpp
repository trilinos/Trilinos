//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrog_impl_hpp
#define Tempus_StepperLeapfrog_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_StepperLeapfrogModifierDefault.hpp"

namespace Tempus {

template <class Scalar>
StepperLeapfrog<Scalar>::StepperLeapfrog()
{
  this->setStepperName("Leapfrog");
  this->setStepperType("Leapfrog");
  this->setUseFSAL(false);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);

  this->setAppAction(Teuchos::null);
}

template <class Scalar>
StepperLeapfrog<Scalar>::StepperLeapfrog(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    const Teuchos::RCP<StepperLeapfrogAppAction<Scalar> >& stepperLFAppAction)
{
  this->setStepperName("Leapfrog");
  this->setStepperType("Leapfrog");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setAppAction(stepperLFAppAction);
  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperLeapfrog<Scalar>::setAppAction(
    Teuchos::RCP<StepperLeapfrogAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperLFAppAction_ =
        Teuchos::rcp(new StepperLeapfrogModifierDefault<Scalar>());
  }
  else {
    stepperLFAppAction_ = appAction;
  }
}

template <class Scalar>
void StepperLeapfrog<Scalar>::setInitialConditions(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDotDot
  if (initialState->getXDotDot() == Teuchos::null)
    this->setStepperXDotDot(initialState->getX()->clone_v());
  else
    this->setStepperXDotDot(initialState->getXDotDot());

  StepperExplicit<Scalar>::setInitialConditions(solutionHistory);

  if (this->getUseFSAL()) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out, 1, "StepperLeapfrog::setInitialConditions()");
    *out << "Warning -- The First-Same-As-Last (FSAL) principle is not "
         << "used with Leapfrog because of the algorithm's prescribed "
         << "order of solution update. The default is to set useFSAL=false, "
         << "however useFSAL=true will also work but have no affect "
         << "(i.e., no-op).\n"
         << std::endl;
  }
}

template <class Scalar>
void StepperLeapfrog<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperLeapfrog::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperLeapfrog<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for Leapfrog.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\nTry setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n"
            << "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");

    RCP<SolutionState<Scalar> > currentState =
        solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    const Scalar time = currentState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    RCP<StepperLeapfrog<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);

    stepperLFAppAction_->execute(
        solutionHistory, thisStepper,
        StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    // Perform half-step startup if working state is synced
    // (i.e., xDot and x are at the same time level).
    if (workingState->getIsSynced() == true) {
      // Half-step startup: xDot_{n+1/2} = xDot_n + 0.5*dt*xDotDot_n
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
                     *(currentState->getXDot()), 0.5 * dt,
                     *(currentState->getXDotDot()));
    }
    stepperLFAppAction_->execute(
        solutionHistory, thisStepper,
        StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::BEFORE_X_UPDATE);
    // x_{n+1} = x_n + dt*xDot_{n+1/2}
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
                   *(currentState->getX()), dt, *(workingState->getXDot()));

    stepperLFAppAction_->execute(
        solutionHistory, thisStepper,
        StepperLeapfrogAppAction<
            Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);
    auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

    // Evaluate xDotDot = f(x,t).
    this->evaluateExplicitODE(workingState->getXDotDot(), workingState->getX(),
                              Teuchos::null, time + dt, p);
    stepperLFAppAction_->execute(
        solutionHistory, thisStepper,
        StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::BEFORE_XDOT_UPDATE);
    if (workingState->getOutput() == true) {
      // Half-step sync: xDot_{n+1} = xDot_{n+1/2} + 0.5*dt*xDotDot_{n+1}
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
                     *(workingState->getXDot()), 0.5 * dt,
                     *(workingState->getXDotDot()));
      workingState->setIsSynced(true);
    }
    else {
      // Full leapfrog step: xDot_{n+3/2} = xDot_{n+1/2} + dt*xDotDot_{n+1}
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
                     *(workingState->getXDot()), dt,
                     *(workingState->getXDotDot()));
      workingState->setIsSynced(false);
    }

    workingState->setSolutionStatus(Status::PASSED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);

    stepperLFAppAction_->execute(
        solutionHistory, thisStepper,
        StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
StepperLeapfrog<Scalar>::getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperLeapfrog<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperExplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperLeapfrog ---\n";
  out << "  stepperLFAppAction_                = " << stepperLFAppAction_
      << std::endl;
  out << "-----------------------" << std::endl;
}

template <class Scalar>
bool StepperLeapfrog<Scalar>::isValidSetup(Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;
  if (!StepperExplicit<Scalar>::isValidSetup(out)) isValidSetup = false;
  if (stepperLFAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Leapfrog AppAction is not set!\n";
  }

  return isValidSetup;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperLeapfrog<Scalar> > createStepperLeapfrog(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperLeapfrog<Scalar>());
  stepper->setStepperExplicitValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperLeapfrog_impl_hpp
