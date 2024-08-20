//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperTrapezoidal_impl_hpp
#define Tempus_StepperTrapezoidal_impl_hpp

#include "Tempus_StepperTrapezoidalModifierDefault.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"

namespace Tempus {

template <class Scalar>
StepperTrapezoidal<Scalar>::StepperTrapezoidal()
{
  this->setStepperName("Trapezoidal Method");
  this->setStepperType("Trapezoidal Method");
  this->setUseFSAL(true);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}

template <class Scalar>
StepperTrapezoidal<Scalar>::StepperTrapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperTrapezoidalAppAction<Scalar> >&
        stepperTrapAppAction)
{
  this->setStepperName("Trapezoidal Method");
  this->setStepperType("Trapezoidal Method");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);

  this->setAppAction(stepperTrapAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperTrapezoidal<Scalar>::setAppAction(
    Teuchos::RCP<StepperTrapezoidalAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperTrapAppAction_ =
        Teuchos::rcp(new StepperTrapezoidalModifierDefault<Scalar>());
  }
  else {
    stepperTrapAppAction_ = appAction;
  }
  this->isInitialized_ = false;
}

template <class Scalar>
void StepperTrapezoidal<Scalar>::setInitialConditions(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDot
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(initialState->getX()->clone_v());
  else
    this->setStepperXDot(initialState->getXDot());

  StepperImplicit<Scalar>::setInitialConditions(solutionHistory);

  TEUCHOS_TEST_FOR_EXCEPTION(
      !(this->getUseFSAL()), std::logic_error,
      "Error - The First-Same-As-Last (FSAL) principle is required\n"
      "        for the Trapezoidal Stepper (i.e., useFSAL=true)!\n");
  //   There are at least two ways around this, but are not implemented.
  //    - Do a solve for xDotOld, xDot_{n-1}, at each time step as for the
  //      initial conditions.  This is expensive since you would be doing
  //      two solves every time step.
  //    - Use evaluateExplicitODE to get xDot_{n-1} if the application
  //      provides it.  Explicit evaluations are cheaper but requires the
  //      application to implement xDot = f(x,t).
}

template <class Scalar>
void StepperTrapezoidal<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperTrapezoidal::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperTrapezoidal<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for Trapezoidal.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\nTry setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n"
            << "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");
    RCP<StepperTrapezoidal<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperTrapAppAction_->execute(
        solutionHistory, thisStepper,
        StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState =
        solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar> > xOld    = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > xDotOld = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > x             = workingState->getX();
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();

    const Scalar time  = workingState->getTime();
    const Scalar dt    = workingState->getTimeStep();
    const Scalar alpha = getAlpha(dt);
    const Scalar beta  = getBeta(dt);

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer = Teuchos::rcp(
        new StepperTrapezoidalTimeDerivative<Scalar>(alpha, xOld, xDotOld));

    auto p = Teuchos::rcp(
        new ImplicitODEParameters<Scalar>(timeDer, dt, alpha, beta));
    stepperTrapAppAction_->execute(
        solutionHistory, thisStepper,
        StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

    const Thyra::SolveStatus<Scalar> sStatus =
        this->solveImplicitODE(x, xDot, time, p);
    stepperTrapAppAction_->execute(
        solutionHistory, thisStepper,
        StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);

    if (workingState->getXDot() != Teuchos::null) timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperTrapAppAction_->execute(
        solutionHistory, thisStepper,
        StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
StepperTrapezoidal<Scalar>::getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperTrapezoidal<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperTrapezoidal ---\n";
  out << "  stepperTrapAppAction_ = " << stepperTrapAppAction_ << std::endl;
  out << "--------------------------" << std::endl;
}

template <class Scalar>
bool StepperTrapezoidal<Scalar>::isValidSetup(Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;
  if (!StepperImplicit<Scalar>::isValidSetup(out)) isValidSetup = false;

  if (stepperTrapAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Trapezoidal AppAction is not set!\n";
  }

  return isValidSetup;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperTrapezoidal<Scalar> > createStepperTrapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperTrapezoidal<Scalar>());
  stepper->setStepperImplicitValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperTrapezoidal_impl_hpp
