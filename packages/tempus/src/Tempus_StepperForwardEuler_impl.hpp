// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperForwardEuler_impl_hpp
#define Tempus_StepperForwardEuler_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Tempus_StepperForwardEulerModifierDefault.hpp"

namespace Tempus {

template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler()
{
  this->setStepperType(        "Forward Euler");
  this->setUseFSAL(            true);
  this->setICConsistency(      "Consistent");
  this->setICConsistencyCheck( false);
  this->setAppAction(Teuchos::null);
}

template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  const Teuchos::RCP<StepperForwardEulerAppAction<Scalar> >& stepperFEAppAction)
{
  this->setStepperType(        "Forward Euler");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);

  this->setAppAction(stepperFEAppAction);
  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setAppAction(
  Teuchos::RCP<StepperForwardEulerAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperFEAppAction_ =
      Teuchos::rcp(new StepperForwardEulerModifierDefault<Scalar>());
  }
  else {
    stepperFEAppAction_ = appAction;
  }
}


template<class Scalar>
void StepperForwardEuler<Scalar>::setInitialConditions(
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

template<class Scalar>
void StepperForwardEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperForwardEuler::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperForwardEuler<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for Forward Euler.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<StepperForwardEuler<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperFEAppAction_->execute(solutionHistory, thisStepper,
      StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    if (currentState->getXDot() != Teuchos::null)
      this->setStepperXDot(currentState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();
    const Scalar dt = workingState->getTimeStep();

    if (!(this->getUseFSAL()) || workingState->getNConsecutiveFailures() != 0) {
      // Need to compute XDotOld.
      stepperFEAppAction_->execute(solutionHistory, thisStepper,
        StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);

      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, currentState->getX(),
                                currentState->getTime(), p);

      // For UseFSAL=false, x and xDot are now sync'ed or consistent
      // at the same time level for the currentState.
      currentState->setIsSynced(true);
    }


    // Forward Euler update, x^n = x^{n-1} + dt^n * xDot^{n-1}
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
      *(currentState->getX()),dt,*(xDot));


    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    xDot = this->getStepperXDot();

    if (this->getUseFSAL()) {
      // Get consistent xDot^n.
      stepperFEAppAction_->execute(solutionHistory, thisStepper,
        StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);

      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, workingState->getX(),
                                workingState->getTime(), p);

      // For UseFSAL=true, x and xDot are now sync'ed or consistent
      // for the workingState.
      workingState->setIsSynced(true);
    } else {
      assign(xDot.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      workingState->setIsSynced(false);
    }

    workingState->setSolutionStatus(Status::PASSED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperFEAppAction_->execute(solutionHistory, thisStepper,
      StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}


/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> > StepperForwardEuler<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperForwardEuler<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperExplicit<Scalar>::describe(out, verbLevel);
  out << "  stepperFEAppAction_                = "
      << stepperFEAppAction_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperForwardEuler<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperExplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if (stepperFEAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Forward Euler AppAction is not set!\n";
  }
  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<bool>("Use FSAL", true);
  pl->set<std::string>("Initial Condition Consistency", "Consistent");
  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperForwardEuler_impl_hpp
