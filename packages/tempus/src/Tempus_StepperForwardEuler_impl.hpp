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


namespace Tempus {


template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler()
{
  this->setStepperType(        "Forward Euler");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());

  this->setObserver();
}


template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<StepperObserver<Scalar> >& obs,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck)
{
  this->setStepperType(        "Forward Euler");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);

  this->setObserver(obs);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
void StepperForwardEuler<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (this->stepperObserver_ == Teuchos::null) {
      stepperFEObserver_ =
        Teuchos::rcp(new StepperForwardEulerObserver<Scalar>());
      this->stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >(stepperFEObserver_,true);
    }
  } else {
    this->stepperObserver_ = obs;
    stepperFEObserver_ =
      Teuchos::rcp_dynamic_cast<StepperForwardEulerObserver<Scalar> >
        (this->stepperObserver_,true);
  }

  this->isInitialized_ = false;
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

    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();

    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot(currentState);
    const Scalar dt = workingState->getTimeStep();

    if ( !(this->getUseFSAL()) ) {
      // Need to compute XDotOld.
      if (!Teuchos::is_null(stepperFEObserver_))
        stepperFEObserver_->observeBeforeExplicit(solutionHistory, *this);

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


    xDot = this->getStepperXDot(workingState);

    if (this->getUseFSAL()) {
      // Get consistent xDot^n.
      if (!Teuchos::is_null(stepperFEObserver_))
        stepperFEObserver_->observeBeforeExplicit(solutionHistory, *this);

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
    this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
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

  out << "--- StepperForwardEuler ---\n";
  out << stepperFEObserver_ << std::endl;
  out << "---------------------------" << std::endl;
}


template<class Scalar>
bool StepperForwardEuler<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperExplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (stepperFEObserver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Forward Euler observer is not set!\n";
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
