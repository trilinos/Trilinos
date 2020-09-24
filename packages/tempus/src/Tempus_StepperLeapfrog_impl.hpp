// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperLeapfrog_impl_hpp
#define Tempus_StepperLeapfrog_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Tempus_StepperLeapfrogModifierDefault.hpp"

namespace Tempus {


template<class Scalar>
StepperLeapfrog<Scalar>::StepperLeapfrog()
{
  this->setStepperType(        "Leapfrog");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  this->setObserver();
#endif
  this->setAppAction(Teuchos::null);
}

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
template<class Scalar>
StepperLeapfrog<Scalar>::StepperLeapfrog(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<StepperObserver<Scalar> >& obs,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck)
{
  this->setStepperType(        "Leapfrog");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);

  this->setObserver(obs);
  this->setAppAction(Teuchos::null);

  if (appModel != Teuchos::null) {

    this->setModel(appModel);
    this->initialize();
  }
}
#endif

template<class Scalar>
StepperLeapfrog<Scalar>::StepperLeapfrog(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  const Teuchos::RCP<StepperLeapfrogAppAction<Scalar> >& stepperLFAppAction)
  {
    this->setStepperType(        "Leapfrog");
    this->setUseFSAL(            useFSAL);
    this->setICConsistency(      ICConsistency);
    this->setICConsistencyCheck( ICConsistencyCheck);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    this->setObserver();
#endif
    this->setAppAction(stepperLFAppAction);
    if (appModel != Teuchos::null) {

      this->setModel(appModel);
      this->initialize();
    }
  }


#ifndef TEMPUS_HIDE_DEPRECATED_CODE
template<class Scalar>
void StepperLeapfrog<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (this->stepperObserver_ == Teuchos::null)
    this->stepperObserver_  =
      Teuchos::rcp(new StepperObserverComposite<Scalar>());

  if (obs == Teuchos::null) {
    if (stepperLFObserver_ == Teuchos::null)
      stepperLFObserver_ = Teuchos::rcp(new StepperLeapfrogObserver<Scalar>());
    if (this->stepperObserver_->getSize() == 0)
      this->stepperObserver_->addObserver(stepperLFObserver_);
  } else {
    stepperLFObserver_ =
      Teuchos::rcp_dynamic_cast<StepperLeapfrogObserver<Scalar> >(obs,true);
    this->stepperObserver_->addObserver(stepperLFObserver_);
  }

  this->isInitialized_ = false;
}
#endif

template<class Scalar>
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


template<class Scalar>
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
    Teuchos::OSTab ostab(out,1,"StepperLeapfrog::setInitialConditions()");
    *out << "Warning -- The First-Same-As-Last (FSAL) principle is not "
         << "used with Leapfrog because of the algorithm's prescribed "
         << "order of solution update. The default is to set useFSAL=false, "
         << "however useFSAL=true will also work but have no affect "
         << "(i.e., no-op).\n" << std::endl;
  }
}

template<class Scalar>
void StepperLeapfrog<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperLeapfrog::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperLeapfrog<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for Leapfrog.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    //this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar time = currentState->getTime();
    const Scalar dt   = workingState->getTimeStep();


    RCP<StepperLeapfrog<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
   
    // Perform half-step startup if working state is synced
    // (i.e., xDot and x are at the same time level).
    if (workingState->getIsSynced() == true) {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
      if (!Teuchos::is_null(stepperLFObserver_))
        stepperLFObserver_->observeBeforeXDotUpdateInitialize(
          solutionHistory, *this);
#endif
      stepperLFAppAction_->execute(solutionHistory, thisStepper,
        StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);
      // Half-step startup: xDot_{n+1/2} = xDot_n + 0.5*dt*xDotDot_n
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
        *(currentState->getXDot()),0.5*dt,*(currentState->getXDotDot()));
    }
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    if (!Teuchos::is_null(stepperLFObserver_))
      stepperLFObserver_->observeBeforeXUpdate(solutionHistory, *this);
#endif
    stepperLFAppAction_->execute(solutionHistory, thisStepper,
      StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::BEFORE_X_UPDATE);
    // x_{n+1} = x_n + dt*xDot_{n+1/2}
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
      *(currentState->getX()),dt,*(workingState->getXDot()));
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    if (!Teuchos::is_null(stepperLFObserver_))
      stepperLFObserver_->observeBeforeExplicit(solutionHistory, *this);
#endif
    stepperLFAppAction_->execute(solutionHistory, thisStepper,
      StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);
    auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

    // Evaluate xDotDot = f(x,t).
    this->evaluateExplicitODE(workingState->getXDotDot(),
                              workingState->getX(),
                              Teuchos::null, time+dt, p);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    if (!Teuchos::is_null(stepperLFObserver_))
      stepperLFObserver_->observeBeforeXDotUpdate(solutionHistory, *this);
#endif
    stepperLFAppAction_->execute(solutionHistory, thisStepper,
      StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION::BEFORE_XDOT_UPDATE);
    if (workingState->getOutput() == true) {
      // Half-step sync: xDot_{n+1} = xDot_{n+1/2} + 0.5*dt*xDotDot_{n+1}
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
        *(workingState->getXDot()),0.5*dt,*(workingState->getXDotDot()));
      workingState->setIsSynced(true);
    } else {
      // Full leapfrog step: xDot_{n+3/2} = xDot_{n+1/2} + dt*xDotDot_{n+1}
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
        *(workingState->getXDot()),dt,*(workingState->getXDotDot()));
      workingState->setIsSynced(false);
    }

    workingState->setSolutionStatus(Status::PASSED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    //this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
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
Teuchos::RCP<Tempus::StepperState<Scalar> > StepperLeapfrog<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperLeapfrog<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperExplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperLeapfrog ---\n";
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  out << "  stepperLFObserver_ = " << stepperLFObserver_ << std::endl;
#endif
  out << "  stepperLFAppAction_                = "
      << stepperLFAppAction_ << std::endl;
  out << "-----------------------" << std::endl;
}


template<class Scalar>
bool StepperLeapfrog<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperExplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  if (stepperLFObserver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Leapfrog observer is not set!\n";
  }
#endif
  if (stepperLFAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Leapfrog AppAction is not set!\n";
  }


  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperLeapfrog<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<std::string>("Initial Condition Consistency",
                       this->getICConsistencyDefault());
  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperLeapfrog_impl_hpp
