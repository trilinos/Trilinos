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
#include "Teuchos_TimeMonitor.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Tempus {

// StepperLeapfrog definitions:
template<class Scalar>
StepperLeapfrog<Scalar>::StepperLeapfrog(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}

template<class Scalar>
void StepperLeapfrog<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validExplicitODE(appModel);
  appModel_ = appModel;

  inArgs_  = appModel_->getNominalValues();
  outArgs_ = appModel_->createOutArgs();
}

template<class Scalar>
void StepperLeapfrog<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}

template<class Scalar>
void StepperLeapfrog<Scalar>::setSolver(std::string solverName)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperLeapfrog::setSolver()");
  *out << "Warning -- No solver to set for StepperLeapfrog "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperLeapfrog<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperLeapfrog::setSolver()");
  *out << "Warning -- No solver to set for StepperLeapfrog "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperLeapfrog<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperLeapfrog::setSolver()");
  *out << "Warning -- No solver to set for StepperLeapfrog "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperLeapfrog<Scalar>::setObserver(
  Teuchos::RCP<StepperLeapfrogObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperLFObserver_ == Teuchos::null) {
      stepperLFObserver_ =
        Teuchos::rcp(new StepperLeapfrogObserver<Scalar>());
    }
  } else {
    stepperLFObserver_ = obs;
  }
}

template<class Scalar>
void StepperLeapfrog<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperLeapfrog::takeStep()");
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    stepperLFObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    // Perform half-step startup if working state is synced
    // (i.e., xDot and x are at the same time level).
    if (workingState->getIsSynced() == true) {
      if (getIsXDotXDotInitialized() == false) {
        inArgs_.set_x(currentState->getX());
        if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);

        // For model evaluators whose state function f(x, x_dot, x_dot_dot, t)
        // describes an implicit ODE, and which accept the optional input
        // arguments, x_dot and x_dot_dot, make sure they are set to null in
        // order to request the evaluation of a state function corresponding
        // to the explicit ODE formulation x_dot_dot = f(x, t) for leapfrog.
        if (inArgs_.supports(MEB::IN_ARG_x_dot))
          inArgs_.set_x_dot(Teuchos::null);
        if (inArgs_.supports(MEB::IN_ARG_x_dot_dot))
          inArgs_.set_x_dot_dot(Teuchos::null);
        outArgs_.set_f(currentState->getXDotDot());

        stepperLFObserver_->observeBeforeExplicitInitialize(
          solutionHistory, *this);
        appModel_->evalModel(inArgs_,outArgs_);
        setIsXDotXDotInitialized(true);
      }

      stepperLFObserver_->observeBeforeXDotUpdateInitialize(
          solutionHistory, *this);
      // Half-step startup: xDot_{n+1/2} = xDot_n + 0.5*dt*xDotDot_n
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
        *(currentState->getXDot()),0.5*dt,*(currentState->getXDotDot()));
    }

    stepperLFObserver_->observeBeforeXUpdate(solutionHistory, *this);
    // x_{n+1} = x_n + dt*xDot_{n+1/2}
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
      *(currentState->getX()),dt,*(workingState->getXDot()));

    inArgs_.set_x(workingState->getX());
    if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time+dt);

    // For model evaluators whose state function f(x, x_dot, x_dot_dot, t)
    // describes an implicit ODE, and which accept the optional input
    // arguments, x_dot and x_dot_dot, make sure they are set to null in
    // order to request the evaluation of a state function corresponding
    // to the explicit ODE formulation x_dot_dot = f(x, t) for leapfrog.
    if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
      if (inArgs_.supports(MEB::IN_ARG_x_dot_dot))
        inArgs_.set_x_dot_dot(Teuchos::null);
    outArgs_.set_f(workingState->getXDotDot());

    stepperLFObserver_->observeBeforeExplicit(solutionHistory, *this);
    appModel_->evalModel(inArgs_,outArgs_);

    stepperLFObserver_->observeBeforeXDotUpdate(solutionHistory, *this);
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

    workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    workingState->setOrder(this->getOrder());
    stepperLFObserver_->observeEndTakeStep(solutionHistory, *this);
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
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperLeapfrog<Scalar>::description() const
{
  std::string name = "Leapfrog";
  return(name);
}


template<class Scalar>
void StepperLeapfrog<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "appModel_ = " << appModel_->description() << std::endl;
}


template <class Scalar>
void StepperLeapfrog<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) stepperPL_ = this->getDefaultParameters();
  else stepperPL_ = pList;
  stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = stepperPL_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Leapfrog",
    std::logic_error,
       "Error - Stepper Type is not 'Leapfrog'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperLeapfrog<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", "Leapfrog",
          "'Stepper Type' must be 'Leapfrog'.");
  pl->set<bool>("Is xDotDot Initialized", 0,
    "At the beginning of an integration, the solution may or may not "
    "be initialized.  If false, the Leapfrog steppers will initialize "
    "xDotDot during the first timestep.");

  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperLeapfrog<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  *pl = *(this->getValidParameters());
  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperLeapfrog<Scalar>::getNonconstParameterList()
{
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperLeapfrog<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperLeapfrog_impl_hpp
