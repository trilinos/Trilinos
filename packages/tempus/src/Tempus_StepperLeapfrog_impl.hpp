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
void StepperLeapfrog<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (this->stepperObserver_ == Teuchos::null) {
      stepperLFObserver_ =
        Teuchos::rcp(new StepperLeapfrogObserver<Scalar>());
      this->stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >(stepperLFObserver_);
     }
  } else {
    this->stepperObserver_ = obs;
    stepperLFObserver_ =
      Teuchos::rcp_dynamic_cast<StepperLeapfrogObserver<Scalar> >
        (this->stepperObserver_);
  }
}

template<class Scalar>
void StepperLeapfrog<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->appModel_ == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperLeapfrog::initialize()\n");

  this->setParameterList(this->stepperPL_);
  this->setObserver();
}

template<class Scalar>
void StepperLeapfrog<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperLeapfrog::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x    = initialState->getX();
  RCP<Thyra::VectorBase<Scalar> > xDot = initialState->getXDot();

  // If initialState has x and xDot set, treat them as the initial conditions.
  // Otherwise use the x and xDot from getNominalValues() as the ICs.
  TEUCHOS_TEST_FOR_EXCEPTION(
    !((x != Teuchos::null && xDot != Teuchos::null) ||
      (this->inArgs_.get_x() != Teuchos::null &&
       this->inArgs_.get_x_dot() != Teuchos::null)), std::logic_error,
    "Error - We need to set the initial conditions for x and xDot from\n"
    "        either initialState or appModel_->getNominalValues::InArgs\n"
    "        (but not from a mixture of the two).\n");

  this->inArgs_ = this->appModel_->getNominalValues();
  using Teuchos::rcp_const_cast;
  // Use the x and xDot from getNominalValues() as the ICs.
  if ( initialState->getX() == Teuchos::null ||
       initialState->getXDot() == Teuchos::null ) {
    TEUCHOS_TEST_FOR_EXCEPTION( (this->inArgs_.get_x() == Teuchos::null) ||
      (this->inArgs_.get_x_dot() == Teuchos::null), std::logic_error,
      "Error - setInitialConditions() needs the ICs from the initialState\n"
      "        or getNominalValues()!\n");
    x    =rcp_const_cast<Thyra::VectorBase<Scalar> >(this->inArgs_.get_x());
    initialState->setX(x);
    xDot =rcp_const_cast<Thyra::VectorBase<Scalar> >(this->inArgs_.get_x_dot());
    initialState->setXDot(xDot);
  }

  // Check if we need Stepper storage for xDotDot
  if (initialState->getXDotDot() == Teuchos::null)
    initialState->setXDotDot(initialState->getX()->clone_v());

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (initialState->getXDotDot() == Teuchos::null) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperForwardEuler::setInitialConditions()");
      *out << "Warning -- Requested IC consistency of 'None' but\n"
           << "           initialState does not have an xDotDot.\n"
           << "           Setting a 'Consistent' xDotDot!\n" << std::endl;
      this->evaluateExplicitODE(initialState->getXDotDot(), x,
                                Teuchos::null, initialState->getTime());
      initialState->setIsSynced(true);
    }
  }
  else if (icConsistency == "Zero")
    Thyra::assign(initialState->getXDotDot().ptr(), Scalar(0.0));
  else if (icConsistency == "App") {
    auto xDotDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                     this->inArgs_.get_x_dot_dot());
    TEUCHOS_TEST_FOR_EXCEPTION(xDotDot == Teuchos::null, std::logic_error,
      "Error - setInitialConditions() requested 'App' for IC consistency,\n"
      "        but 'App' returned a null pointer for xDotDot!\n");
    Thyra::assign(initialState->getXDotDot().ptr(), *xDotDot);
  }
  else if (icConsistency == "Consistent") {
    // Evaluate xDotDot = f(x,t).
    this->evaluateExplicitODE(initialState->getXDotDot(), x,
                              Teuchos::null, initialState->getTime());

    // At this point, x, xDot and xDotDot are sync'ed or consistent
    // at the same time level for the initialState.
    initialState->setIsSynced(true);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Error - setInitialConditions() invalid IC consistency, "
      << icConsistency << ".\n");
  }

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    auto xDotDot = initialState->getXDotDot();
    auto f       = initialState->getX()->clone_v();
    this->evaluateExplicitODE(f, x, Teuchos::null, initialState->getTime());
    Thyra::Vp_StV(f.ptr(), Scalar(-1.0), *(xDotDot));
    Scalar reldiff = Thyra::norm(*f)/Thyra::norm(*xDotDot);

    Scalar eps = Scalar(100.0)*std::abs(Teuchos::ScalarTraits<Scalar>::eps());
    if (reldiff > eps) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperForwardEuler::setInitialConditions()");
      *out << "Warning -- Failed consistency check but continuing!\n"
         << "  ||xDotDot-f(x,t)||/||xDotDot|| > eps" << std::endl
         << "  ||xDotDot-f(x,t)||             = " << Thyra::norm(*f)
         << std::endl
         << "  ||xDotDot||                    = " << Thyra::norm(*xDotDot)
         << std::endl
         << "  ||xDotDot-f(x,t)||/||xDotDot|| = " << reldiff << std::endl
         << "                             eps = " << eps     << std::endl;
    }
  }

  if (this->getUseFSAL()) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperLeapfrog::setInitialConditions()");
    *out << "Warning -- The First-Step-As-Last (FSAL) principle is not "
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

    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar time = currentState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    // Perform half-step startup if working state is synced
    // (i.e., xDot and x are at the same time level).
    if (workingState->getIsSynced() == true) {
      if (!Teuchos::is_null(stepperLFObserver_))
        stepperLFObserver_->observeBeforeXDotUpdateInitialize(
          solutionHistory, *this);
      // Half-step startup: xDot_{n+1/2} = xDot_n + 0.5*dt*xDotDot_n
      Thyra::V_VpStV(Teuchos::outArg(*(workingState->getXDot())),
        *(currentState->getXDot()),0.5*dt,*(currentState->getXDotDot()));
    }

    if (!Teuchos::is_null(stepperLFObserver_))
      stepperLFObserver_->observeBeforeXUpdate(solutionHistory, *this);
    // x_{n+1} = x_n + dt*xDot_{n+1/2}
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
      *(currentState->getX()),dt,*(workingState->getXDot()));

    if (!Teuchos::is_null(stepperLFObserver_))
      stepperLFObserver_->observeBeforeExplicit(solutionHistory, *this);

    // Evaluate xDotDot = f(x,t).
    this->evaluateExplicitODE(workingState->getXDotDot(),
                              workingState->getX(),
                              Teuchos::null, time+dt);

    if (!Teuchos::is_null(stepperLFObserver_))
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

    workingState->setSolutionStatus(Status::PASSED);
    workingState->setOrder(this->getOrder());
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
      << "appModel_ = " << this->appModel_->description() << std::endl;
}


template <class Scalar>
void StepperLeapfrog<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (this->stepperPL_ == Teuchos::null)
      this->stepperPL_ = this->getDefaultParameters();
  } else {
    this->stepperPL_ = pList;
  }
  this->stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType =
    this->stepperPL_->template get<std::string>("Stepper Type");
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
  pl->set<std::string>("Stepper Type", "Leapfrog",
                       "'Stepper Type' must be 'Leapfrog'.");
  this->getValidParametersBasic(pl);
  pl->set<std::string>("Initial Condition Consistency", "Consistent");
  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperLeapfrog<Scalar>::getDefaultParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::rcp_const_cast;

  RCP<ParameterList> pl =
    rcp_const_cast<ParameterList>(this->getValidParameters());

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperLeapfrog<Scalar>::getNonconstParameterList()
{
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperLeapfrog<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperLeapfrog_impl_hpp
