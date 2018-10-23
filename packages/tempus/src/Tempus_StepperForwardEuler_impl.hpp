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

// StepperForwardEuler definitions:
template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validExplicitODE(appModel);
  this->appModel_ = appModel;

  this->inArgs_  = this->appModel_->getNominalValues();
  this->outArgs_ = this->appModel_->createOutArgs();
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
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
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >(stepperFEObserver_);
    }
  } else {
    this->stepperObserver_ = obs;
    stepperFEObserver_ =
      Teuchos::rcp_dynamic_cast<StepperForwardEulerObserver<Scalar> >
        (this->stepperObserver_);
  }
}

template<class Scalar>
void StepperForwardEuler<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->appModel_ == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperForwardEuler::initialize()\n");

  this->setParameterList(this->stepperPL_);
  this->setObserver();
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,
      "StepperForwardEuler::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x = initialState->getX();

  // Use x from inArgs as ICs, if needed.
  this->inArgs_ = this->appModel_->getNominalValues();
  if (x == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION( (x == Teuchos::null) &&
      (this->inArgs_.get_x() == Teuchos::null), std::logic_error,
      "Error - setInitialConditions() needs the ICs from the SolutionHistory\n"
      "        or getNominalValues()!\n");

    x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
          this->inArgs_.get_x());
    initialState->setX(x);
  }

  // Check if we need Stepper storage for xDot
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(initialState->getX()->clone_v());

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (initialState->getXDot() == Teuchos::null) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperForwardEuler::setInitialConditions()");
      *out << "Warning -- Requested IC consistency of 'None' but\n"
           << "           initialState does not have an xDot.\n"
           << "           Setting a 'Consistent' xDot!\n" << std::endl;
      this->evaluateExplicitODE(this->getStepperXDot(initialState), x,
                                initialState->getTime());
      initialState->setIsSynced(true);
    }
  }
  else if (icConsistency == "Zero")
    Thyra::assign(this->getStepperXDot(initialState).ptr(), Scalar(0.0));
  else if (icConsistency == "App") {
    auto xDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                  this->inArgs_.get_x_dot());
    TEUCHOS_TEST_FOR_EXCEPTION(xDot == Teuchos::null, std::logic_error,
      "Error - setInitialConditions() requested 'App' for IC consistency,\n"
      "        but 'App' returned a null pointer for xDot!\n");
    Thyra::assign(this->getStepperXDot(initialState).ptr(), *xDot);
  }
  else if (icConsistency == "Consistent") {
    // Evaluate xDot = f(x,t).
    this->evaluateExplicitODE(this->getStepperXDot(initialState), x,
                              initialState->getTime());

    // At this point, x, and xDot sync'ed or consistent
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
    auto xDot = this->getStepperXDot(initialState);
    auto f    = initialState->getX()->clone_v();
    this->evaluateExplicitODE(f, x, initialState->getTime());
    Thyra::Vp_StV(f.ptr(), Scalar(-1.0), *(xDot));
    Scalar reldiff = Thyra::norm(*f)/Thyra::norm(*xDot);

    Scalar eps = Scalar(100.0)*std::abs(Teuchos::ScalarTraits<Scalar>::eps());
    if (reldiff > eps) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperForwardEuler::setInitialConditions()");
      *out << "Warning -- Failed consistency check but continuing!\n"
         << "  ||xDot-f(x,t)||/||xDot|| > eps" << std::endl
         << "  ||xDot-f(x,t)||          = " << Thyra::norm(*f)    << std::endl
         << "  ||xDot||                 = " << Thyra::norm(*xDot) << std::endl
         << "  ||xDot-f(x,t)||/||xDot|| = " << reldiff            << std::endl
         << "                       eps = " << eps                << std::endl;
    }
  }
}

template<class Scalar>
void StepperForwardEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
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

    if ( !(this->getUseFSAL()) ) {
      // Need to compute XDotOld.
      if (!Teuchos::is_null(stepperFEObserver_))
        stepperFEObserver_->observeBeforeExplicit(solutionHistory, *this);

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, currentState->getX(),
                                currentState->getTime());

      // For UseFSAL=false, x and xDot are now sync'ed or consistent
      // at the same time level for the currentState.
      currentState->setIsSynced(true);
    }


    // Forward Euler update, x^n = x^{n-1} + dt^n * xDot^{n-1}
    const Scalar dt = workingState->getTimeStep();
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
      *(currentState->getX()),dt,*(xDot));


    xDot = this->getStepperXDot(workingState);

    if (this->getUseFSAL()) {
      // Get consistent xDot^n.
      if (!Teuchos::is_null(stepperFEObserver_))
        stepperFEObserver_->observeBeforeExplicit(solutionHistory, *this);

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, workingState->getX(),
                                workingState->getTime());

      // For UseFSAL=true, x and xDot are now sync'ed or consistent
      // for the workingState.
      workingState->setIsSynced(true);
    } else {
      assign(xDot.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
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
Teuchos::RCP<Tempus::StepperState<Scalar> > StepperForwardEuler<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperForwardEuler<Scalar>::description() const
{
  std::string name = "Forward Euler";
  return(name);
}


template<class Scalar>
void StepperForwardEuler<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "appModel_ = " << this->appModel_->description() << std::endl;
}


template <class Scalar>
void StepperForwardEuler<Scalar>::setParameterList(
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
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Forward Euler",
    std::logic_error,
       "Error - Stepper Type is not 'Forward Euler'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", "Forward Euler",
                       "'Stepper Type' must be 'Forward Euler'.");
  this->getValidParametersBasic(pl);
  pl->set<std::string>("Initial Condition Consistency",
                       "Consistent");   // Default for explicit ODEs.
  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getDefaultParameters() const
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
StepperForwardEuler<Scalar>::getNonconstParameterList()
{
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperForwardEuler<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperForwardEuler_impl_hpp
