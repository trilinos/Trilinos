// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBDF2_impl_hpp
#define Tempus_StepperBDF2_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


template<class Scalar>
StepperBDF2<Scalar>::StepperBDF2()
{
  this->setStepperType(        "BDF2");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());
  this->setZeroInitialGuess(   false);

  this->setObserver();
}


template<class Scalar>
StepperBDF2<Scalar>::StepperBDF2(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<StepperObserver<Scalar> >& obs,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  const Teuchos::RCP<Stepper<Scalar> >& startUpStepper,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess)

{
  this->setStepperType(        "BDF2");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setZeroInitialGuess(   zeroInitialGuess);

  this->setObserver(obs);

  if (appModel != Teuchos::null) {

    this->setModel(appModel);
    this->setSolver(solver);
    this->setStartUpStepper(startUpStepper);
    this->initialize();
  }
}


/// Set the startup stepper to a default stepper.
template<class Scalar>
void StepperBDF2<Scalar>::setStartUpStepper(std::string startupStepperType)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->wrapperModel_->getAppModel() == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperBDF2::setStartUpStepper()\n");

  using Teuchos::RCP;
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
  startUpStepper_ =
    sf->createStepper(startupStepperType, this->wrapperModel_->getAppModel());
}


/// Set the start up stepper.
template<class Scalar>
void StepperBDF2<Scalar>::setStartUpStepper(
  Teuchos::RCP<Stepper<Scalar> > startUpStepper)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->wrapperModel_->getAppModel() == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperBDF2::setStartUpStepper()\n");

  startUpStepper_ = startUpStepper;
  startUpStepper_->setModel(this->wrapperModel_->getAppModel());
  startUpStepper_->initialize();
}


template<class Scalar>
void StepperBDF2<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{

  if (this->stepperObserver_ == Teuchos::null)
    this->stepperObserver_  =
      Teuchos::rcp(new StepperObserverComposite<Scalar>());

  if (( obs == Teuchos::null ) and (this->stepperObserver_->getSize() == 0) )
    obs = Teuchos::rcp(new StepperBDF2Observer<Scalar>());

  this->stepperObserver_->addObserver(
      Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> > (obs, true) );

}


template<class Scalar>
void StepperBDF2<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperBDF2::initialize()\n");

  this->setObserver();
  order_ = Scalar(2.0);
}


template<class Scalar>
void StepperBDF2<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDot
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(initialState->getX()->clone_v());

  StepperImplicit<Scalar>::setInitialConditions(solutionHistory);

  if (this->getUseFSAL()) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperBDF2::setInitialConditions()");
    *out << "\nWarning -- The First-Step-As-Last (FSAL) principle is not "
         << "needed with Backward Euler.  The default is to set useFSAL=false, "
         << "however useFSAL=true will also work but have no affect "
         << "(i.e., no-op).\n" << std::endl;
  }
}


template<class Scalar>
void StepperBDF2<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperBDF2::takeStep()");
  {
    int numStates = solutionHistory->getNumStates();

    RCP<Thyra::VectorBase<Scalar> > xOld;
    RCP<Thyra::VectorBase<Scalar> > xOldOld;

    // If there are less than 3 states (e.g., first time step), call
    // startup stepper and return.
    if (numStates < 3) {
      computeStartUp(solutionHistory);
      return;
    }
    TEUCHOS_TEST_FOR_EXCEPTION( (numStates < 3), std::logic_error,
      "Error in Tempus::StepperBDF2::takeStep(): numStates after \n"
        << "startup stepper must be at least 3, whereas numStates = "
        << numStates <<"!\n" << "If running with Storage Type = Static, "
        << "make sure Storage Limit > 2.\n");

    //IKT, FIXME: add error checking regarding states being consecutive and
    //whether interpolated states are OK to use.

    //this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot(workingState);

    //get time, dt and dtOld
    const Scalar time  = workingState->getTime();
    const Scalar dt    = workingState->getTimeStep();
    const Scalar dtOld = currentState->getTimeStep();

    xOld    = solutionHistory->getStateTimeIndexNM1()->getX();
    xOldOld = solutionHistory->getStateTimeIndexNM2()->getX();
    order_ = Scalar(2.0);

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
      Teuchos::rcp(new StepperBDF2TimeDerivative<Scalar>(dt, dtOld, xOld, xOldOld));

    const Scalar alpha = getAlpha(dt, dtOld);
    const Scalar beta  = getBeta (dt);

    auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
      timeDer, dt, alpha, beta));

    if (!Teuchos::is_null(stepperBDF2Observer_))
      stepperBDF2Observer_->observeBeforeSolve(solutionHistory, *this);

    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveImplicitODE(x, xDot, time, p);

    if (!Teuchos::is_null(stepperBDF2Observer_))
      stepperBDF2Observer_->observeAfterSolve(solutionHistory, *this);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(getOrder());
    workingState->computeNorms(currentState);
    //this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
  }
  return;
}

template<class Scalar>
void StepperBDF2<Scalar>::computeStartUp(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperBDF2::computeStartUp()");
  *out << "Warning -- Taking a startup step for BDF2 using '"
       << startUpStepper_->getStepperType()<<"'!" << std::endl;

  //Take one step using startUpStepper_
  startUpStepper_->takeStep(solutionHistory);

  order_ = startUpStepper_->getOrder();
}

/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperBDF2<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperBDF2<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      /* verbLevel */) const
{
  out << this->getStepperType() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperBDF2<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<bool>("Initial Condition Consistency Check",
                this->getICConsistencyCheckDefault());
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<bool>("Zero Initial Guess", false);
  pl->set<std::string>("Start Up Stepper Type", "DIRK 1 Stage Theta Method");
  Teuchos::RCP<Teuchos::ParameterList> solverPL = defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperBDF2_impl_hpp
