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
#include "Tempus_StepperBDF2ModifierDefault.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


template<class Scalar>
StepperBDF2<Scalar>::StepperBDF2()
{
  this->setStepperType(        "BDF2");
  this->setUseFSAL(            false);
  this->setICConsistency(      "None");
  this->setICConsistencyCheck( false);
  this->setZeroInitialGuess(   false);

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
  this->setStartUpStepper("DIRK 1 Stage Theta Method");
}

template<class Scalar>
StepperBDF2<Scalar>::StepperBDF2(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  const Teuchos::RCP<Stepper<Scalar> >& startUpStepper,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperBDF2AppAction<Scalar> >& stepperBDF2AppAction)
{
  this->setStepperType(        "BDF2");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setZeroInitialGuess(   zeroInitialGuess);

  this->setAppAction(stepperBDF2AppAction);
  this->setSolver(solver);
  this->setStartUpStepper(startUpStepper);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template<class Scalar>
void StepperBDF2<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperImplicit<Scalar>::setModel(appModel);
  if (startUpStepper_->getModel() == Teuchos::null) {
    startUpStepper_->setModel(appModel);
    startUpStepper_->initialize();
  }

  this->isInitialized_ = false;
}


/// Set the startup stepper to a default stepper.
template<class Scalar>
void StepperBDF2<Scalar>::setStartUpStepper(std::string startupStepperType)
{
  using Teuchos::RCP;
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
  if (this->wrapperModel_ != Teuchos::null &&
      this->wrapperModel_->getAppModel() != Teuchos::null) {
  startUpStepper_ =
    sf->createStepper(startupStepperType, this->wrapperModel_->getAppModel());
  } else {
    startUpStepper_ = sf->createStepper(startupStepperType);
  }

  this->isInitialized_ = false;
}


/// Set the start up stepper.
template<class Scalar>
void StepperBDF2<Scalar>::setStartUpStepper(
     Teuchos::RCP<Stepper<Scalar> > startUpStepper)
{
  startUpStepper_ = startUpStepper;

  if (this->wrapperModel_ != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(
                               this->wrapperModel_->getAppModel() == Teuchos::null, std::logic_error,
                               "Error - Can not set the startUpStepper to Teuchos::null.\n");

    if (startUpStepper->getModel() == Teuchos::null  &&
        this->wrapperModel_->getAppModel() != Teuchos::null) {
      startUpStepper_->setModel(this->wrapperModel_->getAppModel());
      startUpStepper_->initialize();
    }
  }

  this->isInitialized_ = false;
}

template<class Scalar>
void StepperBDF2<Scalar>::setAppAction(
  Teuchos::RCP<StepperBDF2AppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperBDF2AppAction_ =
      Teuchos::rcp(new StepperBDF2ModifierDefault<Scalar>());
  }
  else {
    stepperBDF2AppAction_ = appAction;
  }
}



template<class Scalar>
void StepperBDF2<Scalar>::initialize()
{
  StepperImplicit<Scalar>::initialize();
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
  else
    this->setStepperXDot(initialState->getXDot());

  StepperImplicit<Scalar>::setInitialConditions(solutionHistory);
}


template<class Scalar>
void StepperBDF2<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

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

    RCP<StepperBDF2<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperBDF2AppAction_->execute(solutionHistory, thisStepper,
                                   StepperBDF2AppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);


    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();

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
    stepperBDF2AppAction_->execute(solutionHistory, thisStepper,
      StepperBDF2AppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveImplicitODE(x, xDot, time, p);

    stepperBDF2AppAction_->execute(solutionHistory, thisStepper,
      StepperBDF2AppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);


    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(getOrder());
    workingState->computeNorms(currentState);
    stepperBDF2AppAction_->execute(solutionHistory, thisStepper,
                                   StepperBDF2AppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
  const Teuchos::EVerbosityLevel      verbLevel ) const
{
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperBDF2 ---\n";
  if (startUpStepper_ != Teuchos::null) {
    out << "  startup stepper type             = "
        << startUpStepper_->description() << std::endl;
  }
  out << "  startUpStepper_                  = "
      << startUpStepper_ << std::endl;
  out << "  startUpStepper_->isInitialized() = "
      << Teuchos::toString(startUpStepper_->isInitialized()) << std::endl;
  out << "  stepperBDF2AppAction_            = "
      << stepperBDF2AppAction_ << std::endl;
  out << "----------------------------" << std::endl;
  out << "  order_                           = " << order_ << std::endl;
  out << "-------------------" << std::endl;
}


template<class Scalar>
bool StepperBDF2<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if ( !this->startUpStepper_->isInitialized() ) {
    isValidSetup = false;
    out << "The startup stepper is not initialized!\n";
  }
  if (stepperBDF2AppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The BDF2 AppAction is not set!\n";
  }
  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperBDF2<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<bool>("Initial Condition Consistency Check", false);
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<bool>("Zero Initial Guess", false);
  pl->set<std::string>("Start Up Stepper Type", "DIRK 1 Stage Theta Method");
  Teuchos::RCP<Teuchos::ParameterList> solverPL = defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperBDF2_impl_hpp
