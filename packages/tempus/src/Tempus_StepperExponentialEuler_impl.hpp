// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExponentialEuler_impl_hpp
#define Tempus_StepperExponentialEuler_impl_hpp

#include "Tempus_StepperExponentialEulerModifierDefault.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Tempus_StepperFactory.hpp"


namespace Tempus {


template<class Scalar>
StepperExponentialEuler<Scalar>::StepperExponentialEuler()
{
  this->setStepperName(        "Exponential Euler");
  this->setStepperType(        "Exponential Euler");
  this->setUseFSAL(            false);
  this->setICConsistency(      "None");
  this->setICConsistencyCheck( false);
  this->setZeroInitialGuess(   false);

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
  this->setPredictor();
}


template<class Scalar>
StepperExponentialEuler<Scalar>::StepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  const Teuchos::RCP<Stepper<Scalar> >& predictorStepper,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> >& stepperBEAppAction)
{
  this->setStepperName(        "Exponential Euler");
  this->setStepperType(        "Exponential Euler");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setZeroInitialGuess(   zeroInitialGuess);

  this->setAppAction(stepperBEAppAction);
  this->setSolver(solver);
  this->setPredictor(predictorStepper);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::setPredictor(std::string predictorType)
{
  if (predictorType == "None") {
    predictorStepper_ = Teuchos::null;
    return;
  }

  auto sf = Teuchos::rcp(new StepperFactory<Scalar>());
  if (this->wrapperModel_ != Teuchos::null &&
      this->wrapperModel_->getAppModel() != Teuchos::null) {
    predictorStepper_ = sf->createStepper(predictorType,
                                          this->wrapperModel_->getAppModel());
  } else {
    predictorStepper_ = sf->createStepper(predictorType);
  }

  this->isInitialized_ = false;
}


/// Set the predictor.
template<class Scalar>
void StepperExponentialEuler<Scalar>::setPredictor(
  Teuchos::RCP<Stepper<Scalar> > predictorStepper)
{
  predictorStepper_ = predictorStepper;
  if (predictorStepper_ == Teuchos::null) return;

  TEUCHOS_TEST_FOR_EXCEPTION(
    predictorStepper_->getModel() == Teuchos::null &&
    this->wrapperModel_->getAppModel() == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperExponentialEuler::setPredictor()\n");

  if (predictorStepper_->getModel() == Teuchos::null)
    predictorStepper_->setModel(this->wrapperModel_->getAppModel());
  predictorStepper_->initialize();

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::setAppAction(
  Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperBEAppAction_ =
      Teuchos::rcp(new StepperExponentialEulerModifierDefault<Scalar>());
  } else {
    stepperBEAppAction_ = appAction;
  }
  this->isInitialized_ = false;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperImplicit<Scalar>::setModel(appModel);

  if (predictorStepper_ != Teuchos::null) {
    // If predictor's model is not set, set it to the stepper model.
    if (predictorStepper_->getModel() == Teuchos::null) {
      predictorStepper_->setModel(appModel);
      predictorStepper_->initialize();
    }
  }

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::setInitialConditions(
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
void StepperExponentialEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperExponentialEuler::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperExponentialEuler<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for Backward Euler.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<StepperExponentialEuler<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperBEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();

    computePredictor(solutionHistory);
    if (workingState->getSolutionStatus() == Status::FAILED)
      return;

    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
      Teuchos::rcp(new StepperExponentialEulerTimeDerivative<Scalar>(
        Scalar(1.0)/dt,xOld));

    const Scalar alpha = Scalar(1.0)/dt;
    const Scalar beta  = Scalar(1.0);
    auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
      timeDer, dt, alpha, beta));

    stepperBEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveImplicitODE(x, xDot, time, p);

    stepperBEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperBEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}

template<class Scalar>
void StepperExponentialEuler<Scalar>::computePredictor(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  if (predictorStepper_ == Teuchos::null) return;
  predictorStepper_->takeStep(solutionHistory);

  if (solutionHistory->getWorkingState()->getSolutionStatus()==Status::FAILED) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperExponentialEuler::computePredictor");
    *out << "Warning - predictorStepper has failed." << std::endl;
  } else {
    // Reset status to WORKING since this is the predictor
    solutionHistory->getWorkingState()->setSolutionStatus(Status::WORKING);
  }
}


/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperExponentialEuler<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperExponentialEuler ---\n";
  out << "  predictorStepper_                  = "
      << predictorStepper_ << std::endl;
  if (predictorStepper_ != Teuchos::null) {
    out << "  predictorStepper_->isInitialized() = "
        << Teuchos::toString(predictorStepper_->isInitialized()) << std::endl;
    out << "  predictor stepper type             = "
        << predictorStepper_->description() << std::endl;
  }
  out << "  stepperBEAppAction_                = "
      << stepperBEAppAction_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperExponentialEuler<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (predictorStepper_ != Teuchos::null) {
    if ( !predictorStepper_->isInitialized() ) {
      isValidSetup = false;
      out << "The predictor stepper is not initialized!\n";
    }
  }

  if (stepperBEAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Backward Euler AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperExponentialEuler<Scalar>::getValidParameters() const
{
  auto pl = this->getValidParametersBasicImplicit();
  if (predictorStepper_ == Teuchos::null)
    pl->set("Predictor Stepper Type", "None");
  else
    pl->set("Predictor Stepper Type", predictorStepper_->getStepperType());
  return pl;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template<class Scalar>
Teuchos::RCP<StepperExponentialEuler<Scalar> >
createStepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponentialEuler<Scalar>());

  stepper->setStepperImplicitValues(pl);

  if (pl != Teuchos::null) {
    std::string predictorName =
      pl->get<std::string>("Predictor Stepper Type", "None");
    stepper->setPredictor(predictorName);
  }

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


} // namespace Tempus
#endif // Tempus_StepperExponentialEuler_impl_hpp
