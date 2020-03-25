// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBackwardEuler_impl_hpp
#define Tempus_StepperBackwardEuler_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


template<class Scalar>
StepperBackwardEuler<Scalar>::StepperBackwardEuler()
{
  this->setStepperType(        "Backward Euler");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());
  this->setZeroInitialGuess(   false);

  this->setObserver();
  this->setDefaultSolver();
  this->setPredictor("None");
}


template<class Scalar>
StepperBackwardEuler<Scalar>::StepperBackwardEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<StepperObserver<Scalar> >& obs,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  const Teuchos::RCP<Stepper<Scalar> >& predictorStepper,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess)

{
  this->setStepperType(        "Backward Euler");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setZeroInitialGuess(   zeroInitialGuess);

  this->setObserver(obs);
  this->setSolver(solver);
  this->setPredictor(predictorStepper);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperImplicit<Scalar>::setModel(appModel);
  if (predictorStepper_ != Teuchos::null) {
    if (predictorStepper_->getModel() == Teuchos::null) {
      predictorStepper_->setModel(appModel);
      predictorStepper_->initialize();
    }
  }

  this->isInitialized_ = false;
}


/// Set the predictor to a Stepper with default settings.
template<class Scalar>
void StepperBackwardEuler<Scalar>::setPredictor(std::string predictorType)
{
  if (predictorType == "None") {
    predictorStepper_ = Teuchos::null;
    return;
  }

  using Teuchos::RCP;
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
  if (this->wrapperModel_ != Teuchos::null &&
      this->wrapperModel_->getAppModel() != Teuchos::null) {
    predictorStepper_ =
      sf->createStepper(predictorType, this->wrapperModel_->getAppModel());
  } else {
    predictorStepper_ = sf->createStepper(predictorType);
  }

  this->isInitialized_ = false;
}


/// Set the predictor.
template<class Scalar>
void StepperBackwardEuler<Scalar>::setPredictor(
  Teuchos::RCP<Stepper<Scalar> > predictorStepper)
{
  predictorStepper_ = predictorStepper;

  if (predictorStepper != Teuchos::null &&
      this->wrapperModel_ != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->wrapperModel_->getAppModel() == Teuchos::null, std::logic_error,
      "Error - Need to set the model, setModel(), before calling "
      "StepperBackwardEuler::setPredictor()\n");

    if (predictorStepper_->getModel() == Teuchos::null) {
      predictorStepper_->setModel(this->wrapperModel_->getAppModel());
      predictorStepper_->initialize();
    }
  }

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (this->stepperObserver_ == Teuchos::null) {
      stepperBEObserver_ =
        Teuchos::rcp(new StepperBackwardEulerObserver<Scalar>());
      this->stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >(stepperBEObserver_,true);
    }
  } else {
    this->stepperObserver_ = obs;
    stepperBEObserver_ =
      Teuchos::rcp_dynamic_cast<StepperBackwardEulerObserver<Scalar> >
        (this->stepperObserver_,true);
  }

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::setInitialConditions(
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
    Teuchos::OSTab ostab(out,1,"StepperBackwardEuler::setInitialConditions()");
    *out << "\nWarning -- The First-Step-As-Last (FSAL) principle is not "
         << "needed with Backward Euler.  The default is to set useFSAL=false, "
         << "however useFSAL=true will also work but have no affect "
         << "(i.e., no-op).\n" << std::endl;
  }
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperBackwardEuler::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperBackwardEuler<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for Backward Euler.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();

    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot(workingState);

    computePredictor(solutionHistory);
    if (workingState->getSolutionStatus() == Status::FAILED)
      return;

    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
      Teuchos::rcp(new StepperBackwardEulerTimeDerivative<Scalar>(
        Scalar(1.0)/dt,xOld));

    const Scalar alpha = Scalar(1.0)/dt;
    const Scalar beta  = Scalar(1.0);
    auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
      timeDer, dt, alpha, beta));

    if (!Teuchos::is_null(stepperBEObserver_))
      stepperBEObserver_->observeBeforeSolve(solutionHistory, *this);

    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveImplicitODE(x, xDot, time, p);

    if (!Teuchos::is_null(stepperBEObserver_))
      stepperBEObserver_->observeAfterSolve(solutionHistory, *this);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
  }
  return;
}

template<class Scalar>
void StepperBackwardEuler<Scalar>::computePredictor(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  if (predictorStepper_ == Teuchos::null) return;
  predictorStepper_->takeStep(solutionHistory);

  if (solutionHistory->getWorkingState()->getSolutionStatus()==Status::FAILED) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperBackwardEuler::computePredictor");
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
StepperBackwardEuler<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperBackwardEuler ---\n";
  if (predictorStepper_ != Teuchos::null) {
    out << "  predictor stepper type             = "
        << predictorStepper_->description() << std::endl;
  }
  out << "  predictorStepper_                  = "
      << predictorStepper_ << std::endl;
  out << "  predictorStepper_->isInitialized() = "
      << Teuchos::toString(predictorStepper_->isInitialized()) << std::endl;
  out << "  stepperBEObserver_             = "
      << stepperBEObserver_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperBackwardEuler<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (predictorStepper_ != Teuchos::null) {
    if ( !predictorStepper_->isInitialized() ) {
      isValidSetup = false;
      out << "The predictor stepper is not initialized!\n";
    }
  }

  if (stepperBEObserver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Backward Euler observer is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<bool>       ("Zero Initial Guess", false);
  pl->set<std::string>("Predictor Stepper Type", "None");
  Teuchos::RCP<Teuchos::ParameterList> solverPL = defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


template <class Scalar>
int
StepperBackwardEuler<Scalar>::stencilLength() const
{
  return 2;
}


template <class Scalar>
void
StepperBackwardEuler<Scalar>::computeStepResidual(
  Thyra::VectorBase<Scalar>& residual,
  const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
  const Teuchos::Array<Scalar>& t,
  const Thyra::VectorBase<Scalar>& p,
  const int param_index) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
  outArgs.set_f(Teuchos::rcpFromRef(residual));
  computeStepResidDerivImpl(outArgs, x, t, p, param_index);
}

template <class Scalar>
void
StepperBackwardEuler<Scalar>::computeStepJacobian(
  Thyra::LinearOpBase<Scalar>& jacobian,
  const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
  const Teuchos::Array<Scalar>& t,
  const Thyra::VectorBase<Scalar>& p,
  const int param_index,
  const int deriv_index) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
  TEUCHOS_ASSERT(outArgs.supports(MEB::OUT_ARG_W_op));
  outArgs.set_W_op(Teuchos::rcpFromRef(jacobian));
  computeStepResidDerivImpl(outArgs, x, t, p, param_index, deriv_index);
}

template <class Scalar>
void
StepperBackwardEuler<Scalar>::computeStepParamDeriv(
  Thyra::LinearOpBase<Scalar>& deriv,
  const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
  const Teuchos::Array<Scalar>& t,
  const Thyra::VectorBase<Scalar>& p,
  const int param_index) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
  TEUCHOS_ASSERT(outArgs.supports(MEB::OUT_ARG_DfDp, param_index).supports(MEB::DERIV_LINEAR_OP));
  outArgs.set_DfDp(param_index,
                   MEB::Derivative<Scalar>(Teuchos::rcpFromRef(deriv)));
  computeStepResidDerivImpl(outArgs, x, t, p, param_index);
}

template <class Scalar>
void
StepperBackwardEuler<Scalar>::computeStepSolver(
  Thyra::LinearOpWithSolveBase<Scalar>& jacobian_solver,
  const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
  const Teuchos::Array<Scalar>& t,
  const Thyra::VectorBase<Scalar>& p,
  const int param_index) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
  TEUCHOS_ASSERT(outArgs.supports(MEB::OUT_ARG_W));
  outArgs.set_W(Teuchos::rcpFromRef(jacobian_solver));
  computeStepResidDerivImpl(outArgs, x, t, p, param_index, 0);
}

template <class Scalar>
void
StepperBackwardEuler<Scalar>::computeStepResidDerivImpl(
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs,
  const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
  const Teuchos::Array<Scalar>& t,
  const Thyra::VectorBase<Scalar>& p,
  const int param_index,
  const int deriv_index) const
{
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  TEUCHOS_ASSERT(x.size() == 2);
  TEUCHOS_ASSERT(t.size() == 2);
  RCP<const Thyra::VectorBase<Scalar> > xn = x[0];
  RCP<const Thyra::VectorBase<Scalar> > xo = x[1];
  const Scalar tn = t[0];
  const Scalar to = t[1];
  const Scalar dt = tn-to;

  // compute x_dot
  RCP<Thyra::VectorBase<Scalar> > x_dot = xn->clone_v();
  Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
    Teuchos::rcp(new StepperBackwardEulerTimeDerivative<Scalar>(Scalar(1.0)/dt,xo));
  timeDer->compute(xn, x_dot);

  // evaluate model
  MEB::InArgs<Scalar> inArgs = this->wrapperModel_->getInArgs();
  inArgs.set_x(xn);
  if (inArgs.supports(MEB::IN_ARG_x_dot         )) inArgs.set_x_dot    (x_dot);
  if (inArgs.supports(MEB::IN_ARG_t             )) inArgs.set_t        (tn);
  if (inArgs.supports(MEB::IN_ARG_step_size     )) inArgs.set_step_size(dt);
  inArgs.set_p(param_index, Teuchos::rcpFromRef(p));
  TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_alpha));
  TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_beta));
  if (deriv_index == 0) {
    // df/dx_n = df/dx_dot * dx_dot/dx_n + df/dx_n = 1/dt*df/dx_dot + df/dx_n
    inArgs.set_alpha(Scalar(1.0)/dt);
    inArgs.set_beta(Scalar(1.0));
  }
  else if (deriv_index == 1) {
    // df/dx_{n-1} = df/dx_dot * dx_dot/dx_{n-1} = -1/dt*df/dx_dot
    inArgs.set_alpha(Scalar(-1.0)/dt);
    inArgs.set_beta(Scalar(0.0));
  }
  this->wrapperModel_->getAppModel()->evalModel(inArgs, outArgs);
}

} // namespace Tempus
#endif // Tempus_StepperBackwardEuler_impl_hpp
