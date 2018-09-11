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


// StepperBackwardEuler definitions:
template<class Scalar>
StepperBackwardEuler<Scalar>::StepperBackwardEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


/** \brief Set the predictor to a pre-defined predictor in the ParameterList.
 *  The predictor is set to predictorName sublist in the Stepper's
 *  ParameterList.  The predictorName sublist should already be defined
 *  in the Stepper's ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperBackwardEuler<Scalar>::setPredictor(std::string predictorName)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> predPL =
    Teuchos::sublist(this->stepperPL_, predictorName, true);
  this->stepperPL_->set("Predictor Name", predictorName);
  if (predictorStepper_ != Teuchos::null) predictorStepper_ = Teuchos::null;
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
}


/** \brief Set the predictor to the supplied Parameter sublist.
 *  This adds a new predictor Parameter sublist to the Stepper's ParameterList.
 *  If the predictor sublist is null, it tests if the predictor is set in
 *  the Stepper's ParameterList.
 */
template<class Scalar>
void StepperBackwardEuler<Scalar>::setPredictor(
  Teuchos::RCP<Teuchos::ParameterList> predPL)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  Teuchos::RCP<Teuchos::ParameterList> stepperPL = this->stepperPL_;
  std::string predictorName =
    stepperPL->get<std::string>("Predictor Name","None");
  if (is_null(predPL)) {
    if (predictorName != "None") {
      RCP<ParameterList> predPL =
        Teuchos::sublist(stepperPL, predictorName, true);
      RCP<StepperFactory<Scalar> > sf =
        Teuchos::rcp(new StepperFactory<Scalar>());
      predictorStepper_ =
        sf->createStepper(this->wrapperModel_->getAppModel(), predPL);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( predictorName == predPL->name(),
      std::logic_error,
         "Error - Trying to add a predictor that is already in ParameterList!\n"
      << "  Stepper Type = " << stepperPL->get<std::string>("Stepper Type")
      << "\n" << "  Predictor Name  = "<<predictorName<<"\n");
    predictorName = predPL->name();
    stepperPL->set("Predictor Name", predictorName);
    stepperPL->set(predictorName, *predPL);           // Add sublist
    if (predictorStepper_ != Teuchos::null) predictorStepper_ = Teuchos::null;
    RCP<StepperFactory<Scalar> > sf =
      Teuchos::rcp(new StepperFactory<Scalar>());
    predictorStepper_ =
      sf->createStepper(this->wrapperModel_->getAppModel(), predPL);
  }
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperObserver_ == Teuchos::null) {
      stepperBEObserver_ =
        Teuchos::rcp(new StepperBackwardEulerObserver<Scalar>());
      stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >(stepperBEObserver_);
    }
  } else {
    stepperObserver_ = obs;
    stepperBEObserver_ =
      Teuchos::rcp_dynamic_cast<StepperBackwardEulerObserver<Scalar> >
        (stepperObserver_);
  }
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->wrapperModel_ == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperBackwardEuler::initialize()\n");

  this->setParameterList(this->stepperPL_);
  this->setSolver();
  this->setPredictor();
  this->setObserver();
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
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

    stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar> > xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > xDot = workingState->getXDot();
    if (xDot == Teuchos::null) xDot = getXDotTemp(x);

    computePredictor(solutionHistory);
    if (workingState->getSolutionStatus() == Status::FAILED)
      return;

    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
      Teuchos::rcp(new StepperBackwardEulerTimeDerivative<Scalar>(1.0/dt,xOld));

    // Setup InArgs and OutArgs
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
    MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
    inArgs.set_x(x);
    if (inArgs.supports(MEB::IN_ARG_x_dot    )) inArgs.set_x_dot    (xDot);
    if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time);
    if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(dt);
    if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (1.0/dt);
    if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (1.0);

    this->wrapperModel_->setForSolve(timeDer, inArgs, outArgs);

    if (!Teuchos::is_null(stepperBEObserver_))
      stepperBEObserver_->observeBeforeSolve(solutionHistory, *this);

    const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(x);

    if (!Teuchos::is_null(stepperBEObserver_))
      stepperBEObserver_->observeAfterSolve(solutionHistory, *this);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED)
      workingState->setSolutionStatus(Status::PASSED);
    else
      workingState->setSolutionStatus(Status::FAILED);
    workingState->setOrder(this->getOrder());
    stepperObserver_->observeEndTakeStep(solutionHistory, *this);
  }
  return;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperBackwardEuler<Scalar>::
getXDotTemp(Teuchos::RCP<const Thyra::VectorBase<Scalar> > x) const
{
  if (xDotTemp_ == Teuchos::null) {
    xDotTemp_ = x->clone_v();
    Thyra::assign(xDotTemp_.ptr(), Scalar(0.0));
  }
  return xDotTemp_;
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
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperBackwardEuler<Scalar>::description() const
{
  std::string name = "Backward Euler";
  return(name);
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
}


template <class Scalar>
void StepperBackwardEuler<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  Teuchos::RCP<Teuchos::ParameterList> stepperPL = this->stepperPL_;
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (stepperPL == Teuchos::null) stepperPL = this->getDefaultParameters();
  } else {
    stepperPL = pList;
  }
  if (!(stepperPL->isParameter("Solver Name"))) {
    stepperPL->set<std::string>("Solver Name", "Default Solver");
    Teuchos::RCP<Teuchos::ParameterList> solverPL =
      this->defaultSolverParameters();
    stepperPL->set("Default Solver", *solverPL);
  }
  // Can not validate because of optional Parameters (e.g., Solver Name).
  // stepperPL->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = stepperPL->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Backward Euler",
    std::logic_error,
       "Error - Stepper Type is not 'Backward Euler'!\n"
    << "  Stepper Type = "<<stepperPL->get<std::string>("Stepper Type")<<"\n");

  this->stepperPL_ = stepperPL;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", this->description());
  pl->set("Zero Initial Guess", false);
  pl->set("Solver Name", "",
    "Name of ParameterList containing the solver specifications.");

  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getDefaultParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", this->description());
  pl->set<bool>       ("Zero Initial Guess", false);
  pl->set<std::string>("Solver Name", "Default Solver");

  RCP<ParameterList> solverPL = this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::getNonconstParameterList()
{
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBackwardEuler<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
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
  RCP<Thyra::VectorBase<Scalar> > x_dot = getXDotTemp(xn);
  Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
    Teuchos::rcp(new StepperBackwardEulerTimeDerivative<Scalar>(1.0/dt,xo));
  timeDer->compute(xn, x_dot);

  // evaluate model
  MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
  inArgs.set_x(xn);
  if (inArgs.supports(MEB::IN_ARG_x_dot         )) inArgs.set_x_dot    (x_dot);
  if (inArgs.supports(MEB::IN_ARG_t             )) inArgs.set_t        (tn);
  if (inArgs.supports(MEB::IN_ARG_step_size     )) inArgs.set_step_size(dt);
  inArgs.set_p(param_index, Teuchos::rcpFromRef(p));
  TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_alpha));
  TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_beta));
  if (deriv_index == 0) {
    // df/dx_n = df/dx_dot * dx_dot/dx_n + df/dx_n = 1/dt*df/dx_dot + df/dx_n
    inArgs.set_alpha(1.0/dt);
    inArgs.set_beta(1.0);
  }
  else if (deriv_index == 1) {
    // df/dx_{n-1} = df/dx_dot * dx_dot/dx_{n-1} = -1/dt*df/dx_dot
    inArgs.set_alpha(-1.0/dt);
    inArgs.set_beta(0.0);
  }
  this->wrapperModel_->getAppModel()->evalModel(inArgs, outArgs);
}

} // namespace Tempus
#endif // Tempus_StepperBackwardEuler_impl_hpp
