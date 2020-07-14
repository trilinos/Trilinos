// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperDIRK_impl_hpp
#define Tempus_StepperDIRK_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;

template<class Scalar>
void StepperDIRK<Scalar>::setupDefault()
{
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());
  this->setUseEmbedded(        this->getUseEmbeddedDefault());
  this->setZeroInitialGuess(   false);

  this->setStageNumber(-1);

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  stepperObserver_ = Teuchos::rcp(new StepperRKObserverComposite<Scalar>());
#endif
  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}


#ifndef TEMPUS_HIDE_DEPRECATED_CODE
template<class Scalar>
void StepperDIRK<Scalar>::setup(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<StepperRKObserver<Scalar> >& obs,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool useEmbedded,
  bool zeroInitialGuess)
{
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setUseEmbedded(        useEmbedded);
  this->setZeroInitialGuess(   zeroInitialGuess);

  this->setStageNumber(-1);

  stepperObserver_ = Teuchos::rcp(new StepperRKObserverComposite<Scalar>());
  this->setObserver(obs);
  this->setAppAction(Teuchos::null);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}
#endif
template<class Scalar>
void StepperDIRK<Scalar>::setup(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool useEmbedded,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
{
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setUseEmbedded(        useEmbedded);
  this->setZeroInitialGuess(   zeroInitialGuess);

  this->setStageNumber(-1);

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  stepperObserver_ = Teuchos::rcp(new StepperRKObserverComposite<Scalar>());
  this->setObserver(Teuchos::null);
#endif
  this->setAppAction(stepperRKAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
void StepperDIRK<Scalar>::getValidParametersBasicDIRK(
  Teuchos::RCP<Teuchos::ParameterList> pl) const
{
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<bool>("Use Embedded", false,
    "'Whether to use Embedded Stepper (if available) or not\n"
    "  'true' - Stepper will compute embedded solution and is adaptive.\n"
    "  'false' - Stepper is not embedded(adaptive).\n");
  pl->set<std::string>("Description", this->getDescription());
  pl->set<std::string>("Solver Name", "Default Solver",
    "Name of ParameterList containing the solver specifications.");
  pl->set<bool>("Zero Initial Guess", false);
  pl->set<bool>("Reset Initial Guess", true);
  Teuchos::RCP<Teuchos::ParameterList> solverPL = defaultSolverParameters();
  pl->set("Default Solver", *solverPL);
}


#ifndef TEMPUS_HIDE_DEPRECATED_CODE
template<class Scalar>
void StepperDIRK<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (this->stepperObserver_ == Teuchos::null)
    this->stepperObserver_  =
      Teuchos::rcp(new StepperRKObserverComposite<Scalar>());

  if (( obs == Teuchos::null ) and (this->stepperObserver_->getSize() >0 ) )
    return;

  if (( obs == Teuchos::null ) and (this->stepperObserver_->getSize() == 0) )
     obs = Teuchos::rcp(new StepperRKObserver<Scalar>());

  // Check that this casts to prevent a runtime error if it doesn't
  if (Teuchos::rcp_dynamic_cast<StepperRKObserver<Scalar> > (obs) != Teuchos::null) {
    this->stepperObserver_->addObserver(
         Teuchos::rcp_dynamic_cast<StepperRKObserver<Scalar> > (obs, true) );
  } else {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,0,"setObserver");
    *out << "Tempus::StepperIMEX_RK::setObserver: Warning: An observer has been provided that";
    *out << " does not support Tempus::StepperRKObserver. This observer WILL NOT be added.";
    *out << " In the future, this will result in a runtime error!" << std::endl;
  }

  this->isInitialized_ = false;
}
#endif


template<class Scalar>
void StepperDIRK<Scalar>::initialize()
{
  // Initialize the stage vectors
  const int numStages = this->tableau_->numStages();
  this->stageX_    = this->wrapperModel_->getNominalValues().get_x()->clone_v();
  stageXDot_.resize(numStages);
  for (int i=0; i<numStages; ++i) {
    stageXDot_[i] = Thyra::createMember(this->wrapperModel_->get_f_space());
    assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }
  xTilde_    = Thyra::createMember(this->wrapperModel_->get_x_space());
  assign(xTilde_.ptr(),    Teuchos::ScalarTraits<Scalar>::zero());

  if (this->tableau_->isEmbedded() and this->getUseEmbedded()) {
    ee_    = Thyra::createMember(this->wrapperModel_->get_f_space());
    abs_u0 = Thyra::createMember(this->wrapperModel_->get_f_space());
    abs_u  = Thyra::createMember(this->wrapperModel_->get_f_space());
    sc     = Thyra::createMember(this->wrapperModel_->get_f_space());
  }

  StepperImplicit<Scalar>::initialize();
}


template<class Scalar>
void StepperDIRK<Scalar>::setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDot
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(stageXDot_.back());

  StepperImplicit<Scalar>::setInitialConditions(solutionHistory);
}


template<class Scalar>
void StepperDIRK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperDIRK::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperDIRK<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for DIRK.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
#endif
    RCP<StepperDIRK<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
      StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int numStages = this->tableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = this->tableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = this->tableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = this->tableau_->c();

    // Reset non-zero initial guess.
    if ( this->getResetInitialGuess() && (!this->getZeroInitialGuess()) )
      Thyra::assign(this->stageX_.ptr(), *(currentState->getX()));

    // Compute stage solutions
    bool pass = true;
    Thyra::SolveStatus<Scalar> sStatus;
    for (int i=0; i < numStages; ++i) {
      this->setStageNumber(i);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
      this->stepperObserver_->observeBeginStage(solutionHistory, *this);

      // ???: is it a good idea to leave this (no-op) here?
      this->stepperObserver_
          ->observeBeforeImplicitExplicitly(solutionHistory, *this);
#endif

      // Check if stageXDot_[i] is needed.
      bool isNeeded = false;
      for (int k=i+1; k<numStages; ++k) if (A(k,i) != 0.0) isNeeded = true;
      if (b(i) != 0.0) isNeeded = true;
      if (this->tableau_->isEmbedded() && this->getUseEmbedded() &&
          this->tableau_->bstar()(i) != 0.0)
        isNeeded = true;
      if (isNeeded == false) {
        assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
        continue;
      }

      Thyra::assign(xTilde_.ptr(), *(currentState->getX()));
      for (int j=0; j < i; ++j) {
        if (A(i,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          Thyra::Vp_StV(xTilde_.ptr(), dt*A(i,j), *(stageXDot_[j]));
        }
      }

      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STAGE);

      Scalar ts = time + c(i)*dt;
      if (A(i,i) == Teuchos::ScalarTraits<Scalar>::zero()) {
        // Explicit stage for the ImplicitODE_DAE
        if (i == 0 && this->getUseFSAL() &&
            workingState->getNConsecutiveFailures() == 0) {
          // Reuse last evaluation for first step
          RCP<Thyra::VectorBase<Scalar> > tmp = stageXDot_[0];
          stageXDot_[0] = stageXDot_.back();
          stageXDot_.back() = tmp;
        } else {
          // Calculate explicit stage
          typedef Thyra::ModelEvaluatorBase MEB;
          MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
          MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
          inArgs.set_x(xTilde_);
          if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(ts);
          if (inArgs.supports(MEB::IN_ARG_x_dot))
            inArgs.set_x_dot(Teuchos::null);
          outArgs.set_f(stageXDot_[i]);

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
          this->stepperObserver_->observeBeforeExplicit(solutionHistory, *this);
#endif
          this->wrapperModel_->getAppModel()->evalModel(inArgs,outArgs);
        }
      } else {
        // Implicit stage for the ImplicitODE_DAE
        const Scalar alpha = 1.0/(dt*A(i,i));
        const Scalar beta  = 1.0;

        // Setup TimeDerivative
        Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
          Teuchos::rcp(new StepperDIRKTimeDerivative<Scalar>(
            alpha,xTilde_.getConst()));

        auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
          timeDer, dt, alpha, beta, SOLVE_FOR_X, i));

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
        this->stepperObserver_->observeBeforeSolve(solutionHistory, *this);
#endif
        this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
          StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

        sStatus = this->solveImplicitODE(this->stageX_, stageXDot_[i], ts, p);

        if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass=false;

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
        this->stepperObserver_->observeAfterSolve(solutionHistory, *this);
#endif
        this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
          StepperRKAppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);

        timeDer->compute(this->stageX_, stageXDot_[i]);
      }
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
      this->stepperObserver_->observeEndStage(solutionHistory, *this);
#endif
      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);
      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::END_STAGE);
    }

    // Sum for solution: x_n = x_n-1 + Sum{ dt*b(i) * f(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0; i < numStages; ++i) {
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), dt*b(i), *(stageXDot_[i]));
      }
    }

    if (this->tableau_->isEmbedded() and this->getUseEmbedded()) {
      const Scalar tolRel = workingState->getTolRel();
      const Scalar tolAbs = workingState->getTolAbs();

      // just compute the error weight vector
      // (all that is needed is the error, and not the embedded solution)
      Teuchos::SerialDenseVector<int,Scalar> errWght = b ;
      errWght -= this->tableau_->bstar();

      // compute local truncation error estimate: | u^{n+1} - \hat{u}^{n+1} |
      // Sum for solution: ee_n = Sum{ (b(i) - bstar(i)) * dt*f(i) }
      assign(ee_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      for (int i=0; i < numStages; ++i) {
         if (errWght(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
            Thyra::Vp_StV(ee_.ptr(), dt*errWght(i), *(stageXDot_[i]));
         }
      }

      // compute: Atol + max(|u^n|, |u^{n+1}| ) * Rtol
      Thyra::abs( *(currentState->getX()), abs_u0.ptr());
      Thyra::abs( *(workingState->getX()), abs_u.ptr());
      Thyra::pair_wise_max_update(tolRel, *abs_u0, abs_u.ptr());
      Thyra::add_scalar(tolAbs, abs_u.ptr());

      // compute: || ee / sc ||
      assign(sc.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *ee_, *abs_u, sc.ptr());
      Scalar err = std::abs(Thyra::norm_inf(*sc));
      workingState->setErrorRel(err);

      // test if step should be rejected
      if (std::isinf(err) || std::isnan(err) || err > Teuchos::as<Scalar>(1.0))
        pass = false;
    }

    if (pass) workingState->setSolutionStatus(Status::PASSED);
    else      workingState->setSolutionStatus(Status::FAILED);

    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
#endif
    this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
      StepperRKAppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  // reset the stage number
  this->setStageNumber(-1);
  return;
}

/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperDIRK<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperDIRK<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperDIRK ---\n";
  out << "  tableau_            = " << this->tableau_ << std::endl;
  if (this->tableau_ != Teuchos::null) this->tableau_->describe(out, verbLevel);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  out << "  stepperObserver_    = " << stepperObserver_ << std::endl;
#endif
  out << "  stepperRKAppAction_= " << this->stepperRKAppAction_ << std::endl;
  out << "  xTilde_             = " << xTilde_ << std::endl;
  out << "  stageX_             = " << this->stageX_ << std::endl;
  out << "  stageXDot_.size()   = " << stageXDot_.size() << std::endl;
  const int numStages = stageXDot_.size();
  for (int i=0; i<numStages; ++i)
    out << "    stageXDot_["<<i<<"] = " << stageXDot_[i] << std::endl;
  out << "  useEmbedded_        = "
      << Teuchos::toString(useEmbedded_) << std::endl;
  out << "  ee_                 = " << ee_ << std::endl;
  out << "  abs_u0              = " << abs_u0 << std::endl;
  out << "  abs_u               = " << abs_u << std::endl;
  out << "  sc                  = " << sc << std::endl;
  out << "-------------------" << std::endl;
}


template<class Scalar>
bool StepperDIRK<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (this->tableau_ == Teuchos::null) {
    isValidSetup = false;
    out << "The tableau is not set!\n";
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  if (stepperObserver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The observer is not set!\n";
  }
#endif
  if (this->stepperRKAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperDIRK<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  this->getValidParametersBasicDIRK(pl);

  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperDIRK_impl_hpp
