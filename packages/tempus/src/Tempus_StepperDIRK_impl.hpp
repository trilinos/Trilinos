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

  stepperObserver_ = Teuchos::rcp(new StepperRKObserverComposite<Scalar>());
}


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

  stepperObserver_ = Teuchos::rcp(new StepperRKObserverComposite<Scalar>());
  this->setObserver(obs);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->setSolver(solver);
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

}


template<class Scalar>
void StepperDIRK<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( tableau_ == Teuchos::null, std::logic_error,
    "Error - Need to set the tableau, before calling "
    "StepperDIRK::initialize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperDIRK::initialize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION( this->solver_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the solver, setSolver(), before calling "
    "StepperDIRK::initialize()\n");

  this->setObserver();

  TEUCHOS_TEST_FOR_EXCEPTION( this->stepperObserver_ == Teuchos::null,
    std::logic_error,
    "Error - StepperRKObserver is null!\n");

  // Initialize the stage vectors
  const int numStages = tableau_->numStages();
  stageX_    = this->wrapperModel_->getNominalValues().get_x()->clone_v();
  stageXDot_.resize(numStages);
  for (int i=0; i<numStages; ++i) {
    stageXDot_[i] = Thyra::createMember(this->wrapperModel_->get_f_space());
    assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }
  xTilde_    = Thyra::createMember(this->wrapperModel_->get_x_space());
  assign(xTilde_.ptr(),    Teuchos::ScalarTraits<Scalar>::zero());

  if (tableau_->isEmbedded() and this->getUseEmbedded()) {
    ee_    = Thyra::createMember(this->wrapperModel_->get_f_space());
    abs_u0 = Thyra::createMember(this->wrapperModel_->get_f_space());
    abs_u  = Thyra::createMember(this->wrapperModel_->get_f_space());
    sc     = Thyra::createMember(this->wrapperModel_->get_f_space());
  }
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

    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int numStages = tableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = tableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = tableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = tableau_->c();

    // Reset non-zero initial guess.
    if ( this->getResetInitialGuess() && (!this->getZeroInitialGuess()) )
      Thyra::assign(stageX_.ptr(), *(currentState->getX()));

    // Compute stage solutions
    bool pass = true;
    Thyra::SolveStatus<Scalar> sStatus;
    for (int i=0; i < numStages; ++i) {
        this->stepperObserver_->observeBeginStage(solutionHistory, *this);

        // ???: is it a good idea to leave this (no-op) here?
        this->stepperObserver_
            ->observeBeforeImplicitExplicitly(solutionHistory, *this);

      if ( i == 0 && this->getUseFSAL() &&
           workingState->getNConsecutiveFailures() == 0 ) {

        RCP<Thyra::VectorBase<Scalar> > tmp = stageXDot_[0];
        stageXDot_[0] = stageXDot_.back();
        stageXDot_.back() = tmp;

      } else {

        Thyra::assign(xTilde_.ptr(), *(currentState->getX()));
        for (int j=0; j < i; ++j) {
          if (A(i,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
            Thyra::Vp_StV(xTilde_.ptr(), dt*A(i,j), *(stageXDot_[j]));
          }
        }

        Scalar ts = time + c(i)*dt;
        if (A(i,i) == Teuchos::ScalarTraits<Scalar>::zero()) {
          // Explicit stage for the ImplicitODE_DAE
          bool isNeeded = false;
          for (int k=i+1; k<numStages; ++k) if (A(k,i) != 0.0) isNeeded = true;
          if (b(i) != 0.0) isNeeded = true;
          if (isNeeded == false) {
            // stageXDot_[i] is not needed.
            assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
          } else {
            typedef Thyra::ModelEvaluatorBase MEB;
            MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
            MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
            inArgs.set_x(xTilde_);
            if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(ts);
            if (inArgs.supports(MEB::IN_ARG_x_dot))
              inArgs.set_x_dot(Teuchos::null);
            outArgs.set_f(stageXDot_[i]);

            this->stepperObserver_->observeBeforeExplicit(solutionHistory, *this);
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

          this->stepperObserver_->observeBeforeSolve(solutionHistory, *this);

          sStatus = this->solveImplicitODE(stageX_, stageXDot_[i], ts, p);

          if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass=false;

          this->stepperObserver_->observeAfterSolve(solutionHistory, *this);

          timeDer->compute(stageX_, stageXDot_[i]);
        }
      }

      this->stepperObserver_->observeEndStage(solutionHistory, *this);
    }

    // Sum for solution: x_n = x_n-1 + Sum{ dt*b(i) * f(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0; i < numStages; ++i) {
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), dt*b(i), *(stageXDot_[i]));
      }
    }

    if (tableau_->isEmbedded() and this->getUseEmbedded()) {
      RCP<SolutionStateMetaData<Scalar> > metaData=workingState->getMetaData();
      const Scalar tolAbs = metaData->getTolRel();
      const Scalar tolRel = metaData->getTolAbs();

      // just compute the error weight vector
      // (all that is needed is the error, and not the embedded solution)
      Teuchos::SerialDenseVector<int,Scalar> errWght = b ;
      errWght -= tableau_->bstar();

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
      metaData->setErrorRel(err);

      // test if step should be rejected
      if (std::isinf(err) || std::isnan(err) || err > Teuchos::as<Scalar>(1.0))
        pass = false;
    }

    if (pass) workingState->setSolutionStatus(Status::PASSED);
    else      workingState->setSolutionStatus(Status::FAILED);

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
   const Teuchos::EVerbosityLevel      /* verbLevel */) const
{
  out << this->getStepperType() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
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
