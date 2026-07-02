// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExponentialEuler_impl_hpp
#define Tempus_StepperExponentialEuler_impl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_StepperExponentialEulerModifierDefault.hpp"
#include "Tempus_StepperExponentialEuler_decl.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperExponentialEuler.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

namespace Tempus {


template<class Scalar>
StepperExponentialEuler<Scalar>::StepperExponentialEuler()
{
  this->setStepperName("Exponential Euler");
  this->setStepperType("Exponential Euler");
  this->setUseFSAL(false);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);

  this->setAppAction(Teuchos::null);
}


template<class Scalar>
StepperExponentialEuler<Scalar>::StepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  const Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> >& stepperEEAppAction)
{
  this->setStepperName("Exponential Euler");
  this->setStepperType("Exponential Euler");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);

  this->setAppAction(stepperEEAppAction);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template<class Scalar>
void StepperExponentialEuler<Scalar>::setAppAction(
  Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperEEAppAction_ =
      Teuchos::rcp(new StepperExponentialEulerModifierDefault<Scalar>());
  } else {
    stepperEEAppAction_ = appAction;
  }
  this->isInitialized_ = false;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperExponentialEuler::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperExponentialEuler<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for Exponential Euler.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<StepperExponentialEuler<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    Thyra::SolveStatus<Scalar> sStatus;

    RCP<SolutionState<Scalar>> workingState = solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar>> currentState = solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar>> xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar>> x = workingState->getX();

    // we always use the memory of workingState to compute the new solution
    // and initialize it to xOld at the beginning of takeStep
    // (this is needed for application models (appModel) that set boundary conditions during evalModel)
    Thyra::copy(*xOld, x.ptr());

    // from now, assume that x is equal to xOld (potentially modified for BC by evalModel)
    // only update x to the next step once all right-hand sides and matrices have been assembled

    // either get a fresh xDot vector from the workingState, or use the temporary one from getStepperXDot
    // save an RCP to the correct choice to xDot and this->getStepperXDot()
    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    RCP<Thyra::VectorBase<Scalar>> xDot = this->getStepperXDot();

    // xDot will be set to an all zero vector, to signal to the ModelEvaluator that we desire the implicit mode
    // it must remain zero until the last call to ModelEvaluator, until the end of this method
    Thyra::assign(xDot.ptr(), ST::zero());

    const Scalar time = workingState->getTime();
    const Scalar t0   = currentState->getTime();
    const Scalar dt = workingState->getTimeStep();
    auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));

    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXP);

    // setup system Jacobian (and mass) at the current time t0
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = this->createInArgsExponentialODE(x, xDot, t0, p);
    this->getPhiEvaluator()->setLinearizationPoint(inArgs, PhiInitialization::JACOBIAN_AND_MASS);

    // if requested, update any hyperpameters of the phiEvaluator
    const int adaptInterval = this->getAdaptPhiEvaluator();
    if ((adaptInterval > 0)
        && (workingState->getIndex() < 2 || workingState->getIndex() % adaptInterval == 0))
    {
      this->getPhiEvaluator()->adaptEvaluator();
    }

    // compute the right hand side Mf at x (which is still equal to xOld)
    // and potentially set the correct Dirichlet BC to x
    RCP<Thyra::VectorBase<Scalar>> Mf = Thyra::createMember(x->space());

    RCP<Thyra::VectorBase<Scalar>> xDotOld = currentState->getXDot();
    if (this->getUseFSAL() && currentState->getIsSynced() && xDotOld != Teuchos::null) {
      // Get Mf = -M*f(xOld, t0) from xDotOld = f(xOld, t0)
      assign(Mf.ptr(), ST::zero());
      this->getPhiEvaluator()->applyMass(Mf.ptr(), xDotOld);
      Thyra::scale(Scalar(-1.0), Mf.ptr());
    }
    else {
      // Evaluate Mf = -M*f(xOld, t0).
      // NB: this will also set the correct BC at time t0 to x
      this->evaluateExponentialODE(Mf, x, xDot, t0, p);

      // If the old state currentState has memory for xDot, save it now
      // TODO: should we do this, or not touch currentState, and rely on the user to set FSAL?
      if (xDotOld != Teuchos::null) {
        // Compute xDotOld = f(xOld, t0)
        assign(xDotOld.ptr(), ST::zero());
        this->getPhiEvaluator()->solveMass(xDotOld.ptr(), Mf);
        Thyra::scale(Scalar(-1.0), xDotOld.ptr());
        // x and xDot are now sync'ed or consistent at the same time level for the currentState.
        currentState->setIsSynced(true);
        // TODO: for time-dependent Dirichlet conditions, we should also overwrite x0 with x?
      }
    }

    // initialize vector for the update
    RCP<Thyra::VectorBase<Scalar>> vphi = Thyra::createMember(x->space());
    assign(vphi.ptr(), ST::zero());  // Must initialize to a guess before solve!

    // call the PhiEvaluator to compute vphi = - \phi( dt * J ) f = - \phi( dt * M^{-1} MJ ) M^{-1} Mf
    sStatus = this->getPhiEvaluator()->computePhi(vphi.ptr(), 1, dt, Mf);

    // TODO: make this configurable
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    // Teuchos::OSTab ostab(out, 1, "StepperEPI::takeStep");
    int current_iters = -1;
    if(!sStatus.extraParameters.is_null()) {
      current_iters = sStatus.extraParameters->get("Iteration Count", 0);
    }
    Scalar achieved_tol = sStatus.achievedTol;

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED && current_iters >=0) {
      *out << "Phi converged: iters: " << current_iters << " tol: " << achieved_tol << std::endl;
    }
    else {
      *out << sStatus.message << std::endl;
    }

    // compute the final update of x (which has the correct BC) with vphi
    // this can only happen at the end of the function,
    // after any calls to the ModelEvaluator that rely on x being equal to xOld in the non-BC degrees of freedom
    // TODO: should this only happen if PhiEvaluator converged?
    Thyra::Vp_StV(x.ptr(), Scalar(-1.0)*dt, *vphi);

    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::AFTER_EXP);

    // If FSAL is on, and there is memory for xDot, compute new xDot and sync
    if (this->getUseFSAL() && workingState->getXDot() != Teuchos::null) {
      // Get xDot = f(x, t) from Mf = -M*f(x, t)
      // reevaluate f at current time: this will also set the correct BC at time t to x
      this->evaluateExponentialODE(Mf, x, xDot, time, p);

      // in the (untested) case that the mass matrix depends on time, we should recompute it here
      // if the matrix is constant in time and cached properly, this should be a NOOP
      Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs
        = this->createInArgsExponentialODE(x, xDot, time, p);
      this->getPhiEvaluator()->setLinearizationPoint(inArgs, PhiInitialization::ONLY_MASS);

      // solve the mass matrix and scale to obtain final xDot
      this->getPhiEvaluator()->solveMass(xDot.ptr(), Mf);
      Thyra::scale(Scalar(-1.0), xDot.ptr());

      // xDot is a pointer to workingState->getXDot()
      workingState->setIsSynced(true);
    }
    else {
      // workingState->getXDot() is either a null pointer, or will be left as a zero vector
      workingState->setIsSynced(false);
    }

    // give the working state the status of the PhiEvaluator
    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperEEAppAction_->execute(solutionHistory, thisStepper,
      StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}


template<class Scalar>
void StepperExponentialEuler<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  StepperExponential<Scalar>::describe(out, verbLevel);
  out << "--- StepperExponentialEuler ---\n";
  out << "  stepperEEAppAction_                = "
      << stepperEEAppAction_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperExponentialEuler<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperExponential<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (stepperEEAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Backward Euler AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperExponentialEuler<Scalar>::getValidParameters() const
{
  auto pl = this->getValidParametersBasicExponential();
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

  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


} // namespace Tempus
#endif // Tempus_StepperExponentialEuler_impl_hpp
