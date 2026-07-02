//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI_impl_hpp
#define Tempus_StepperEPI_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_StepperEPIModifierDefault.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperEPI.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

namespace Tempus {

template <class Scalar>
StepperEPI<Scalar>::StepperEPI()
{
  this->setStepperName("EPI");
  this->setStepperType("EPI");
  this->setUseFSAL(false);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setOrder(2.0);

  this->setAppAction(Teuchos::null);
}


template<class Scalar>
StepperEPI<Scalar>::StepperEPI(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
{
  this->setStepperName("EPI");
  this->setStepperType("EPI");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setOrder(2.0);

  this->setAppAction(stepperEPIAppAction);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template<class Scalar>
void StepperEPI<Scalar>::setAppAction(
  Teuchos::RCP<StepperEPIAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperEPIAppAction_ =
      Teuchos::rcp(new StepperEPIModifierDefault<Scalar>());
  } else {
    stepperEPIAppAction_ = appAction;
  }
  this->isInitialized_ = false;
}

template<class Scalar>
void StepperEPI<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperEPI::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION((solutionHistory->getNumStates() < 2 && order_ <= 2.0),
      std::logic_error,
      "Error - StepperEPI<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for EPI2.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\"\n"
      "  \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\" for EPI2.\n"
      "  or \"Storage Type\" = \"Unlimited\"\n");

    TEUCHOS_TEST_FOR_EXCEPTION((solutionHistory->getStorageLimit() < 3 && order_ > 2.0),
      std::logic_error,
      "Error - StepperEPI<Scalar>::takeStep(...)\n"
      "Need at least three SolutionStates for EPI3.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\"\n"
      "  \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"3\" for EPI3.\n"
      "  or \"Storage Type\" = \"Unlimited\"\n");

    // use EPI3 if order is 3 and we have enough history, fall back to EPI2 in first step
    bool useEPI3 = (solutionHistory->getNumStates() >= 3 && order_ > 2.0);

    Thyra::SolveStatus<Scalar> sStatus;

    RCP<StepperEPI<Scalar>> thisStepper = Teuchos::rcpFromRef(*this);
    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar>> workingState = solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar>> currentState = solutionHistory->getCurrentState();

    // TODO: remove this.
    const bool very_verbose = false;
    if (very_verbose) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      this->describe(*out, Teuchos::VERB_EXTREME);
      *out << "current: " << *currentState << std::endl;
      *out << "x: " << currentState->getX() << std::endl;
      *out << "xdot: " << currentState->getXDot() << std::endl;
      *out << "working: " << *workingState << std::endl;
      *out << "x: " << workingState->getX() << std::endl;
      *out << "xdot: " << workingState->getXDot() << std::endl;
      currentState->describe(*out, Teuchos::VERB_EXTREME);
      workingState->describe(*out, Teuchos::VERB_EXTREME);
      *out << "stepper_xdot: " << this->getStepperXDot() << std::endl;
    }

    RCP<const Thyra::VectorBase<Scalar>> xOld = currentState->getX();
    RCP<Thyra::VectorBase<Scalar>> x = workingState->getX();

    // we always use the memory of workingState to compute the new solution
    // and initialize it to xOld at the beginning of takeStep
    // (this is needed for application models (appModel) that set boundary conditions during evalModel)
    //
    // TODO: maybe we should allow the user to choose if not to overwrite it with xOld at the beginning,
    //       in case there are some meaningful changes to x done by an observer/modifier.
    Thyra::copy(*xOld, x.ptr());

    // from now, assume that x is equal to xOld (potentially modified for BC by evalModel)
    // only update x to the next step, once all right-hand sides and matrices have been assembled

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

    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
                                  StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXP);

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

    // First RHS evaluation
    // 
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

    Teuchos::ArrayRCP<RCP<const Thyra::VectorBase<Scalar>>> Mrhs_B(3);
    // Mrhs_B is default initialized with Teuchos::null
    Mrhs_B[1] = Mf;

    // if requested, compute the time derivative for the nonaotonomous correction
    RCP<Thyra::VectorBase<Scalar>> dt_Mf_deriv = Teuchos::null;
    if (this->getTemporalDerivative()) {
      dt_Mf_deriv = Thyra::createMember(Mf->space());
      this->computeTemporalFD(dt_Mf_deriv, x, t0, dt, Mf);

      // compute rhs for the phi_2 function.
      Mrhs_B[2] = dt_Mf_deriv;
    }

    // if requested, use the EPI3 3rd order update
    if (useEPI3) {
      // Retrieve x_{n-1} when available (used by EPI3 formula, added in next step).
      // Note this last solution is obtained with getStateTimeIndexNM2()
      // On the first step only two states exist, so the method runs as pure EPI2.
      RCP<SolutionState<Scalar>> NM2State = solutionHistory->getStateTimeIndexNM2();
      Scalar tOldOld = NM2State->getTime();
      RCP<const Thyra::VectorBase<Scalar>> xOldOld = NM2State->getX();
      RCP<const Thyra::VectorBase<Scalar>> xDotOldOld = NM2State->getXDot();

      // either retrieve from xDotOldOld, or compute remainder with new RHS eval
      RCP<Thyra::VectorBase<Scalar>> remf = Thyra::createMember(Mf->space());
      if (xDotOldOld != Teuchos::null && NM2State->getIsSynced()) {
        RCP<Thyra::VectorBase<Scalar>> MfOld = Thyra::createMember(Mf->space());
        assign(MfOld.ptr(), ST::zero());
        this->getPhiEvaluator()->applyMass(MfOld.ptr(), xDotOldOld);
        Thyra::scale(Scalar(-1.0), MfOld.ptr());
        this->computeRemf(remf, xOldOld, tOldOld, xOld, t0, dt, Mf, dt_Mf_deriv, MfOld);
      }
      else {
        this->computeRemf(remf, xOldOld, tOldOld, xOld, t0, dt, Mf, dt_Mf_deriv);
      }

      // Subtract (2/3)*R from phi_2 term
      Thyra::scale(Scalar(-2.0 / 3.0), remf.ptr());
      if (Mrhs_B[2] != Teuchos::null) {
        // add previous term, and set RCP to the sum
        Thyra::Vp_V(remf.ptr(), *Mrhs_B[2]);
      }
      Mrhs_B[2] = remf;
    }

    // initialize vector for the update
    RCP<Thyra::VectorBase<Scalar>> vphi = Thyra::createMember(x->space());
    assign(vphi.ptr(), ST::zero());  // Must initialize to a guess before solve!

    // if temporalDerivative is included or we use EPI3,
    // use the extension formula.
    if (this->getTemporalDerivative() | useEPI3) {
      sStatus = this->getPhiEvaluator()->computePhis(vphi.ptr(), dt, Mrhs_B());
    }
    else {
      sStatus = this->getPhiEvaluator()->computePhi(vphi.ptr(), 1, dt, Mf);
    }

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

    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::AFTER_EXP);

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
    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
      StepperEPIAppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}

template<class Scalar>
void StepperEPI<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  StepperExponential<Scalar>::describe(out, verbLevel);

  out << "--- StepperEPI ---\n";
  out << "  stepperEPIAppAction_                = "
      << stepperEPIAppAction_ << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperEPI<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if ( !StepperExponential<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (stepperEPIAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The EPI AppAction is not set!\n";
  }

  return isValidSetup;
}

} // namespace Tempus
#endif // Tempus_StepperEPI_impl_hpp
