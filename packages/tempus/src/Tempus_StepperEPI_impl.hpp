// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPI_impl_hpp
#define Tempus_StepperEPI_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_StepperEPIModifierDefault.hpp"
#include "Tempus_StepperEPI_decl.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Tempus_StepperFactory.hpp"

// TODO: have to include this header to get LSP to work.
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
  this->setZeroInitialGuess(false);
  this->setOrder(2.0);

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}


template<class Scalar>
StepperEPI<Scalar>::StepperEPI(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
{
  this->setStepperName("EPI");
  this->setStepperType("EPI");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);
  this->setOrder(2.0);

  this->setAppAction(stepperEPIAppAction);
  this->setSolver(solver);

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
void StepperEPI<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperImplicit<Scalar>::setModel(appModel);

  phiEvaluator_->setModel(appModel);
  phiEvaluator_->initialize();

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperEPI<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using magScalar = typename ST::magnitudeType;
  typedef Teuchos::ScalarTraits<magScalar> STM;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(
      numStates < 1, std::logic_error,
      "Error - setInitialConditions() needs at least one SolutionState\n"
      "        to set the initial condition.  Number of States = "
          << numStates);
  TEUCHOS_TEST_FOR_EXCEPTION(
      (this->getOrderODE() == SECOND_ORDER_ODE), std::logic_error,
      "Error - StepperEPI does not support SECOND_ORDER_ODE.\n");

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  RCP<Thyra::VectorBase<Scalar> > x        = initialState->getX();
  RCP<Thyra::VectorBase<Scalar> > xDot     = initialState->getXDot();

  auto inArgs = this->getModel()->getNominalValues();
  // Use x from inArgs as ICs if null.
  if (x == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        (x == Teuchos::null) && (inArgs.get_x() == Teuchos::null),
        std::logic_error,
        "Error - setInitialConditions needs the ICs from the "
        "SolutionHistory\n"
        "        or getNominalValues()!\n");

    x = inArgs.get_x()->clone_v();
    initialState->setX(x);
  }

  // check and remember if the solutionHistory stores xdot
  bool xDotHistoryStored = !(xDot == Teuchos::null);

  // Check if we need to create local Stepper storage for xDot
  if (!xDotHistoryStored) {
    xDot = Thyra::createMember(x->space());
    Thyra::assign(xDot.ptr(), ST::zero());
    this->setStepperXDot(xDot);
    }
  else
    this->setStepperXDot(xDot);

  magScalar reldiff = STM::zero();
  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    // Mark xDot as not synced
    initialState->setIsSynced(false);
  }
  else if (icConsistency == "Zero") {
    // assign zero to xDot, even if something was provided in currentState
    Thyra::assign(xDot.ptr(), ST::zero());

    // TODO: Mark xDot as not synced or synced?
    //initialState->setIsSynced(false);
  }
  else if (icConsistency == "App") {
    auto x_dot = inArgs.get_x_dot();
    TEUCHOS_TEST_FOR_EXCEPTION(
        x_dot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDot!\n");
    Thyra::assign(xDot.ptr(), *x_dot);

    // TODO: Mark xDot as not synced or synced?
    //initialState->setIsSynced(false);
  }
  else if (icConsistency == "Consistent") {
    TEUCHOS_TEST_FOR_EXCEPTION(
        !xDotHistoryStored, std::logic_error,
        "Error - icConsistency = 'Consistent' requested, but xDot not provided in history.\n");
    // TODO: since the user requested a "Consistent" initial condition,
    // should we make sure that xDot will be saved to the history, even if Teuchos::zero initially?
    // initialState->setXDot(xDot);
    // xDotHistoryStored = true;

    // xDot will be set to an all zero vector, to signal to the ModelEvaluator that we desire the implicit mode
    // it must remain zero until the last call to ModelEvaluator, until the end of this method
    Thyra::assign(xDot.ptr(), ST::zero());

    const Scalar t0 = initialState->getTime();
    const Scalar dt = initialState->getTimeStep();

    auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
    // compute the right hand side for x
    // and potentially set the correct Dirichlet BC to x
    RCP<Thyra::VectorBase<Scalar>> Mf = Thyra::createMember(x->space());
    this->evaluateExponentialODE(Mf, x, xDot, t0, p);

    // setup system Jacobian computation
    // TODO: only using the mass matrix here, fix this.
    phiEvaluator_->setLinearizationPoint(inArgs);

    // overwrite xDot with f = (-M) \ Mf
    Thyra::scale(Scalar(-1.0), Mf.ptr());
    phiEvaluator_->solveMass(xDot.ptr(), Mf);
    reldiff = STM::zero();  // TODO: make solveMass return the tolerance.

    // At this point, x and xDot are sync'ed or consistent
    // at the same time level for the initialState.
    initialState->setIsSynced(true);

    // TODO: we do not track potential failures in mass solve
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //    sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED,
    //    std::logic_error,
    //    "Error - Solver failed while determining the initial conditions.\n"
    //    "        Solver status is "
    //        << Thyra::toString(sStatus.solveStatus) << ".\n");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error - setInitialConditions() invalid IC consistency, "
            << icConsistency << ".\n");
  }

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        !xDotHistoryStored, std::logic_error,
        "Error - ICConsistenceCheck requested, but xDot not provided in history.\n");

    RCP<Thyra::VectorBase<Scalar>> xDot_temp = Thyra::createMember(x->space());
    Thyra::assign(xDot_temp.ptr(), ST::zero());

    const Scalar t0 = initialState->getTime();
    const Scalar dt = initialState->getTimeStep();

    auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
    // compute the right hand side for x
    // and potentially set the correct Dirichlet BC to x
    RCP<Thyra::VectorBase<Scalar>> Mf = Thyra::createMember(x->space());
    this->evaluateExponentialODE(Mf, x, xDot_temp, t0, p);

    // setup system Jacobian computation
    phiEvaluator_->setLinearizationPoint(inArgs);

    // overwrite xDot_temp with Mf - (-M) * xDot
    phiEvaluator_->applyMass(xDot_temp.ptr(), xDot);
    Thyra::Vp_V(xDot_temp.ptr(), *Mf);

    magScalar normX = Thyra::norm(*x);
    if (normX == STM::zero())
      reldiff = Thyra::norm(*xDot_temp);
    else
      reldiff = Thyra::norm(*xDot_temp) / normX;

    magScalar eps = magScalar(100.0) * STM::eps();
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out, 1, "StepperImplicit::setInitialConditions()");
    if (reldiff < eps) {
      *out << "\n---------------------------------------------------\n"
           << "Info -- Stepper = " << this->getStepperType() << "\n"
           << "  Initial condition PASSED consistency check!\n"
           << "  (||f(x,xDot,t)||/||x|| = " << reldiff << ") < "
           << "(eps = " << eps << ")" << std::endl
           << "---------------------------------------------------\n"
           << std::endl;
    }
    else {
      *out << "\n---------------------------------------------------\n"
           << "Info -- Stepper = " << this->getStepperType() << "\n"
           << "  Initial condition FAILED consistency check but continuing!\n"
           << "  (||f(x,xDot,t)||/||x|| = " << reldiff << ") > "
           << "(eps = " << eps << ")" << std::endl
           << "  ||f(x,xDot,t)|| = " << Thyra::norm(*xDot_temp) << std::endl
           << "  ||x||           = " << Thyra::norm(*x) << std::endl
           << "---------------------------------------------------\n"
           << std::endl;
    }
  }
}


template <class Scalar>
void StepperEPI<Scalar>::evaluateExponentialODE(
    Teuchos::RCP<Thyra::VectorBase<Scalar> >& f,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
    const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
  MEB::InArgs<Scalar> inArgs = appModel->createInArgs();
  inArgs.set_x(x);
  TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_x_dot))
  inArgs.set_x_dot(xDot);

  // TODO: we could also rely on the user (takeStep) to set this, but double setting may be more robust
  Thyra::assign(xDot.ptr(), ST::zero());

  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time);
  if (inArgs.supports(MEB::IN_ARG_step_size))
    inArgs.set_step_size(p->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha)) inArgs.set_alpha(Scalar(0.0));
  if (inArgs.supports(MEB::IN_ARG_beta)) inArgs.set_beta(Scalar(1.0));
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(p->stageNumber_);

  MEB::OutArgs<Scalar> outArgs = appModel->createOutArgs();
  outArgs.set_f(f);

  appModel->evalModel(inArgs, outArgs);
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

    RCP<SolutionState<Scalar>> workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar>> currentState=solutionHistory->getCurrentState();

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

    stepperEPIAppAction_->execute(solutionHistory, thisStepper,
                                  StepperEPIAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXP);

    // Use the appModel to configure Jacobian and Mass
    RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = appModel->createInArgs();
    // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
    inArgs.set_x(x);  // we need to set x here, sice this vector may get boundary conditions set
    inArgs.set_t(t0);
    // set x_dot == 0 to signal to some model evaluators that we want the implicit version
    inArgs.set_x_dot(xDot);  // we need to make sure that this vector is equal to all zeros, ensured by evaluateImplicitODE

    // setup system Jacobian computation at the current time
    phiEvaluator_->setLinearizationPoint(inArgs);

    // if requested, update any hyperpameters of the phiEvaluator
    if ( (workingState->getIndex() < 2 || workingState->getIndex() % adapt_phi_evaluator_interval_ == 0)
      && (adapt_phi_evaluator_interval_ > 0))
    {
      phiEvaluator_->adaptEvaluator();
    }

    // compute the first RHS evaluation
    auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
    // compute the right hand side for at x (which is still equal to xOld)
    // and potentially set the correct Dirichlet BC to x
    RCP<Thyra::VectorBase<Scalar>> Mf = Thyra::createMember(x->space());

    RCP<Thyra::VectorBase<Scalar>> xDotOld = currentState->getXDot();
    if (this->getUseFSAL() && currentState->getIsSynced() && xDotOld != Teuchos::null) {
      // Get Mf = -M*f(xOld, t0) from xDotOld = f(xOld, t0)
      phiEvaluator_->applyMass(Mf.ptr(), xDotOld);
      Thyra::scale(Scalar(-1.0), Mf.ptr());

      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      *out << "read xDotOld at t=" << t0 << " with fsal=" << this->getUseFSAL() << std::endl;
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
        phiEvaluator_->solveMass(xDotOld.ptr(), Mf);
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
    if (temporal_finite_difference_eps_ > 0.0) {
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
      Scalar tOldOld = solutionHistory->getStateTimeIndexNM2()->getTime();
      RCP<const Thyra::VectorBase<Scalar>> xOldOld = solutionHistory->getStateTimeIndexNM2()->getX();
      RCP<Thyra::VectorBase<Scalar>> xDotOldOld = solutionHistory->getStateTimeIndexNM2()->getXDot();

      // TODO: In the future xDotOldOld may contain the residual
      // and be used to save an additional residual evaluation.
      RCP<Thyra::VectorBase<Scalar>> remf = Thyra::createMember(Mf->space());
      computeRemf(remf, xOldOld, tOldOld, xOld, t0, dt, Mf, dt_Mf_deriv);

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

    // if temporal_finite_difference_eps_ is positive or we use EPI3,
    // use the extension formula.
    if (temporal_finite_difference_eps_ > 0.0 || useEPI3) {
      sStatus = phiEvaluator_->computePhis(vphi.ptr(), dt, Mrhs_B());
    }
    else {
      sStatus = phiEvaluator_->computePhi(vphi.ptr(), 1, dt, Mf);
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
      // solve the mass matrix and scale to obtain final xDot
      phiEvaluator_->solveMass(xDot.ptr(), Mf);
      Thyra::scale(Scalar(-1.0), xDot.ptr());

      // xDot is a pointer to workingState->getXDot()
      workingState->setIsSynced(true);

      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      *out << "evaled f at t=" << time << " with fsal=" << this->getUseFSAL() << std::endl;
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


template <class Scalar>
void StepperEPI<Scalar>::computeTemporalFD(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf
)
{
  // evaluate Mf at t + dt * eps
  auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
  Teuchos::RCP<Thyra::VectorBase<Scalar>> x_BC = x->clone_v();
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xDot = this->getStepperXDot();
  this->evaluateExponentialODE(
      dt_Mf_deriv, x_BC, xDot, t0 + dt * this->temporal_finite_difference_eps_, p);
  // this sets the Dirichlet BC in x_BC to the value at t0 + dt * eps
  // thus we need to use a temporary x_BC instead of x here

  // compute dt times the temporal finite difference of Mf:
  // dt_Mf_deriv = (Mf(t + eps*dt) - Mf(t)) / eps
  Scalar one_over_eps = Scalar(1. / this->temporal_finite_difference_eps_);
  Thyra::linear_combination<Scalar>(Teuchos::tuple(-one_over_eps),
                                    Teuchos::tuple(Mf.getConst().ptr()),
                                    one_over_eps, dt_Mf_deriv.ptr());

  // also compute the finite difference (x_BC - x) / eps, reusing storage
  Teuchos::RCP<Thyra::VectorBase<Scalar>> dt_x_BC = x_BC;
  Thyra::linear_combination<Scalar>(Teuchos::tuple(-one_over_eps),
                                    Teuchos::tuple(x.getConst().ptr()),
                                    one_over_eps, dt_x_BC.ptr());

  // add dt_x_BC to dt_Mf_deriv:
  // this simple logic should take care of _either_ inhomogeneous BC or RHS.
  //
  // TODO: In case both are inhomogeneous, dt_Mf_deriv contains two contributions:
  // dt_Mf_deriv = (Mf(t + eps*dt, x + eps*dt_x_BC) - Mf(t, x)) / eps
  //             = (Mf(t + eps*dt, x) - Mf(t, x)) / eps
  //             + (Mf(t + eps*dt, x + eps*dt_x_BC) - Mf(t + eps*dt, x)) / eps
  // The second term in the sum is approximately J(t + eps*dt, x) * dt_x_BC,
  // which we may have to correct for
  Thyra::Vp_V(dt_Mf_deriv.ptr(), *dt_x_BC);
}


template <class Scalar>
void StepperEPI<Scalar>::computeRemf(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& remf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& xr,
    const Scalar tr,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x0,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mfr
)
{
  auto p = Teuchos::rcp(new ExponentialODEParameters<Scalar>(dt));
  // Evaluate the rhs -M*F at (xr, tr)
  //   Mf_old = -M*F(xr, tr) = F_impl(xDot = 0, xr, tr)
  // this will reset xDot to zero and re-set the BC in xr (that is why we need a non_const temporary)
  // if we save Mfr in the previous iteration correctly, we can reuse it.
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xr_temp = xr->clone_v();
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xDot = this->getStepperXDot();
  this->evaluateExponentialODE(remf, xr_temp, xDot, tr, p);
  // remf = Mf_old

  // update remf = (Mf - Mf_old)
  // Mf is rhs at the current time: -M*F(x0, t0)
  Thyra::linear_combination<Scalar>(Teuchos::tuple(Scalar(1.0)),
                                    Teuchos::tuple(Mf.getConst().ptr()),
                                    Scalar(-1.0), remf.ptr());

  // xd = xr - x0
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xd = xr_temp;  // reuse xr_temp memory
  Thyra::V_VpStV(xd.ptr(), *xr, Scalar(-1.0), *x0);

  // J_xd = -M*J * (xr - x0)
  Teuchos::RCP<Thyra::VectorBase<Scalar>> MJ_xd = Thyra::createMember(x0->space());
  phiEvaluator_->applyJacobian(MJ_xd.ptr(), xd);
  // TODO: could save one temp vector by extending the interface of apply Jacobian

  // update remf = (Mf - Mf_old) + J_xd
  Thyra::Vp_V(remf.ptr(), *MJ_xd);

  // add time derivative remainder term; only nonzero in nonautonomous case
  // add  -(M * F') * (tr - t0) = -(dt * M * F') * ((tr - t0) / dt)
  if (dt_Mf_deriv != Teuchos::null) {
    Thyra::Vp_StV(remf.ptr(), Scalar((tr - t0) / dt), *dt_Mf_deriv);
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
StepperEPI<Scalar>::getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperEPI<Scalar>::describe(
  Teuchos::FancyOStream               &out,
  const Teuchos::EVerbosityLevel      verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

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
  if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (stepperEPIAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The EPI AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperEPI<Scalar>::getValidParameters() const
{
  auto pl = this->getValidParametersBasicImplicit();
  return pl;
}


template <class Scalar>
void StepperEPI<Scalar>::setStepperExponentialValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<Teuchos::ParameterList> phiPL = Teuchos::null;

  temporal_finite_difference_eps_ = pl->get<double>("Epsilon for RHS finite difference", 1e-4);
  adapt_phi_evaluator_interval_ = pl->get<int>("Adapt PhiEvaluator Interval", -1);

  if (pl != Teuchos::null) {
    // TODO read in the pl for the exponential solver
    phiPL = sublist(pl, "PhiEvaluator");
  }

  // TODO: is this the right place to initialize the PhiEvaluator?
  auto phif = Teuchos::rcp(new PhiEvaluatorFactory<Scalar>());
  phiEvaluator_ = phif->createPhiEvaluator(phiPL);
}


// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template<class Scalar>
Teuchos::RCP<StepperEPI<Scalar> >
createStepperEPI(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperEPI<Scalar>());

  stepper->setStepperImplicitValues(pl);

  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


} // namespace Tempus
#endif // Tempus_StepperEPI_impl_hpp
