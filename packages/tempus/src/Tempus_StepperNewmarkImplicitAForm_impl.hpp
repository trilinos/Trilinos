//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitAForm_impl_hpp
#define Tempus_StepperNewmarkImplicitAForm_impl_hpp

#include "Tempus_StepperNewmarkImplicitAFormModifierDefault.hpp"

//#define VERBOSE_DEBUG_OUTPUT
//#define DEBUG_OUTPUT

namespace Tempus {

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::predictVelocity(
    Thyra::VectorBase<Scalar>& vPred, const Thyra::VectorBase<Scalar>& v,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // vPred = v + dt*(1.0-gamma_)*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0, v, dt * (1.0 - gamma_), a);
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::predictDisplacement(
    Thyra::VectorBase<Scalar>& dPred, const Thyra::VectorBase<Scalar>& d,
    const Thyra::VectorBase<Scalar>& v, const Thyra::VectorBase<Scalar>& a,
    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > tmp =
      Thyra::createMember<Scalar>(dPred.space());
  // dPred = dt*v + dt*dt/2.0*(1.0-2.0*beta_)*a
  Scalar aConst = dt * dt / 2.0 * (1.0 - 2.0 * beta_);
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), dt, v, aConst, a);
  // dPred += d;
  Thyra::Vp_V(Teuchos::ptrFromRef(dPred), d, 1.0);
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::correctVelocity(
    Thyra::VectorBase<Scalar>& v, const Thyra::VectorBase<Scalar>& vPred,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // v = vPred + dt*gamma_*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(v), 1.0, vPred, dt * gamma_, a);
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::correctDisplacement(
    Thyra::VectorBase<Scalar>& d, const Thyra::VectorBase<Scalar>& dPred,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // d = dPred + beta_*dt*dt*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(d), 1.0, dPred, beta_ * dt * dt, a);
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setBeta(Scalar beta)
{
  if (schemeName_ != "User Defined") {
    *out_ << "\nWARNING: schemeName != 'User Defined' (=" << schemeName_
          << ").\n"
          << " Not setting beta, and leaving as beta = " << beta_ << "!\n";
    return;
  }

  beta_ = beta;

  if (beta_ == 0.0) {
    *out_ << "\nWARNING: Running (implicit implementation of) Newmark "
          << "Implicit a-Form Stepper with Beta = 0.0, which \n"
          << "specifies an explicit scheme.  Mass lumping is not possible, "
          << "so this will be slow!  To run explicit \n"
          << "implementation of Newmark Implicit a-Form Stepper, please "
          << "re-run with 'Stepper Type' = 'Newmark Explicit a-Form'.\n"
          << "This stepper allows for mass lumping when called through "
          << "Piro::TempusSolver.\n";
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      (beta_ > 1.0) || (beta_ < 0.0), std::logic_error,
      "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Beta = "
          << beta_ << ".  Please select Beta >= 0 and <= 1. \n");
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setGamma(Scalar gamma)
{
  if (schemeName_ != "User Defined") {
    *out_ << "\nWARNING: schemeName != 'User Defined' (=" << schemeName_
          << ").\n"
          << " Not setting gamma, and leaving as gamma = " << gamma_ << "!\n";
    return;
  }

  gamma_ = gamma;

  TEUCHOS_TEST_FOR_EXCEPTION(
      (gamma_ > 1.0) || (gamma_ < 0.0), std::logic_error,
      "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Gamma ="
          << gamma_ << ".  Please select Gamma >= 0 and <= 1. \n");
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setSchemeName(std::string schemeName)
{
  schemeName_ = schemeName;

  if (schemeName_ == "Average Acceleration") {
    beta_  = 0.25;
    gamma_ = 0.5;
  }
  else if (schemeName_ == "Linear Acceleration") {
    beta_  = 0.25;
    gamma_ = 1.0 / 6.0;
  }
  else if (schemeName_ == "Central Difference") {
    beta_  = 0.0;
    gamma_ = 0.5;
  }
  else if (schemeName_ == "User Defined") {
    beta_  = 0.25;
    gamma_ = 0.5;  // Use defaults until setBeta and setGamma calls.
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "\nError in Tempus::StepperNewmarkImplicitAForm!  "
            << "Invalid Scheme Name = " << schemeName_ << ".  \n"
            << "Valid Scheme Names are: 'Average Acceleration', "
            << "'Linear Acceleration', \n"
            << "'Central Difference' and 'User Defined'.\n");
  }

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setAppAction(
    Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperNewmarkImpAppAction_ =
        Teuchos::rcp(new StepperNewmarkImplicitAFormModifierDefault<Scalar>());
  }
  else {
    stepperNewmarkImpAppAction_ = appAction;
  }

  this->isInitialized_ = false;
}

template <class Scalar>
StepperNewmarkImplicitAForm<Scalar>::StepperNewmarkImplicitAForm()
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  this->setStepperName("Newmark Implicit a-Form");
  this->setStepperType("Newmark Implicit a-Form");
  this->setUseFSAL(true);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);
  this->setSchemeName("Average Acceleration");

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}

template <class Scalar>
StepperNewmarkImplicitAForm<Scalar>::StepperNewmarkImplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    bool zeroInitialGuess, std::string schemeName, Scalar beta, Scalar gamma,
    const Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> >&
        stepperAppAction)
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  this->setStepperName("Newmark Implicit a-Form");
  this->setStepperType("Newmark Implicit a-Form");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);
  this->setSchemeName(schemeName);
  this->setBeta(beta);
  this->setGamma(gamma);
  this->setAppAction(stepperAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  validSecondOrderODE_DAE(appModel);
  this->wrapperModel_ =
      Teuchos::rcp(new WrapperModelEvaluatorSecondOrder<Scalar>(
          appModel, "Newmark Implicit a-Form"));

  TEUCHOS_TEST_FOR_EXCEPTION(this->getSolver() == Teuchos::null,
                             std::logic_error, "Error - Solver is not set!\n");
  this->getSolver()->setModel(this->wrapperModel_);

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setInitialConditions(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(
      numStates < 1, std::logic_error,
      "Error - setInitialConditions() needs at least one SolutionState\n"
      "        to set the initial condition.  Number of States = "
          << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1,
                         "StepperNewmarkImplicitAForm::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"
         << std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x        = initialState->getX();
  RCP<Thyra::VectorBase<Scalar> > xDot     = initialState->getXDot();

  auto inArgs = this->wrapperModel_->getNominalValues();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !((x != Teuchos::null && xDot != Teuchos::null) ||
        (inArgs.get_x() != Teuchos::null &&
         inArgs.get_x_dot() != Teuchos::null)),
      std::logic_error,
      "Error - We need to set the initial conditions for x and xDot from\n"
      "        either initialState or appModel_->getNominalValues::InArgs\n"
      "        (but not from a mixture of the two).\n");

  // Use x and xDot from inArgs as ICs, if needed.
  if (x == Teuchos::null || xDot == Teuchos::null) {
    using Teuchos::rcp_const_cast;
    TEUCHOS_TEST_FOR_EXCEPTION(
        (inArgs.get_x() == Teuchos::null) ||
            (inArgs.get_x_dot() == Teuchos::null),
        std::logic_error,
        "Error - setInitialConditions() needs the ICs from the initialState\n"
        "        or getNominalValues()!\n");
    x = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
    initialState->setX(x);
    xDot = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot());
    initialState->setXDot(xDot);
  }

  // Check if we need Stepper storage for xDotDot
  if (initialState->getXDotDot() == Teuchos::null)
    initialState->setXDotDot(initialState->getX()->clone_v());
  else
    this->setStepperXDotDot(initialState->getXDotDot());

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (initialState->getXDotDot() == Teuchos::null) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(
          out, 1, "StepperNewmarkImplicitAForm::setInitialConditions()");
      *out << "Warning -- Requested IC consistency of 'None' but\n"
           << "           initialState does not have an xDot.\n"
           << "           Setting a 'Zero' xDot!\n"
           << std::endl;

      Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
    }
  }
  else if (icConsistency == "Zero")
    Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
  else if (icConsistency == "App") {
    auto xDotDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
        inArgs.get_x_dot_dot());
    TEUCHOS_TEST_FOR_EXCEPTION(
        xDotDot == Teuchos::null, std::logic_error,
        "Error - setInitialConditions() requested 'App' for IC consistency,\n"
        "        but 'App' returned a null pointer for xDotDot!\n");
    Thyra::assign(this->getStepperXDotDot(initialState).ptr(), *xDotDot);
  }
  else if (icConsistency == "Consistent") {
    // Solve f(x, xDot, xDotDot, t) = 0.
    const Scalar time = initialState->getTime();
    auto xDotDot      = this->getStepperXDotDot(initialState);

    // Compute initial acceleration using initial displacement
    // and initial velocity.
    if (this->initialGuess_ != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          !((xDotDot->space())->isCompatible(*this->initialGuess_->space())),
          std::logic_error,
          "Error - User-provided initial guess for Newton is not compatible\n"
          "        with solution vector!\n");
      Thyra::copy(*this->initialGuess_, xDotDot.ptr());
    }
    else {
      Thyra::put_scalar(0.0, xDotDot.ptr());
    }

    auto wrapperModel =
        Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
            this->wrapperModel_);

    wrapperModel->initializeNewmark(xDot, x, 0.0, time, beta_, gamma_);
    const Thyra::SolveStatus<Scalar> sStatus =
        (*(this->solver_)).solve(&*xDotDot);

    TEUCHOS_TEST_FOR_EXCEPTION(
        sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED, std::logic_error,
        "Error - Solver failed while determining the initial conditions.\n"
        "        Solver status is "
            << Thyra::toString(sStatus.solveStatus) << ".\n");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error - setInitialConditions() invalid IC consistency, "
            << icConsistency << ".\n");
  }

  // At this point, x, xDot and xDotDot are sync'ed or consistent
  // at the same time level for the initialState.
  initialState->setIsSynced(true);

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    auto f       = initialState->getX()->clone_v();
    auto xDotDot = this->getStepperXDotDot(initialState);

    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Scalar> appInArgs =
        this->wrapperModel_->getAppModel()->createInArgs();
    MEB::OutArgs<Scalar> appOutArgs =
        this->wrapperModel_->getAppModel()->createOutArgs();

    appInArgs.set_x(x);
    appInArgs.set_x_dot(xDot);
    appInArgs.set_x_dot_dot(xDotDot);

    appOutArgs.set_f(appOutArgs.get_f());

    appInArgs.set_W_x_dot_dot_coeff(Scalar(0.0));  // da/da
    appInArgs.set_alpha(Scalar(0.0));              // dv/da
    appInArgs.set_beta(Scalar(0.0));               // dd/da

    appInArgs.set_t(initialState->getTime());

    this->wrapperModel_->getAppModel()->evalModel(appInArgs, appOutArgs);

    Scalar reldiff = Thyra::norm(*f);
    Scalar normx   = Thyra::norm(*x);
    Scalar eps     = Scalar(100.0) * std::abs(Teuchos::ScalarTraits<Scalar>::eps());
    if (normx > eps * reldiff) reldiff /= normx;

    if (reldiff > eps) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(
          out, 1, "StepperNewmarkImplicitAForm::setInitialConditions()");
      *out << "Warning -- Failed consistency check but continuing!\n"
           << "  ||f(x,xDot,xDotDot,t)||/||x|| > eps" << std::endl
           << "  ||f(x,xDot,xDotDot,t)||       = " << Thyra::norm(*f)
           << std::endl
           << "  ||x||                         = " << Thyra::norm(*x)
           << std::endl
           << "  ||f(x,xDot,xDotDot,t)||/||x|| = " << reldiff << std::endl
           << "                            eps = " << eps << std::endl;
    }
  }

  if (!(this->getUseFSAL())) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1,
                         "StepperNewmarkImplicitAForm::setInitialConditions()");
    *out << "\nWarning -- The First-Same-As-Last (FSAL) principle is "
         << "part of the Newmark Implicit A-Form.  The default is to "
         << "set useFSAL=true, and useFSAL=false will be ignored." << std::endl;
  }
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperNewmarkImplicitAForm::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperNewmarkImplicitAForm<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for NewmarkImplicitAForm.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\nTry setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n"
            << "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");

    auto thisStepper = Teuchos::rcpFromRef(*this);
    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitAFormAppAction<
            Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState =
        solutionHistory->getCurrentState();

    auto wrapperModel =
        Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
            this->wrapperModel_);

    // Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar> > d_old = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > v_old = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_old       = currentState->getXDotDot();

    // Get new values of d, v and a from workingState
    RCP<Thyra::VectorBase<Scalar> > d_new = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > v_new = workingState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_new = workingState->getXDotDot();

    // Get time and dt
    const Scalar time = currentState->getTime();
    const Scalar dt   = workingState->getTimeStep();
    Scalar t          = time + dt;

    // Compute displacement and velocity predictors
    predictDisplacement(*d_new, *d_old, *v_old, *a_old, dt);
    predictVelocity(*v_new, *v_old, *a_old, dt);

    // Inject d_new, v_new, a and other relevant data into wrapperModel
    wrapperModel->initializeNewmark(v_new, d_new, dt, t, beta_, gamma_);

    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitAFormAppAction<
            Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

    if (this->getZeroInitialGuess())
      Thyra::assign(a_new.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

    // Solve nonlinear system with a_new as initial guess
    const Thyra::SolveStatus<Scalar> sStatus =
        (*(this->solver_)).solve(&*a_new);

    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitAFormAppAction<
            Scalar>::ACTION_LOCATION::AFTER_SOLVE);

    // Correct velocity, displacement.
    correctVelocity(*v_new, *v_new, *a_new, dt);
    correctDisplacement(*d_new, *d_new, *a_new, dt);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);

    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitAFormAppAction<
            Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}

/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template <class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperNewmarkImplicitAForm<Scalar>::getDefaultStepperState()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  out.setOutputToRootOnly(0);
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperNewmarkImplicitAForm ---\n";
  out << "  schemeName_ = " << schemeName_ << std::endl;
  out << "  beta_       = " << beta_ << std::endl;
  out << "  gamma_      = " << gamma_ << std::endl;
  out << "-----------------------------------" << std::endl;
}

template <class Scalar>
bool StepperNewmarkImplicitAForm<Scalar>::isValidSetup(
    Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;
  out.setOutputToRootOnly(0);

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;

  // if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if (this->wrapperModel_->getAppModel() == Teuchos::null) {
    isValidSetup = false;
    out << "The application ModelEvaluator is not set!\n";
  }

  if (this->wrapperModel_ == Teuchos::null) {
    isValidSetup = false;
    out << "The wrapper ModelEvaluator is not set!\n";
  }

  if (this->solver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The solver is not set!\n";
  }

  if (this->stepperNewmarkImpAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Newmark Implicit A-Form AppAction is not set!\n";
  }

  return isValidSetup;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperNewmarkImplicitAForm<Scalar>::getValidParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  auto pl = this->getValidParametersBasicImplicit();

  auto newmarkPL = Teuchos::parameterList("Newmark Parameters");
  newmarkPL->set<std::string>("Scheme Name", schemeName_);
  newmarkPL->set<double>("Beta", beta_);
  newmarkPL->set<double>("Gamma", gamma_);
  pl->set("Newmark Parameters", *newmarkPL);

  return pl;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> >
createStepperNewmarkImplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperNewmarkImplicitAForm<Scalar>());
  stepper->setStepperImplicitValues(pl);

  if (pl != Teuchos::null) {
    if (pl->isSublist("Newmark Parameters")) {
      auto newmarkPL = pl->sublist("Newmark Parameters", true);
      std::string schemeName =
          newmarkPL.get<std::string>("Scheme Name", "Average Acceleration");
      stepper->setSchemeName(schemeName);
      if (schemeName == "User Defined") {
        stepper->setBeta(newmarkPL.get<double>("Beta", 0.25));
        stepper->setGamma(newmarkPL.get<double>("Gamma", 0.5));
      }
    }
    else {
      stepper->setSchemeName("Average Acceleration");
    }
  }

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperNewmarkImplicitAForm_impl_hpp
