//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlpha_impl_hpp
#define Tempus_StepperHHTAlpha_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Tempus_StepperHHTAlphaModifierDefault.hpp"

//#define VERBOSE_DEBUG_OUTPUT
//#define DEBUG_OUTPUT

namespace Tempus {

template <class Scalar>
void StepperHHTAlpha<Scalar>::predictVelocity(
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
void StepperHHTAlpha<Scalar>::predictDisplacement(
    Thyra::VectorBase<Scalar>& dPred, const Thyra::VectorBase<Scalar>& d,
    const Thyra::VectorBase<Scalar>& v, const Thyra::VectorBase<Scalar>& a,
    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // dPred = dt*v + dt*dt/2.0*(1.0-2.0*beta_)*a
  Scalar aConst = dt * dt / 2.0 * (1.0 - 2.0 * beta_);
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), dt, v, aConst, a);
  // dPred += d;
  Thyra::Vp_V(Teuchos::ptrFromRef(dPred), d, 1.0);
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::predictVelocity_alpha_f(
    Thyra::VectorBase<Scalar>& vPred, const Thyra::VectorBase<Scalar>& v) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // vPred = (1-alpha_f)*vPred + alpha_f*v
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0 - alpha_f_, vPred, alpha_f_,
                   v);
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::predictDisplacement_alpha_f(
    Thyra::VectorBase<Scalar>& dPred, const Thyra::VectorBase<Scalar>& d) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // dPred = (1-alpha_f)*dPred + alpha_f*d
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), 1.0 - alpha_f_, dPred, alpha_f_,
                   d);
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::correctAcceleration(
    Thyra::VectorBase<Scalar>& a_n_plus1,
    const Thyra::VectorBase<Scalar>& a_n) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Scalar c = 1.0 / (1.0 - alpha_m_);
  // a_n_plus1 = 1.0/(1.0-alpha_m_)*a_n_plus1 - alpha_m/(1.0-alpha_m)*a_n =
  // (1-alpha_f)*vPred + alpha_f*v
  Thyra::V_StVpStV(Teuchos::ptrFromRef(a_n_plus1), c, a_n_plus1, -c * alpha_m_,
                   a_n);
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::correctVelocity(
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
void StepperHHTAlpha<Scalar>::correctDisplacement(
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
void StepperHHTAlpha<Scalar>::setBeta(Scalar beta)
{
  if (schemeName_ != "Newmark Beta User Defined") {
    out_->setOutputToRootOnly(0);
    *out_ << "\nWARNING: schemeName != 'Newmark Beta User Defined' (= '"
          << schemeName_ << "').\n"
          << " Leaving as beta = " << beta_ << "!\n";
    return;
  }

  beta_ = beta;

  if (beta_ == 0.0) {
    out_->setOutputToRootOnly(0);
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

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::setGamma(Scalar gamma)
{
  if (schemeName_ != "Newmark Beta User Defined") {
    out_->setOutputToRootOnly(0);
    *out_ << "\nWARNING: schemeName != 'Newmark Beta User Defined' (= '"
          << schemeName_ << "').\n"
          << " Leaving as gamma = " << gamma_ << "!\n";
    return;
  }

  gamma_ = gamma;

  TEUCHOS_TEST_FOR_EXCEPTION(
      (gamma_ > 1.0) || (gamma_ < 0.0), std::logic_error,
      "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Gamma ="
          << gamma_ << ".  Please select Gamma >= 0 and <= 1. \n");

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::setAlphaF(Scalar alpha_f)
{
  alpha_f_ = alpha_f;

  TEUCHOS_TEST_FOR_EXCEPTION(
      (alpha_f_ > 1.0) || (alpha_f_ < 0.0), std::logic_error,
      "\nError in 'HHT-Alpha' stepper: invalid value of Alpha_f = "
          << alpha_f_ << ".  Please select Alpha_f >= 0 and <= 1. \n");

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::setAlphaM(Scalar alpha_m)
{
  alpha_m_ = alpha_m;

  TEUCHOS_TEST_FOR_EXCEPTION(
      (alpha_m_ >= 1.0) || (alpha_m_ < 0.0), std::logic_error,
      "\nError in 'HHT-Alpha' stepper: invalid value of Alpha_m = "
          << alpha_m_ << ".  Please select Alpha_m >= 0 and < 1. \n");

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::setSchemeName(std::string schemeName)
{
  schemeName_ = schemeName;

  if (schemeName_ == "Newmark Beta Average Acceleration") {
    beta_  = 0.25;
    gamma_ = 0.5;
  }
  else if (schemeName_ == "Newmark Beta Linear Acceleration") {
    beta_  = 0.25;
    gamma_ = 1.0 / 6.0;
  }
  else if (schemeName_ == "Newmark Beta Central Difference") {
    beta_  = 0.0;
    gamma_ = 0.5;
  }
  else if (schemeName_ == "Newmark Beta User Defined") {
    beta_  = 0.25;
    gamma_ = 0.5;  // Use defaults until setBeta and setGamma calls.
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "\nError in Tempus::StepperHHTAlpha!  "
            << "Invalid Scheme Name = " << schemeName_ << ".  \n"
            << "Valid Scheme Names are: 'Newmark Beta Average Acceleration', "
            << "'Newmark Beta Linear Acceleration', \n"
            << "'Newmark Beta Central Difference' and 'Newmark Beta User "
               "Defined'.\n");
  }

  this->isInitialized_ = false;
}

template <class Scalar>
StepperHHTAlpha<Scalar>::StepperHHTAlpha()
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  this->setStepperName("HHT-Alpha");
  this->setStepperType("HHT-Alpha");
  this->setUseFSAL(false);
  this->setICConsistency("None");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);
  this->setSchemeName("Newmark Beta Average Acceleration");
  this->setAlphaF(0.0);
  this->setAlphaM(0.0);
  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}

template <class Scalar>
StepperHHTAlpha<Scalar>::StepperHHTAlpha(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    bool zeroInitialGuess, std::string schemeName, Scalar beta, Scalar gamma,
    Scalar alpha_f, Scalar alpha_m,
    const Teuchos::RCP<StepperHHTAlphaAppAction<Scalar> >&
        stepperHHTAlphaAppAction)
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  this->setStepperName("HHT-Alpha");
  this->setStepperType("HHT-Alpha");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);
  this->setSchemeName(schemeName);
  if (schemeName == "Newmark Beta User Defined") {
    this->setBeta(beta);
    this->setGamma(gamma);
  }
  this->setAlphaF(alpha_f);
  this->setAlphaM(alpha_m);
  this->setAppAction(stepperHHTAlphaAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::setAppAction(
    Teuchos::RCP<StepperHHTAlphaAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperHHTAlphaAppAction_ =
        Teuchos::rcp(new StepperHHTAlphaModifierDefault<Scalar>());
  }
  else {
    stepperHHTAlphaAppAction_ = appAction;
  }
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  validSecondOrderODE_DAE(appModel);
  Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
      Teuchos::rcp(
          new WrapperModelEvaluatorSecondOrder<Scalar>(appModel, "HHT-Alpha"));
  this->wrapperModel_ = wrapperModel;

  TEUCHOS_TEST_FOR_EXCEPTION(this->solver_ == Teuchos::null, std::logic_error,
                             "Error - Solver is not set!\n");
  if (this->wrapperModel_ != Teuchos::null)
    this->solver_->setModel(this->wrapperModel_);

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperHHTAlpha::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperHHTAlpha<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for HHTAlpha.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\nTry setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");

    RCP<StepperHHTAlpha<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperHHTAlphaAppAction_->execute(
        solutionHistory, thisStepper,
        StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState =
        solutionHistory->getCurrentState();

    Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
        Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
            this->wrapperModel_);

    // Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar> > d_old = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > v_old = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_old       = currentState->getXDotDot();

#ifdef DEBUG_OUTPUT
    // IKT, 3/21/17, debug output: pring d_old, v_old, a_old to check for
    //  correctness.
    *out_ << "IKT d_old = " << Thyra::max(*d_old) << "\n";
    *out_ << "IKT v_old = " << Thyra::max(*v_old) << "\n";
#endif

    // Get new values of d, v and a from current workingState
    //(to be updated here)
    RCP<Thyra::VectorBase<Scalar> > d_new = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > v_new = workingState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_new = workingState->getXDotDot();

    // Get time and dt
    const Scalar time = currentState->getTime();
    const Scalar dt   = workingState->getTimeStep();
    // Update time
    Scalar t = time + dt;

    // Compute initial acceleration, a_old, using initial displacement (d_old)
    // and initial velocity (v_old) if in 1st time step
    if (time == solutionHistory->minTime()) {
      RCP<Thyra::VectorBase<Scalar> > d_init =
          Thyra::createMember(d_old->space());
      RCP<Thyra::VectorBase<Scalar> > v_init =
          Thyra::createMember(v_old->space());
      RCP<Thyra::VectorBase<Scalar> > a_init =
          Thyra::createMember(a_old->space());
      Thyra::copy(*d_old, d_init.ptr());
      Thyra::copy(*v_old, v_init.ptr());
      if (this->initialGuess_ !=
          Teuchos::null) {  // set initial guess for Newton, if provided
        // Throw an exception if initial_guess is not compatible with solution
        bool is_compatible =
            (a_init->space())->isCompatible(*this->initialGuess_->space());
        TEUCHOS_TEST_FOR_EXCEPTION(
            is_compatible != true, std::logic_error,
            "Error in Tempus::NemwarkImplicitAForm takeStep(): user-provided "
            "initial guess'!\n"
                << "for Newton is not compatible with solution vector!\n");
        Thyra::copy(*this->initialGuess_, a_init.ptr());
      }
      else {  // if no initialGuess_ provide, set 0 initial guess
        Thyra::put_scalar(0.0, a_init.ptr());
      }
      wrapperModel->initializeNewmark(v_init, d_init, 0.0, time, beta_, gamma_);
      const Thyra::SolveStatus<Scalar> sStatus =
          (*(this->solver_)).solve(&*a_init);

      workingState->setSolutionStatus(sStatus);  // Converged --> pass.
      Thyra::copy(*a_init, a_old.ptr());
    }
#ifdef DEBUG_OUTPUT
    // IKT, 3/30/17, debug output: pring a_old to check for correctness.
    *out_ << "IKT a_old = " << Thyra::max(*a_old) << "\n";
#endif

    // allocate d and v predictors
    RCP<Thyra::VectorBase<Scalar> > d_pred =
        Thyra::createMember(d_old->space());
    RCP<Thyra::VectorBase<Scalar> > v_pred =
        Thyra::createMember(v_old->space());

    // compute displacement and velocity predictors
    predictDisplacement(*d_pred, *d_old, *v_old, *a_old, dt);
    predictVelocity(*v_pred, *v_old, *a_old, dt);

    // compute second displacement and velocity predictors (those that are
    // functions of alpha_f)
    predictDisplacement_alpha_f(*d_pred, *d_old);
    predictVelocity_alpha_f(*v_pred, *v_old);

    // inject d_pred, v_pred, a and other relevant data into wrapperModel
    wrapperModel->initializeNewmark(v_pred, d_pred, dt, t, beta_, gamma_);

    stepperHHTAlphaAppAction_->execute(
        solutionHistory, thisStepper,
        StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

    // Solve for new acceleration
    const Thyra::SolveStatus<Scalar> sStatus =
        (*(this->solver_)).solve(&*a_new);

    stepperHHTAlphaAppAction_->execute(
        solutionHistory, thisStepper,
        StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);

    // correct acceleration (function of alpha_m)
    correctAcceleration(*a_new, *a_old);

    // correct velocity and displacement
    correctVelocity(*v_new, *v_pred, *a_new, dt);
    correctDisplacement(*d_new, *d_pred, *a_new, dt);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);

    stepperHHTAlphaAppAction_->execute(
        solutionHistory, thisStepper,
        StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
StepperHHTAlpha<Scalar>::getDefaultStepperState()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperHHTAlpha<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  *l_out << std::endl;
  Stepper<Scalar>::describe(*l_out, verbLevel);
  StepperImplicit<Scalar>::describe(*l_out, verbLevel);

  *l_out << "--- StepperHHTAlpha ---\n";
  *l_out << "  schemeName_ = " << schemeName_ << std::endl;
  *l_out << "  beta_       = " << beta_ << std::endl;
  *l_out << "  gamma_      = " << gamma_ << std::endl;
  *l_out << "  alpha_f_    = " << alpha_f_ << std::endl;
  *l_out << "  alpha_m_    = " << alpha_m_ << std::endl;
  *l_out << "-----------------------" << std::endl;
}

template <class Scalar>
bool StepperHHTAlpha<Scalar>::isValidSetup(Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

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

  return isValidSetup;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperHHTAlpha<Scalar>::getValidParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  auto pl = this->getValidParametersBasicImplicit();

  auto hhtalphaPL = Teuchos::parameterList("HHT-Alpha Parameters");
  hhtalphaPL->set<std::string>("Scheme Name", schemeName_);
  hhtalphaPL->set<double>("Beta", beta_);
  hhtalphaPL->set<double>("Gamma", gamma_);
  hhtalphaPL->set<double>("Alpha_f", alpha_f_);
  hhtalphaPL->set<double>("Alpha_m", alpha_m_);
  pl->set("HHT-Alpha Parameters", *hhtalphaPL);

  return pl;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperHHTAlpha<Scalar> > createStepperHHTAlpha(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperHHTAlpha<Scalar>());
  stepper->setStepperImplicitValues(pl);

  if (pl != Teuchos::null) {
    if (pl->isSublist("HHT-Alpha Parameters")) {
      auto hhtalphaPL        = pl->sublist("HHT-Alpha Parameters", true);
      std::string schemeName = hhtalphaPL.get<std::string>(
          "Scheme Name", "Newmark Beta Average Acceleration");
      stepper->setSchemeName(schemeName);
      if (schemeName == "Newmark Beta User Defined") {
        stepper->setBeta(hhtalphaPL.get<double>("Beta", 0.25));
        stepper->setGamma(hhtalphaPL.get<double>("Gamma", 0.5));
      }
      stepper->setAlphaF(hhtalphaPL.get<double>("Alpha_f", 0.0));
      stepper->setAlphaM(hhtalphaPL.get<double>("Alpha_m", 0.0));
    }
    else {
      stepper->setSchemeName("Newmark Beta Average Acceleration");
      stepper->setAlphaF(0.0);
      stepper->setAlphaM(0.0);
    }
  }

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperHHTAlpha_impl_hpp
